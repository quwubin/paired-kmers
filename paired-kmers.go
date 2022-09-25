package main

/*
Copyright © 2022 Wubin Qu <quwubin@gmail.com>
*/

import (
	"bufio"
	"bytes"
	"context"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"sort"
	"strings"
	"time"

	"github.com/RoaringBitmap/roaring/roaring64"
	"github.com/neilotoole/errgroup"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/kmers"
	flag "github.com/spf13/pflag"
)

func checkErr(err error) {
	if err != nil {
		log.Fatal(err)
		os.Exit(-1)
	}
}

// checkMaxRuns return run number, e.g., return 5 if agcTTTTTagc
func checkMaxRuns(input interface{}) int {
	var seq []byte
	var ok bool
	if seq, ok = input.([]byte); !ok {
		if input_string, ok := input.(string); ok {
			seq = []byte(input_string)
		} else {
			panic("require string or []byte for max run checking")
		}
	}

	if len(seq) == 0 {
		return 0
	}

	if len(seq) == 1 {
		return 1
	}

	last := seq[0]
	runCount := 1
	maxRun := 1
	for i := 1; i < len(seq); i++ {
		if seq[i] == last {
			runCount++
		} else {
			runCount = 0
		}

		if maxRun < runCount {
			maxRun = runCount
		}

		// fmt.Printf("i: %d, char: %s, runCount: %d, maxRun: %d\n", i+1, string(seq[i]), runCount, maxRun)

		last = seq[i]
	}

	return maxRun + 1
}

// https://biopython.org/DIST/docs/api/Bio.SeqUtils-pysrc.html#GC
func calGC(input interface{}) float64 {
	var seq []byte
	var ok bool
	if seq, ok = input.([]byte); !ok {
		if input_string, ok := input.(string); ok {
			seq = []byte(input_string)
		} else {
			panic("require string or []byte for max run checking")
		}
	}

	if len(seq) == 0 {
		return 0
	}

	gc := 0
	for i := 0; i < len(seq); i++ {
		switch seq[i] {
		case 'C', 'c':
			gc++
		case 'G', 'g':
			gc++
		case 'S', 's':
			gc++
		}
	}

	return float64(gc*100.0) / float64(len(seq))
}

// kmer module
func buildKmer(seq []byte, k int, p Para) *roaring64.Bitmap {
	kmersSet := roaring64.New()
	nKmers := len(seq) - k + 1

	if nKmers <= 0 {
		return kmersSet
	}

	for i := 0; i < nKmers; i++ {
		kmer := seq[i : i+k]
		if p.MaxRun < k {
			maxRun := checkMaxRuns(kmer)
			if maxRun > p.MaxRun {
				continue
			}
		}

		if p.MinGC > 0 && p.MaxGC < 100 {
			gc := calGC(kmer)
			if gc < p.MinGC || gc > p.MaxGC {
				continue
			}
		}

		code, err := kmers.Encode(kmer)
		if err != nil {
			log.Printf("kmer to code err: %s", err)
			continue
		}

		kmersSet.Add(code)
	}

	return kmersSet
}

func JaccardContainment(a, b *roaring64.Bitmap) uint64 {
	if a == nil || b == nil {
		return 0
	}

	return a.AndCardinality(b)
}

type recordInData struct {
	Record *fastx.Record
	Index  uint64
}

type recordOutData struct {
	KmerSet roaring64.Bitmap
	Record  *fastx.Record
	Index   uint64
}

// kmer module
func buildKmerInfo(db string, kvalue int, cpu int, p Para) (KmerInfoList, *roaring64.Bitmap, map[uint64]*fastx.Record) {
	g, ctx := errgroup.WithContext(context.Background())
	inChan := make(chan recordInData)

	g.Go(
		func() error {
			defer close(inChan)

			reader, err := fastx.NewDefaultReader(db)
			if err != nil {
				return err
			}

			var record *fastx.Record
			index := uint64(0)
			for {
				record, err = reader.Read()
				if err != nil {
					if err == io.EOF {
						break
					}
					log.Println(err)
					break
				}

				fmt.Printf("record: %d\n", index)

				select {
				case inChan <- recordInData{Record: record.Clone(), Index: index}:
				case <-ctx.Done():
					return ctx.Err()
				}

				index++
			}

			return nil
		},
	)

	outChan := make(chan recordOutData)

	for i := 0; i < cpu; i++ {
		g.Go(
			func() error {
				for t := range inChan {
					t := t

					seq := bytes.ToUpper(t.Record.Seq.Seq)
					t.Record.Seq.Seq = seq

					kmerSet := buildKmer(seq, kvalue, p)

					// fmt.Printf("record: %d, kmer: %d\n", t.Index, kmerSet.GetCardinality())

					select {
					case outChan <- recordOutData{KmerSet: *kmerSet, Record: t.Record, Index: t.Index}:
					case <-ctx.Done():
						return ctx.Err()
					}
				}

				return nil
			},
		)
	}

	go func() {
		g.Wait()
		close(outChan)
	}()

	allRecords := roaring64.New()
	fastaHash := make(map[uint64]*fastx.Record)

	kmerInfo := make(map[uint64]*roaring64.Bitmap)

	var batch []recordOutData
	expKmers := roaring64.New()
	first := true

	if p.RecordsSample < 0 {
		first = false
	}

	for o := range outChan {
		allRecords.Add(o.Index)
		fastaHash[o.Index] = o.Record

		// fmt.Printf("first: %t, batch: %d, expKmers: %d\n", first, len(batch), expKmers.GetCardinality())

		if first {
			if len(batch) >= p.RecordsSample {
				if p.Test {
					// batch1 := clone.Slowly(batch).([]recordOutData)
					// batch2 := clone.Slowly(batch).([]recordOutData)
					unionKmers := reduceKmerByBatch(batch, false)
					unionCount := unionKmers.GetCardinality()
					interKmers := reduceKmerByBatch(batch, true)
					interCount := interKmers.GetCardinality()

					log.Printf("union kmers: %d, intersection kmers: %d, reduced %d for the first %d records", unionCount, interCount, unionCount-interCount, p.RecordsSample)
					os.Exit(0)
				}

				expKmers = reduceKmerByBatch(batch, p.Intersect)
				if p.Intersect {
					batchSet := roaring64.New()
					for _, b := range batch {
						batchSet.Add(b.Index)
					}
					i := expKmers.Iterator()
					for i.HasNext() {
						code := i.Next()
						_, ok := kmerInfo[code]
						if !ok {
							kmerInfo[code] = batchSet
						} else {
							kmerInfo[code].Or(batchSet)
						}
					}
				} else {
					for _, b := range batch {
						i := b.KmerSet.Iterator()
						for i.HasNext() {
							code := i.Next()
							_, ok := kmerInfo[code]
							if !ok {
								kmerInfo[code] = roaring64.New()
							}

							kmerInfo[code].Add(b.Index)
						}
					}
				}

				// fmt.Printf("batch: %d, expKmers: %d\n", len(batch), expKmers.GetCardinality())

				batch = nil

				first = false

				// this record shoud be processed, so NO continue here
			} else {
				batch = append(batch, o)
				continue
			}
		}

		if !expKmers.IsEmpty() {
			o.KmerSet.And(expKmers)
		}

		i := o.KmerSet.Iterator()
		for i.HasNext() {
			code := i.Next()
			_, ok := kmerInfo[code]
			if !ok {
				kmerInfo[code] = roaring64.New()
			}
			kmerInfo[code].Add(o.Index)
		}

		// fmt.Printf("o.Index: %d, o.KmerSet: %d\n", o.Index, o.KmerSet.GetCardinality())
	}

	if len(batch) > 0 {
		batchSet := roaring64.New()
		for _, b := range batch {
			batchSet.Add(b.Index)
		}

		expKmers = reduceKmerByBatch(batch, false) // no need to union even --intersection == false
		i := expKmers.Iterator()
		for i.HasNext() {
			code := i.Next()
			_, ok := kmerInfo[code]
			if !ok {
				kmerInfo[code] = batchSet
			} else {
				kmerInfo[code].Or(batchSet)
			}
		}
	}

	kiList := make(KmerInfoList, len(kmerInfo))

	if err := g.Wait(); err != nil {
		log.Println(err)
		return kiList, allRecords, fastaHash
	}

	i := 0
	for k, v := range kmerInfo {
		kiList[i] = KmerInfo{K: k, RecordSet: v}
		i++
	}

	sort.Slice(kiList, func(i, j int) bool {
		return kiList[i].RecordSet.GetCardinality() > kiList[j].RecordSet.GetCardinality()
	})

	for k := range kmerInfo {
		delete(kmerInfo, k)
	}

	return kiList, allRecords, fastaHash
}

func reduceKmerByBatch(list []recordOutData, intersect bool) *roaring64.Bitmap {
	reducedKmer := roaring64.New()
	if len(list) == 0 {
		return reducedKmer
	}

	if &list[0].KmerSet == nil {
		return reducedKmer
	}

	reducedKmer.Or(&list[0].KmerSet)

	if len(list) == 1 {
		return reducedKmer
	}

	for _, l := range list[1:] {
		if &l.KmerSet == nil {
			continue
		}

		if intersect {
			reducedKmer.And(&l.KmerSet)
		} else {
			reducedKmer.Or(&l.KmerSet)
		}

		// fmt.Printf("i: %d, mergedKmer: %d\n", i, reducedKmer.GetCardinality())
	}

	return reducedKmer
}

func fincCandiRoads(fastaHash map[uint64]*fastx.Record, k1 KmerInfo, k2 KmerInfo, para Para) PairedKmer {
	k1seq := string(kmers.Decode(k1.K, para.Kvalue))
	k2seq := string(kmers.Decode(k2.K, para.Kvalue))
	re1 := regexp.MustCompile(k1seq)
	re2 := regexp.MustCompile(k2seq)
	var candi PairedKmer

	if k1.RecordSet == nil || k2.RecordSet == nil {
		return candi
	}

	k1.RecordSet.And(k2.RecordSet)

	sharedRecords := roaring64.New()
	k1Iter := k1.RecordSet.Iterator()
	for k1Iter.HasNext() {
		recordIndex := k1Iter.Next()
		recordSeq := string(fastaHash[recordIndex].Seq.Seq)
		k1PosList := re1.FindAllStringIndex(recordSeq, -1)
		if k1PosList == nil {
			continue
		}

		k2PosList := re2.FindAllStringIndex(recordSeq, -1)
		if k2PosList == nil {
			continue
		}

		shared := false
		for _, p1 := range k1PosList {
			if shared {
				break
			}

			for _, p2 := range k2PosList {
				var f, r []int
				if p1[0] < p2[0] {
					f = p1
					r = p2
				} else {
					f = p2
					r = p1
				}

				dist := r[0] - f[0] - int(para.Kvalue)
				if dist < 0 {
					continue
				}

				ampSize := dist + 2*int(para.Kvalue)
				if ampSize < int(para.MinSize) || ampSize > int(para.MaxSize) {
					continue
				}

				shared = true

				if candi.Empty() {
					k2seq, err := seq.NewSeq(seq.DNA, []byte(k2seq))
					if err != nil {
						log.Println(err)
						continue
					}

					candi = PairedKmer{
						Index: recordIndex,
						K1:    k1seq,
						K2:    string(k2seq.RevCom().Seq),
						F5:    f[0],     // 包含
						F3:    f[1] - 1, // 包含
						R5:    r[1],     // 不包含
						R3:    r[0],     // 包含
					}

					candi.Amp5 = recordSeq[candi.F5:candi.R5]
					candi.Amp3 = recordSeq[candi.F3+1 : candi.R3]
				}

				if shared {
					break
				}
			}
		}

		if shared {
			sharedRecords.Add(recordIndex)
		}
	}

	candi.Records = sharedRecords

	return candi
}

type PairedKmer struct {
	Index   uint64 // record index
	K1      string
	K2      string
	Dist    int
	F5      int
	F3      int
	R5      int
	R3      int
	Amp5    string
	Amp3    string
	Records *roaring64.Bitmap
}

func (a *PairedKmer) Empty() bool {
	if a.K1 == "" {
		return true
	}

	return false
}

type PairedKmers []PairedKmer

type kiData struct {
	K1 KmerInfo
	KI KmerInfoList
}

func findPairRoads(kiList KmerInfoList, fastaHash map[uint64]*fastx.Record, para Para) PairedKmers {
	g, ctx := errgroup.WithContext(context.Background())
	inChan := make(chan kiData)
	var roadList PairedKmers

	g.Go(
		func() error {
			defer close(inChan)

			for i := 0; i < len(kiList); i++ {
				if i == len(kiList)-1 {
					break
				}

				k1 := kiList[i]

				select {
				case inChan <- kiData{K1: k1, KI: kiList[i+1:]}:
				case <-ctx.Done():
					return ctx.Err()
				}
			}

			return nil
		},
	)

	outChan := make(chan PairedKmer)

	for i := 0; i < para.CPU; i++ {
		g.Go(
			func() error {
				for t := range inChan {
					t := t

					k2 := findMaxContainment(t.K1, t.KI)
					if k2.Empty() {
						continue
					}

					candi := fincCandiRoads(fastaHash, t.K1, k2, para)
					if candi.Empty() {
						continue
					}

					select {
					case outChan <- candi:
					case <-ctx.Done():
						return ctx.Err()
					}
				}

				return nil
			},
		)
	}

	go func() {
		g.Wait()
		close(outChan)
	}()

	for candi := range outChan {
		roadList = append(roadList, candi)
	}

	if err := g.Wait(); err != nil {
		log.Println(err)
	}

	sort.Slice(roadList, func(i, j int) bool {
		return roadList[i].Records.GetCardinality() > roadList[j].Records.GetCardinality()
	})

	return roadList
}

func findMaxContainment(k1 KmerInfo, kiList KmerInfoList) KmerInfo {
	maxJC := uint64(0)
	var maxK2 KmerInfo
	if k1.RecordSet == nil || k1.RecordSet.IsEmpty() {
		return maxK2
	}

	if len(kiList) == 0 {
		return maxK2
	}

	// fmt.Fprintf(os.Stderr, "begin findMaxContainment, left size: %d\n", len(kiList))

	// kiList sorted already
	for _, k2 := range kiList {
		if k2.RecordSet == nil {
			break
		}

		if k2.RecordSet.IsEmpty() {
			break
		}

		if k2.Empty() {
			break
		}

		if maxJC >= k2.RecordSet.GetCardinality() {
			// no need for following items
			break
		}

		jc := JaccardContainment(k1.RecordSet, k2.RecordSet)
		if jc > maxJC {
			maxJC = jc
			maxK2 = k2
		}

		// fmt.Fprintf(os.Stderr, "maxJC: k1: %d, k2: %d, %.2f\n", k1.RecordSet.Cardinality(), k2.RecordSet.Cardinality(), jc)
	}

	return maxK2
}

type KmerInfo struct {
	K         uint64
	RecordSet *roaring64.Bitmap
}

func (a *KmerInfo) Empty() bool {
	if a == nil || a.RecordSet == nil || a.RecordSet.IsEmpty() {
		return true
	}

	return false
}

type KmerInfoList []KmerInfo

type Para struct {
	Out           string
	Input         string
	Kvalue        int
	KmerEval      int
	MinSize       int
	MaxSize       int
	CPU           int
	MaxRun        int
	MinGC         float64
	MaxGC         float64
	Force         bool
	RecordsSample int // kmer sample from the first n records
	Intersect     bool
	Test          bool
}

func outExists(out string) bool {
	if _, err := os.Stat(out); os.IsNotExist(err) {
		return false
	}

	fh, err := os.Open(out)
	if err != nil {
		log.Println(err)
		return false
	}
	defer fh.Close()

	scanner := bufio.NewScanner(fh)

	count := 0
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if len(line) == 0 || line[0] == '#' {
			continue
		}

		count++
	}

	if err := scanner.Err(); err != nil {
		log.Println(err)
		return false
	}

	return count > 0
}

// work flow
func searchPairedKmers(para Para) {
	startTime := time.Now()

	logFh, err := os.OpenFile(para.Input+".paired_kmers.log", os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	checkErr(err)
	defer logFh.Close()

	log.SetOutput(logFh)

	log.Println("Cmd: ", strings.Join(os.Args, " "))

	log.Printf("records sample: sample count: %d, union: %t", para.RecordsSample, !para.Intersect)

	if len(para.Out) == 0 {
		para.Out = fmt.Sprintf("%s.paired_kmers.bed", para.Input)
	}

	if !para.Force {
		if outExists(para.Out) {
			log.Printf("%s file exists and return", para.Out)
			return
		}
	}

	kiList, allRecords, fastaHash := buildKmerInfo(para.Input, int(para.Kvalue), para.CPU, para)

	// 抽样，否则太多了: 理论基础：保守区域应该是排在前列的；潜在问题，可能不够用。
	if para.KmerEval != -1 {
		if len(kiList) > para.KmerEval {
			kiList = kiList[:para.KmerEval]
		}
	}

	allRecordsSize := allRecords.GetCardinality()

	log.Printf("build kmer info done: %d records, %d kiList", allRecordsSize, len(kiList))
	fmt.Fprintf(os.Stderr, "build kmer info done.\n")

	fo, err := os.Create(para.Out)
	checkErr(err)
	defer fo.Close()

	fmid, err := os.Create(para.Out + ".mid.bed")
	checkErr(err)
	defer fmid.Close()

	fmt.Fprintf(fo, "# %s\n", strings.Join(os.Args, " "))
	fmt.Fprintln(fo, "# inner position")

	roadList := findPairRoads(kiList, fastaHash, para)
	for _, road := range roadList {
		fmt.Fprintf(fo, "%s\t%d\t%d\t%d\t%.2f\t%s\t%s\n", fastaHash[road.Index].ID, road.F3, road.R3, road.Records.GetCardinality(), float64(road.Records.GetCardinality())/float64(allRecordsSize)*100, road.K1, road.K2)

		mid := int((road.F3 + road.R3) / 2)
		fmt.Fprintf(fmid, "%s\t%d\t%d\t%d\t%.2f\n", fastaHash[road.Index].ID, mid, mid+1, road.Records.GetCardinality(), float64(road.Records.GetCardinality())/float64(allRecordsSize)*100)
	}

	log.Printf("Total time used: %s", time.Since(startTime).String())
}

func main() {
	var input *string = flag.StringP("in", "i", "", "[*] input virus/bacterial file in fasta format")
	var out *string = flag.StringP("out", "o", "", "output file")
	var kvalue *int = flag.IntP("kvalue", "k", 17, "k value")
	var recordsSample *int = flag.Int("recordsSample", 10, "kmer sample from the first n records, set -1 for all records")
	var intersect *bool = flag.Bool("intersect", false, "intersect (not union) the kmers of the first n records")
	var kmerEval *int = flag.IntP("kmerEval", "n", 2000, "the first n kmers (sorted by hit records) for next evaluation; set -1 for all kmers")
	var minSize *int = flag.IntP("minSize", "s", 60, "min distance between the two paired-kmers")
	var maxSize *int = flag.IntP("maxSize", "S", 130, "max distance between the two paired-kmers")
	var cpu *int = flag.IntP("cpu", "c", 64, "CPU number to use")
	var maxRun *int = flag.Int("maxRun", 5, "runs are repeated nucleotides (e.g. TAAAAAGC has a 5 bp run of Adenine)")
	var minGC *float64 = flag.Float64P("minGC", "g", 25.0, "max gc allowed for each kmer")
	var maxGC *float64 = flag.Float64P("maxGC", "G", 75.0, "min gc allowed for each kmer")
	var force *bool = flag.BoolP("force", "f", false, "force to override exists output files")
	var test *bool = flag.BoolP("test", "t", false, "test model")
	var help *bool = flag.BoolP("help", "h", false, "print this message")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n\n", os.Args[0])

		fmt.Fprintln(os.Stderr, "./paired-kmer -i sars2.fa -o sars2.paired-kmer.bed -k 17 -s 60 -S 130")

		fmt.Fprintln(os.Stderr)

		flag.PrintDefaults()
	}

	flag.CommandLine.SortFlags = false

	flag.Parse()

	if *help || len(*input) == 0 {
		flag.Usage()
		return
	}

	p := Para{
		Input:         *input,
		Out:           *out,
		Kvalue:        *kvalue,
		RecordsSample: *recordsSample,
		Intersect:     *intersect,
		KmerEval:      *kmerEval,
		MinSize:       *minSize,
		MaxSize:       *maxSize,
		CPU:           *cpu,
		MaxRun:        *maxRun,
		MinGC:         *minGC,
		MaxGC:         *maxGC,
		Force:         *force,
		Test:          *test,
	}

	searchPairedKmers(p)
}
