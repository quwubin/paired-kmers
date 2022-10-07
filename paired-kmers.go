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
	"runtime"
	"sort"
	"strings"
	"time"

	"github.com/RoaringBitmap/roaring"
	"github.com/RoaringBitmap/roaring/roaring64"
	"github.com/neilotoole/errgroup"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/kmers"
	"github.com/shirou/gopsutil/mem"
	flag "github.com/spf13/pflag"
)

func checkErr(err error) {
	if err != nil {
		log.Fatal(err)
		os.Exit(-1)
	}
}

func exitErr(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, err.Error())
		os.Exit(0)
	}
}

// checkMaxRuns return run number, e.g., return 5 if agcTTTTTagc
func checkMaxRuns(seq []byte) int {
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
func calGC(seq []byte) float64 {
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

		code, err := kmers.Encode(kmer)
		if err != nil {
			log.Printf("encode kmer err: %s", err)
			continue
		}

		if kmersSet.Contains(code) {
			continue
		}

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

		kmersSet.Add(code)
	}

	return kmersSet
}

func JaccardContainment(a, b *roaring.Bitmap) uint64 {
	if a.IsEmpty() || b.IsEmpty() {
		return 0
	}

	return a.AndCardinality(b)
}

type recordInData struct {
	Record *fastx.Record
	Index  uint32
}

type recordOutData struct {
	KmerSet roaring64.Bitmap
	Index   uint32
}

func statsInputDB(db string) (map[uint32]*fastx.Record, uint64, error) {
	fastaHash := make(map[uint32]*fastx.Record)

	reader, err := fastx.NewDefaultReader(db)
	if err != nil {
		return fastaHash, 0, err
	}

	size := uint64(0)

	var record *fastx.Record
	index := uint32(0)
	for {
		record, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Println(err)
			break
		}

		fastaHash[index] = record.Clone()
		size += uint64(len(record.Seq.Seq))
		index++
	}

	aveSize := uint64(float64(size) / float64(index))

	return fastaHash, aveSize, nil
}

// kmer module
func buildKmerInfo(fastaHash map[uint32]*fastx.Record, kvalue int, cpu int, p Para) KmerInfoList {
	g, ctx := errgroup.WithContext(context.Background())
	inChan := make(chan recordInData)

	g.Go(
		func() error {
			defer close(inChan)

			for index, record := range fastaHash {
				if p.Test {
					fmt.Printf("record: %d\n", index)
				}

				select {
				case inChan <- recordInData{Record: record, Index: index}:
				case <-ctx.Done():
					return ctx.Err()
				}
			}

			return nil
		},
	)

	outChan := make(chan recordOutData)

	for i := 0; i < cpu; i++ {
		g.Go(
			func() error {
				for t := range inChan {
					seq := bytes.ToUpper(t.Record.Seq.Seq)
					t.Record.Seq.Seq = seq

					kmerSet := buildKmer(seq, kvalue, p)

					if p.Test {
						fmt.Printf("record: %d, kmer: %d\n", t.Index, kmerSet.GetCardinality())
					}

					select {
					case outChan <- recordOutData{KmerSet: *kmerSet, Index: t.Index}:
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

	kmerInfoMap := make(map[uint64]*roaring.Bitmap)

	processed := 0
	filteredKmerSet := roaring64.New()
	for o := range outChan {
		ks := o.KmerSet.Iterator()
		for ks.HasNext() {
			code := ks.Next()
			if !filteredKmerSet.IsEmpty() {
				if !filteredKmerSet.Contains(code) {
					continue
				}
			}

			_, ok := kmerInfoMap[code]
			if !ok {
				kmerInfoMap[code] = roaring.New()
			}
			kmerInfoMap[code].Add(o.Index)
		}

		processed++

		if filteredKmerSet.IsEmpty() {
			v, err := mem.VirtualMemory()
			exitErr(err)

			// a temporary break
			if processed >= p.RecordCutoff {
				if p.ForceCutoff || v.UsedPercent >= p.Memory {
					log.Printf("processed %d records, use these records's first %d kmers (sorted by hit records) for following analysis", processed, p.KmerCutoff)
					kiList := sortKmerInfoList(kmerInfoMap)
					if len(kiList) > p.KmerCutoff {
						for i := 0; i < p.KmerCutoff; i++ {
							filteredKmerSet.Add(kiList[i].K)
						}
					}

					for code := range kmerInfoMap {
						if !filteredKmerSet.Contains(code) {
							delete(kmerInfoMap, code)
						}
					}
				}
			}
		}

		if p.Test {
			fmt.Printf("processed: %d, o.KmerSet: %d, kmerInfoMap size: %d\n", processed, o.KmerSet.GetCardinality(), len(kmerInfoMap))
		}
	}

	kiList := make(KmerInfoList, len(kmerInfoMap))

	if err := g.Wait(); err != nil {
		log.Println(err)
		return kiList
	}

	kiList = sortKmerInfoList(kmerInfoMap)

	for k := range kmerInfoMap {
		delete(kmerInfoMap, k)
	}

	return kiList
}

func sortKmerInfoList(kmerInfoMap map[uint64]*roaring.Bitmap) KmerInfoList {
	kiList := make(KmerInfoList, len(kmerInfoMap))

	i := 0
	for k, v := range kmerInfoMap {
		kiList[i] = KmerInfo{K: k, RecordSet: v}
		i++
	}

	sort.Slice(kiList, func(i, j int) bool {
		return kiList[i].RecordSet.GetCardinality() > kiList[j].RecordSet.GetCardinality()
	})

	return kiList
}

func fincCandiRoads(fastaHash map[uint32]*fastx.Record, k1 KmerInfo, k2 KmerInfo, para Para, offset int) PairedKmer {
	k1seq := string(kmers.Decode(k1.K, para.Kvalue))
	k2seq := string(kmers.Decode(k2.K, para.Kvalue))
	re1 := regexp.MustCompile(k1seq)
	re2 := regexp.MustCompile(k2seq)
	var candi PairedKmer

	if k1.IsEmpty() || k2.IsEmpty() {
		return candi
	}

	k1.RecordSet.And(k2.RecordSet)

	sharedRecords := roaring.New()
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
					candi.F5 += offset
					candi.F3 += offset
					candi.R5 += offset
					candi.R3 += offset
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
	Index   uint32 // record index
	K1      string
	K2      string
	Dist    int
	F5      int
	F3      int
	R5      int
	R3      int
	Amp5    string
	Amp3    string
	Records *roaring.Bitmap
}

func (a *PairedKmer) Empty() bool {
	return a.K1 == ""
}

type PairedKmers []PairedKmer

type kiData struct {
	K1 KmerInfo
	KI KmerInfoList
}

func findPairRoads(kiList KmerInfoList, fastaHash map[uint32]*fastx.Record, para Para, offset int) PairedKmers {
	g, ctx := errgroup.WithContext(context.Background())
	inChan := make(chan kiData)
	var roadList PairedKmers

	g.Go(
		func() error {
			defer close(inChan)

			for i := 0; i < len(kiList)-1; i++ {
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
					if k2.IsEmpty() {
						continue
					}

					candi := fincCandiRoads(fastaHash, t.K1.Clone(), k2, para, offset)
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
	if k1.IsEmpty() {
		return maxK2
	}

	if len(kiList) == 0 {
		return maxK2
	}

	// fmt.Fprintf(os.Stderr, "begin findMaxContainment, left size: %d\n", len(kiList))

	// kiList sorted already
	for _, k2 := range kiList {
		if k2.IsEmpty() {
			break
		}

		if maxJC >= k2.RecordSet.GetCardinality() {
			// no need for following items
			return maxK2
		}

		jc := JaccardContainment(k1.RecordSet, k2.RecordSet)
		if jc > maxJC {
			maxJC = jc
			maxK2 = k2
		}

		// fmt.Fprintf(os.Stderr, "maxJC: k1: %d, k2: %d, %.2f\n", k1.RecordSet.Cardinality(), k2.RecordSet.Cardinality(), jc)
	}

	return maxK2.Clone()
}

type KmerInfo struct {
	K         uint64
	RecordSet *roaring.Bitmap
}

func (a *KmerInfo) Clone() KmerInfo {
	var b KmerInfo
	b.K = a.K
	b.RecordSet = a.RecordSet.Clone()

	return b
}

func (a *KmerInfo) IsEmpty() bool {
	if a.K == 0 || a.RecordSet.IsEmpty() {
		return true
	}

	return false
}

type KmerInfoList []KmerInfo

type Para struct {
	Out          string // output file name
	OutputCutoff int    // output the maximum number of paired-kmers
	Input        string // input virus/bacterial file in fasta format
	Kvalue       int    // kmer size, default 17, usually, smaller value for virus and big value for bacterial
	RecordCutoff int
	KmerCutoff   int
	ColumnSize   int
	ForceCutoff  bool
	Memory       float64
	MinSize      int
	MaxSize      int
	CPU          int
	MaxRun       int
	MinGC        float64
	MaxGC        float64
	Force        bool
	Test         bool
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

	if len(para.Out) == 0 {
		para.Out = fmt.Sprintf("%s.paired_kmers.bed", para.Input)
	}

	if !para.Force {
		if outExists(para.Out) {
			log.Printf("%s file exists and return", para.Out)
			return
		}
	}

	fastaHash, aveSize, err := statsInputDB(para.Input)
	exitErr(err)

	columnCount := aveSize/uint64(para.ColumnSize) + 1

	totalRecordsCount := len(fastaHash)

	var roadList PairedKmers
	start := 0
	stop := para.ColumnSize
	for i := 1; i <= int(columnCount); i++ {
		columnFastaHash := make(map[uint32]*fastx.Record, len(fastaHash))

		for index, record := range fastaHash {
			record = record.Clone()

			recordSize := len(record.Seq.Seq)
			if recordSize <= start {
				continue
			}

			recordStop := stop

			if recordSize <= recordStop {
				recordStop = recordSize
			}

			record.Seq.Seq = record.Seq.Seq[start:recordStop]

			columnFastaHash[index] = record
		}

		kiList := buildKmerInfo(columnFastaHash, int(para.Kvalue), para.CPU, para)

		log.Printf("column: %d/%d, range: %d-%d, build kmer info done: %d records, %d kiList", i, columnCount, start, stop, totalRecordsCount, len(kiList))

		subRoadList := findPairRoads(kiList, columnFastaHash, para, start)

		log.Printf("column: %d/%d, range: %d-%d, found %d/%d sub-road list", i, columnCount, start, stop, len(subRoadList), len(roadList))

		roadList = append(roadList, subRoadList...)

		sort.Slice(roadList, func(i, j int) bool {
			return roadList[i].Records.GetCardinality() > roadList[j].Records.GetCardinality()
		})

		if len(roadList) > para.OutputCutoff {
			roadList = roadList[:para.OutputCutoff]
		}

		start = stop
		stop += para.ColumnSize
	}

	sort.Slice(roadList, func(i, j int) bool {
		return roadList[i].Records.GetCardinality() > roadList[j].Records.GetCardinality()
	})

	printOut(roadList, para, fastaHash)

	log.Printf("Total time used: %s", time.Since(startTime).String())
}

func printOut(roadList PairedKmers, p Para, fastaHash map[uint32]*fastx.Record) {
	fo, err := os.Create(p.Out)
	checkErr(err)
	defer fo.Close()

	fmid, err := os.Create(p.Out + ".mid.bed")
	checkErr(err)
	defer fmid.Close()

	fmt.Fprintf(fo, "# %s\n", strings.Join(os.Args, " "))
	fmt.Fprintln(fo, "# inner position")

	totalRecordsCount := len(fastaHash)

	for _, road := range roadList {
		fmt.Fprintf(fo, "%s\t%d\t%d\t%d\t%.2f\t%s\t%s\n", fastaHash[road.Index].ID, road.F3, road.R3, road.Records.GetCardinality(), float64(road.Records.GetCardinality())/float64(totalRecordsCount)*100, road.K1, road.K2)

		mid := int((road.F3 + road.R3) / 2)
		fmt.Fprintf(fmid, "%s\t%d\t%d\t%d\t%.2f\n", fastaHash[road.Index].ID, mid, mid+1, road.Records.GetCardinality(), float64(road.Records.GetCardinality())/float64(totalRecordsCount)*100)
	}

}

func main() {
	var input *string = flag.StringP("in", "i", "", "[*] input virus/bacterial file in fasta format")
	var out *string = flag.StringP("out", "o", "", "output file")
	var outputCutoff *int = flag.Int("outputCutoff", 1000, "output the maximum number of paired-kmers")
	var kvalue *int = flag.IntP("kvalue", "k", 17, "k value")
	var recordCutoff *int = flag.Int("recordCutoff", 100, "the minimum records to be used to sort kmers, if memory is below of the threshold, this limit won't work. Add -F to force use this cutoff.")
	var kmerCutoff *int = flag.IntP("kmerCutoff", "n", 10000, "the first n kmers (sorted by hit records) for next evaluation, if memory is below of the threshold, this limit won't work. Add -F to force use this cutoff.")
	var columnSize *int = flag.Int("columnSize", 100000, "max column size for large size of genome")
	var forceCutoff *bool = flag.BoolP("forceCutoff", "F", false, "force to use --kmerCutoff for kmer/record cutoff even the memory usage is below the threshold")
	var memory *float64 = flag.Float64P("memory", "m", 70, "memory percent to use, e.g., 70 for 70%")
	var minSize *int = flag.IntP("minSize", "s", 60, "min distance between the two paired-kmers")
	var maxSize *int = flag.IntP("maxSize", "S", 130, "max distance between the two paired-kmers")
	var cpu *int = flag.IntP("cpu", "c", 64, "CPU number to use")
	var maxRun *int = flag.Int("maxRun", 5, "runs are repeated nucleotides (e.g. TAAAAAGC has a 5 bp run of Adenine)")
	var minGC *float64 = flag.Float64P("minGC", "g", 25.0, "max gc allowed for each kmer")
	var maxGC *float64 = flag.Float64P("maxGC", "G", 75.0, "min gc allowed for each kmer")
	var force *bool = flag.BoolP("force", "f", false, "force to override exists output files")
	var test *bool = flag.BoolP("test", "t", false, "test model for debugging")
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
		Input:        *input,
		Out:          *out,
		OutputCutoff: *outputCutoff,
		Kvalue:       *kvalue,
		RecordCutoff: *recordCutoff,
		KmerCutoff:   *kmerCutoff,
		ColumnSize:   *columnSize,
		ForceCutoff:  *forceCutoff,
		Memory:       *memory,
		MinSize:      *minSize,
		MaxSize:      *maxSize,
		CPU:          *cpu,
		MaxRun:       *maxRun,
		MinGC:        *minGC,
		MaxGC:        *maxGC,
		Force:        *force,
		Test:         *test,
	}

	if p.ForceCutoff {
		if p.KmerCutoff < 0 {
			fmt.Fprintln(os.Stderr, "-F requires a positive kmer cutoff for -n parameters")
			return
		}
	}

	runtime.GOMAXPROCS(p.CPU)

	searchPairedKmers(p)
}
