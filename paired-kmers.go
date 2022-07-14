package main

/*
Copyright © 2019 Wubin Qu <quwubin@gmail.com>
*/

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"sort"
	"strings"
	"time"

	mapset "github.com/deckarep/golang-set/v2"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	flag "github.com/spf13/pflag"
)

func checkErr(err error) {
	if err != nil {
		log.Fatal(err)
		os.Exit(-1)
	}
}

// kmer module
func BuildKmerSet(seq string, k int) mapset.Set[string] {
	kmers := mapset.NewSet[string]()
	nKmers := len(seq) - k + 1

	if nKmers <= 0 {
		return kmers
	}

	for i := 0; i < nKmers; i++ {
		kmer := seq[i : i+k]
		kmers.Add(kmer)
	}

	return kmers
}

func JaccardContainment(a, b mapset.Set[string]) float64 {
	return float64(a.Intersect(b).Cardinality()) / float64(a.Cardinality())
}

func JaccardContainment2(a, b mapset.Set[string]) int {
	return a.Intersect(b).Cardinality()
}

// kmer module

func buildKmerInfo(db string, kvalue int) (KmerInfoList, mapset.Set[string], map[string]*fastx.Record) {
	kmerInfo := make(map[string]mapset.Set[string])
	allRecords := mapset.NewSet[string]()
	fastaHash := make(map[string]*fastx.Record)

	reader, err := fastx.NewDefaultReader(db)
	checkErr(err)
	var record *fastx.Record
	for {
		record, err = reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkErr(err)
			break
		}

		allRecords.Add(string(record.ID))

		seq := bytes.ToUpper(record.Seq.Seq)

		record.Seq.Seq = seq

		fastaHash[string(record.ID)] = record

		kmers := BuildKmerSet(string(seq), kvalue)

		for k := range kmers.Iterator().C {
			_, ok := kmerInfo[k]
			if !ok {
				kmerInfo[k] = mapset.NewSet[string]()
			}

			kmerInfo[k].Add(string(record.ID))
		}
	}

	kiList := make(KmerInfoList, len(kmerInfo))
	i := 0
	for k, v := range kmerInfo {
		kiList[i] = KmerInfo{K: k, RecordSet: v}
		i++
	}

	sort.Slice(kiList, func(i, j int) bool { return kiList[i].RecordSet.Cardinality() > kiList[j].RecordSet.Cardinality() })

	return kiList, allRecords, fastaHash
}

func fincCandiRoads(fastaHash map[string]*fastx.Record, k1 KmerInfo, k2 KmerInfo, para Para) PairedKmer {
	re1 := regexp.MustCompile(k1.K)
	re2 := regexp.MustCompile(k2.K)
	var candi PairedKmer
	sharedRecords := mapset.NewSet[string]()
	for recordName := range k1.RecordSet.Iterator().C {
		if !k2.RecordSet.Contains(recordName) {
			continue
		}

		recordSeq := string(fastaHash[recordName].Seq.Seq)
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
					k2seq, err := seq.NewSeq(seq.DNA, []byte(k2.K))
					if err != nil {
						log.Println(err)
						continue
					}

					candi = PairedKmer{
						Name: recordName,
						K1:   k1.K,
						K2:   string(k2seq.RevCom().Seq),
						F5:   f[0],     // 包含
						F3:   f[1] - 1, // 包含
						R5:   r[1],     // 不包含
						R3:   r[0],     // 包含
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
			sharedRecords.Add(recordName)
		}
	}

	candi.Records = sharedRecords

	return candi
}

type PairedKmer struct {
	Name    string
	K1      string
	K2      string
	Dist    int
	F5      int
	F3      int
	R5      int
	R3      int
	Amp5    string
	Amp3    string
	Records mapset.Set[string]
}

func (a *PairedKmer) Empty() bool {
	if a.Name == "" {
		return true
	}

	return false
}

type PairedKmers []PairedKmer

func findPairRoads(kiList KmerInfoList, fastaHash map[string]*fastx.Record, para Para) PairedKmers {
	var roadList PairedKmers
	for i := 0; i < len(kiList); i++ {
		if i == len(kiList)-1 {
			break
		}

		k1 := kiList[i]

		k2 := findMaxContainment(k1, kiList, i+1)
		if k2.Empty() {
			continue
		}

		candi := fincCandiRoads(fastaHash, k1, k2, para)
		if candi.Empty() {
			continue
		}

		roadList = append(roadList, candi)
	}

	sort.Slice(roadList, func(i, j int) bool { return roadList[i].Records.Cardinality() > roadList[j].Records.Cardinality() })

	return roadList
}

func findMaxContainment(k1 KmerInfo, kiList KmerInfoList, j int) KmerInfo {
	maxJC := 0
	var maxK2 KmerInfo
	if k1.RecordSet.Cardinality() == 0 {
		return maxK2
	}

	if len(kiList[j:]) == 0 {
		return maxK2
	}

	// fmt.Fprintf(os.Stderr, "begin findMaxContainment, left size: %d\n", len(kiList))

	// kiList sorted already
	for _, k2 := range kiList[j:] {
		if k2.RecordSet.Cardinality() == 0 {
			break
		}

		if maxJC >= k2.RecordSet.Cardinality() {
			// no need for following items
			break
		}

		jc := JaccardContainment2(k1.RecordSet, k2.RecordSet)
		if jc > maxJC {
			maxJC = jc
			maxK2 = k2
		}

		// fmt.Fprintf(os.Stderr, "maxJC: k1: %d, k2: %d, %.2f\n", k1.RecordSet.Cardinality(), k2.RecordSet.Cardinality(), jc)
	}

	return maxK2
}

type KmerInfo struct {
	K         string
	RecordSet mapset.Set[string]
}

func (a *KmerInfo) Empty() bool {
	if a == nil || a.K == "" {
		return true
	}

	return false
}

type KmerInfoList []KmerInfo

func searchPairedKmers(para Para) {
	logFh, err := os.OpenFile(para.Input+".paired_kmers.log", os.O_RDWR|os.O_CREATE|os.O_APPEND, 0666)
	checkErr(err)
	defer logFh.Close()

	log.SetOutput(logFh)

	log.Println("create log file done")

	log.Println("Cmd: ", strings.Join(os.Args, " "))

	kiList, allRecords, fastaHash := buildKmerInfo(para.Input, int(para.Kvalue))

	allRecordsSize := allRecords.Cardinality()

	fmt.Fprintf(os.Stderr, "build kmer info done.\n")

	if len(para.Out) == 0 {
		para.Out = fmt.Sprintf("%s.%d.%d.%d.paired_kmers.bed", para.Input, para.Kvalue, para.MinSize, para.MaxSize)
	}

	fo, err := os.Create(para.Out)
	checkErr(err)
	defer fo.Close()

	startTime := time.Now()

	roadList := findPairRoads(kiList, fastaHash, para)
	for _, road := range roadList {
		fmt.Fprintf(fo, "%s\t%d\t%d\t%d\t%.2f\n", road.Name, road.F3, road.R3, road.Records.Cardinality(), float64(road.Records.Cardinality())/float64(allRecordsSize)*100)
	}

	log.Printf("Total time used: %s", time.Since(startTime).String())
}

type Para struct {
	Out     string
	Input   string
	Kvalue  int
	Number  int
	MinSize int
	MaxSize int
}

func main() {
	var input *string = flag.StringP("in", "i", "", "[*] input virus/bacterial file in fasta format")
	var out *string = flag.StringP("out", "o", "", "output file")
	var kvalue *int = flag.IntP("kvalue", "k", 21, "k value")
	var number *int = flag.IntP("number", "n", 100, "total n paired-kmers will be print out")
	var minSize *int = flag.IntP("minSize", "s", 60, "min size of two paired-kmers")
	var maxSize *int = flag.IntP("maxSize", "S", 130, "max size of two paired-kmers")
	var help *bool = flag.BoolP("help", "h", false, "help")

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n\n", os.Args[0])

		fmt.Fprintln(os.Stderr, "./paired-kmer -i sars2.fa -o sars2.paired-kmer.txt -n 20 -s 90 -S 200")

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
		Input:   *input,
		Out:     *out,
		Kvalue:  *kvalue,
		Number:  *number,
		MinSize: *minSize,
		MaxSize: *maxSize,
	}

	searchPairedKmers(p)
}
