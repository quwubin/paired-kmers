# Paired-Kmers: finding conserved regions in highly diverse genomes

## Introduction

It is a challenge to identify long enough conserved regions for PCR primer design in highly diverse genomes, such as virus and bacteria. Both multiple sequence alignment and k-mer based alignment free methods are utilized for the task. However, with more and more genome sequences are available, the computational and computer memory requirements are stll the bottleneck.

Here, we propose a fast and memory-efficient k-mer based method, which we named Paired-Kmers, to identify two regions that are close to each other but not necessarily contiguous for picking primers. Importantly, we implemented column-based and row-based reduction methods to overcome memory bottleneck when analyzing species with large genome sizes and species with millions of different strains. 

Paired-Kmers is a general tool to identify conserved regions from microbial genomes.

## Installation

Download the latest pre-compiled binaries from: https://github.com/quwubin/paired-kmers/releases 

```bash
tar zxvf paired-kmer-xxx-amd64.tar.bz

./paired-kmers -h
```


## Usage

```bash
./paired-kmers -k 17 -c 16 -s 90 -S 300 -i test/hpv16.fa -o test/hpv16.paired_kmers.bed
```

## Output

There are two output files, one is "hpv16.paired_kmer.bed", and the other is "hpv16.paired_kmer.bed.mid.bed".

For each file, each row is a conserved region represented by one chrosome location. There are 8 columns in file "hpv16.paired_kmer.bed".
1. The 1st column is the represented chromosome;
2. The 2nd column is the inner position (3' end) of left kmer;
3. The 3rd column is the inner position (5' end) of right kmer;
4. The 4th column is the number of covered records of this conserved paired-kmer;
5. The 5th column is the records covered percent of this conserved paired-kmer;
6. The 6th is the left kmer sequence;
7. The 7th is the right kmer sequence;
8. The 8th is the records name list seperated by space.

Instead of the inner location of paired-kmers, file "hpv16.paired_kmer.bed.mid.bed" stores the middle position of the paired-kmers. This file is used for program (or scripts) reading and further processing.

```bash
# paired-kmers -i hpv16.fa -o hpv16.paired_kmer.bed -k 17 -s 90 -S 300 -f
# inner position
# chr   innerStart      innerEnd        hitNumber       hitPercent(%)   leftKmer        rightKmer       hitRecords
JN565303.1      2251    2477    506     100.00  GGTGATTGGAAGCAAAT       CATCTAACATACCTATT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2759    2849    506     100.00  GTCAGGACAAAATACTA       GGATTTCCAGTTCTTAT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2032    2130    506     100.00  GTAAAGGATTGTGCAAC       ATGTCATTATCGTAGGC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      708     846     506     100.00  AGCAGAACCGGACAGAG       GTAGATTATGGTTTCTG       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      427     491     506     100.00  GTCGGTGGACCGGTCGA       GCTTTTGACAGTTAATA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      1862    2063    506     100.00  TGCACAATTGGCAGACA       CAATTTTGGAGGCTCTA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2487    2737    506     100.00  TGCCAAAATAGGTATGT       CCAGTTCTTATCATTAA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2487    2748    506     100.00  AACTGGAAATCCTTTTT       ACATACCTATTTTGGCA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      421     687     506     100.00  GGACAAGCAGAACCGGA       GACAGTTAATACACCTA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      398     552     506     100.00  GCTGTAATCATGCATGG       AAATCACACAACGGTTT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2691    2923    506     100.00  GGTGGTGTTTACATTTC       GCACATTCTAGGCGCAT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      209     410     506     100.00  GTATTAACTGTCAAAAG       TGTTGCTTGCAGTACAC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2620    2756    506     100.00  ATCCTTTTTCTCAAGGA       GTAATTAATAATGGAGG       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2133    2243    506     100.00  GAAGCAAATTGTTATGT       TTACAATTTTTGCCTGT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      1576    1734    506     100.00  TCATGGGGAATGGTTGT       ACCCCGTATAACTCTTT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2689    2831    506     100.00  GCCAACGTTTAAATGTG       AATGTAAACACCACCAA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2050    2232    506     100.00  GGAGGTGATTGGAAGCA       ATTTCACTATCGTCTAC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      616     688     506     100.00  GACAAGCAGAACCGGAC       GTCTCTGGTTGCAAATC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      422     668     506     100.00  TGAAATAGATGGTCCAG       TGACAGTTAATACACCT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2494    2628    506     100.00  ATAGGTATGTTAGATGA       TCTGTACCAGCATTAAT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      1625    1844    506     100.00  GATAGAGCCTCCAAAAT       ACACGTTGATTTATTAC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      1392    1613    506     100.00  TAAATCAACGTGTTGCG       GTGTTTCAGTCTCATGG       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      510     734     506     100.00  TTGCAAGTGTGACTCTA       ACATCGACCGGTCCACC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2490    2746    506     100.00  CAAAATAGGTATGTTAG       AAAGGATTTCCAGTTCT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      392     577     506     100.00  TACAACAAACCGTTGTG       ATATTCATGCAATGTAG       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      1859    2039    506     100.00  CGATAGTGAAATTGCAT       TTTTGGAGGCTCTATCA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      396     596     506     100.00  ACAAACCGTTGTGTGAT       CTGGTTGCAAATCTAAC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2763    2839    506     100.00  TTAAATGTGTGTCAGGA       AAAAGGATTTCCAGTTC       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      3203    3375    506     100.00  GACATATGCAATACAAT       CCTGACCACCCGCATGA       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      2484    2577    506     100.00  CATAGACCATTGGTACA       TACCTATTTTGGCATCT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
JN565303.1      512     767     506     100.00  TGGACCGGTCGATGTAT       CTACGTGTGTGCTTTGT       JN565303.1 KY549224.1 KY549225.1 JQ004093.1 JN565302.1 JQ067943.1 JQ067944.1 KY549166.1 KY549200.1 KY549211.
```

### Author notes about the output

1. Hit number and hit percent means that from the region we can find primers who will amplify the corresponding hit targets. 
2. Primers should selected around the region. If we have primer design workflow, I suggest that pick and check primers from the middle point of the region (also provided by file file suffix ".mid.bed").
3. Paired-kmers has options for limiting kmers with proper GC content and skipping repeated nucleotides. However, these options are used for reducting kmers to control memory when finding conserved regions in some very large genome database. So, the output should not used directly as primers, even they are paired. Good primers should also be checked for dimers, non-specificity etc. My previous work MFEprimer (http://academic.oup.com/nar/article/47/W1/W610/5486745) is designed for primers quality check if you need.


### Note about genome sequences

We compare sequences to find the conserved regions, so we recommend:

1. use complete genome sequences in priority;
2. collect as many as sequences as possible;
3. scaffold or partial sequences are valuable when few of comple genomes are available.

## Issues and Bugs

https://github.com/quwubin/paired-kmers/issues

## Cite us

Manuscript is submitted.

## Contact

Wubin Qu (<quwubin@gmail.com>)

Haoyang Cai (<haoyang.cai@scu.edu.cn>)