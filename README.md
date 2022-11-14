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
paired-kmers -k 17 -c 16 -s 90 -S 300 -i hpv16.fa -o hpv16.paired-kmers.bed
```

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