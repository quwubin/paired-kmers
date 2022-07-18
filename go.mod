module github.com/quwubin/paired-kmers

go 1.18

require (
	github.com/deckarep/golang-set/v2 v2.1.0
	github.com/neilotoole/errgroup v0.1.6
	github.com/quwubin/bio v0.0.0-00010101000000-000000000000
	github.com/shenwei356/bio v0.7.1
	github.com/spf13/pflag v1.0.5
)

require (
	github.com/edsrzf/mmap-go v1.1.0 // indirect
	github.com/ka-weihe/fast-levenshtein v0.0.0-20201227151214-4c99ee36a1ba // indirect
	github.com/klauspost/compress v1.15.0 // indirect
	github.com/klauspost/pgzip v1.2.5 // indirect
	github.com/montanaflynn/stats v0.6.6 // indirect
	github.com/shenwei356/util v0.5.0 // indirect
	github.com/shenwei356/xopen v0.2.1 // indirect
	github.com/ulikunitz/xz v0.5.10 // indirect
	golang.org/x/sync v0.0.0-20210220032951-036812b2e83c // indirect
	golang.org/x/sys v0.0.0-20220503163025-988cb79eb6c6 // indirect
	gopkg.in/check.v1 v1.0.0-20201130134442-10cb98267c6c // indirect
)

replace github.com/quwubin/bio => ../bio
