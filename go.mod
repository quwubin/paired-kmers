module github.com/quwubin/paired-kmers

go 1.18

require (
	github.com/RoaringBitmap/roaring v1.2.1
	github.com/neilotoole/errgroup v0.1.6
	github.com/shenwei356/bio v0.7.1
	github.com/shenwei356/kmers v0.1.0
	github.com/spf13/pflag v1.0.5
)

require (
	github.com/bits-and-blooms/bitset v1.2.0 // indirect
	github.com/davecgh/go-spew v1.1.1 // indirect
	github.com/klauspost/compress v1.15.0 // indirect
	github.com/klauspost/pgzip v1.2.5 // indirect
	github.com/mschoch/smat v0.2.0 // indirect
	github.com/shenwei356/util v0.5.0 // indirect
	github.com/shenwei356/xopen v0.2.1 // indirect
	github.com/ulikunitz/xz v0.5.10 // indirect
	golang.org/x/sync v0.0.0-20210220032951-036812b2e83c // indirect
	gopkg.in/check.v1 v1.0.0-20201130134442-10cb98267c6c // indirect
)

replace github.com/quwubin/bio => ../bio
