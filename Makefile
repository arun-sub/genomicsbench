CXX=g++
CC=gcc
ARCH=avx2
VTUNE_HOME=/opt/intel/oneapi/vtune/2021.1.1
MKLROOT=/opt/intel/oneapi/mkl/2021.1.1
MKL_IOMP5_DIR=/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin

.PHONY: clean

all:
	$(info Starting build..this may take a while..)
	cd tools/htslib && autoreconf -i && ./configure && $(MAKE)
	cd tools/bwa-mem2; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/fmi; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/bsw; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/dbg; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd tools/GKL; ./gradlew test
	cd benchmarks/phmm; $(MAKE) CC=$(CC) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd tools/minimap2; $(MAKE)
	cd benchmarks/chain; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd tools/spoa; mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Release ..; $(MAKE)
	cd benchmarks/poa; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/pileup; $(MAKE) CC=$(CC) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/kmer-cnt; $(MAKE) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME)
	cd benchmarks/grm/2.0/build_dynamic; $(MAKE) CC=$(CC) CXX=$(CXX) arch=$(ARCH) VTUNE_HOME=$(VTUNE_HOME) MKLROOT=$(MKLROOT) MKL_IOMP5_DIR=$(MKL_IOMP5_DIR) #needs MKL
	$(info Successfully built GenomicsBench!)

clean:
	cd tools/bwa-mem2; $(MAKE) clean
	cd benchmarks/fmi; $(MAKE) clean
	cd benchmarks/bsw; $(MAKE) clean
	cd benchmarks/dbg; $(MAKE) clean
	cd tools/GKL; ./gradlew clean
	cd benchmarks/phmm; $(MAKE) clean
	cd tools/minimap2; $(MAKE)
	cd benchmarks/chain; $(MAKE) clean
	cd benchmarks/poa; $(MAKE) clean
	cd benchmarks/pileup; $(MAKE) clean
	cd benchmarks/kmer-cnt; $(MAKE) clean
	cd benchmarks/grm/2.0/build_dynamic; $(MAKE) clean
