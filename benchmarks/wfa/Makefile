###############################################################################
# Flags & Folders
###############################################################################
FOLDER_BIN:=bin
FOLDER_BUILD:=build

# UNAME=$(shell uname)

CC=gcc
CXX=g++

ifeq ($(arch),sse41)
	ARCH_FLAGS=-msse4.1
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else	
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else	
		ARCH_FLAGS=-mavx512bw
	endif
else ifeq ($(arch),native)
		ARCH_FLAGS=-march=native
else ifneq ($(arch),)
	## To provide a different architecture flag like -march=core-avx2.
		ARCH_FLAGS=$(arch)
endif

LD_FLAGS=-lm
CC_FLAGS=-Wall -fopenmp $(ARCH_FLAGS) # -g
# ifeq ($(UNAME), Linux)
#   LD_FLAGS+=-lrt 
# endif

AR=ar
AR_FLAGS=-rsc

###############################################################################
# Compile rules
###############################################################################
SUBDIRS=gap_affine \
        utils
       
LIB_WFA=$(FOLDER_BUILD)/libwfa.a

all: CC_FLAGS+=-O3
all: MODE=all
all: setup
all: $(SUBDIRS) tools $(LIB_WFA)

debug: setup
debug: MODE=all
debug: $(SUBDIRS) tools $(LIB_WFA)

$(LIB_WFA): .FORCE
	$(AR) $(AR_FLAGS) $(LIB_WFA) $(FOLDER_BUILD)/*.o 2> /dev/null

setup:
	@mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD)

clean:
	rm -rf $(FOLDER_BIN) $(FOLDER_BUILD)
	
###############################################################################
# Subdir rule
###############################################################################
export
$(SUBDIRS):
	$(MAKE) --directory=$@ all
	
tools:
	$(MAKE) --directory=$@ $(MODE)

.PHONY: $(SUBDIRS) tools
.FORCE:

