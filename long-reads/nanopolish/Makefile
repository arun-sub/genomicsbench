# dir
DIR = ../../tools/nanopolish
SUBDIRS := src src/hmm src/thirdparty src/thirdparty/scrappie src/common src/alignment src/pore_model src/io
SUBDIRS_INCLUDE = $(foreach dir, $(SUBDIRS), $(addprefix -I${DIR}/, $(SUBDIRS)))
#
# Set libraries, paths, flags and options
#

#Basic flags every build needs
LIBS = -lz      # libz is Zlib, a compression library which can do gzip, deflate, and a few other formats
CXXFLAGS ?= -g -O0      # ?= means only if it's not set/
CXXFLAGS += -std=c++11 -fopenmp -fsigned-char -D_FILE_OFFSET_BITS=64 #D_FILE_OFFSET_BITS=64 makes nanopolish work in 32 bit systems
LDFLAGS ?=
CXX ?= g++

# Check operating system, OSX doesn't have -lrt
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LIBS += -lrt
endif


H5_LIB = ${DIR}/lib/libhdf5.a
H5_INCLUDE = -I${DIR}/include    # adds a directory to the list of places searched by the compiler for a file named on a #include line
LIBS += -ldl

EIGEN_CHECK = ${DIR}/eigen/INSTALL
EIGEN_INCLUDE = -I${DIR}/eigen/

HTS_LIB = ${DIR}/htslib/libhts.a
HTS_INCLUDE = -I${DIR}/htslib

MINIMAP2_LIB = ${DIR}/minimap2/libminimap2.a
MINIMAP2_INCLUDE = -I${DIR}/minimap2

ifeq ($(ARM), 1)
    MINIMAP2_OPT=arm_neon=1
else
    MINIMAP2_OPT=
endif

# Include the header-only fast5 library
FAST5_INCLUDE = -I${DIR}/fast5/include

# Include the vtune lib
VTUNE_HOME = /opt/intel/vtune_profiler_2020
VTUNE_INCLUDE = -I$(VTUNE_HOME)/include
VTUNE_LIB = -littnotify -L$(VTUNE_HOME)/lib64
CPPFLAGS += $(VTUNE_INCLUDE)
LIBS += $(VTUNE_LIB)

# Add include flags
CPPFLAGS += $(H5_INCLUDE) $(HTS_INCLUDE) $(MINIMAP2_INCLUDE) $(FAST5_INCLUDE) $(EIGEN_INCLUDE) ${SUBDIRS_INCLUDE}

# Main programs to build
RUN_ALIGN = run_align

.PHONY: all
all: depend $(RUN_ALIGN)

#
# Source files
#

# Find the source files by searching subdirectories
CPP_SRC := $(foreach dir, $(SUBDIRS), $(wildcard ${DIR}/$(dir)/*.cpp))
C_SRC := $(foreach dir, $(SUBDIRS), $(wildcard ${DIR}$(dir)/*.c))

# Automatically generated object names
CPP_OBJ = $(CPP_SRC:.cpp=.o)
C_OBJ = $(C_SRC:.c=.o)

# Generate dependencies
.PHONY: depend
depend: .depend

.depend: $(CPP_SRC)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM run_align.cpp > ./.depend;

# Compile objects
.cpp.o:
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) -fPIC $<

$(RUN_ALIGN): run_align.o $(HTS_LIB) $(MINIMAP2_LIB) $(H5_LIB) $(EIGEN_CHECK)
	echo $(CPP_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) -fPIC $< $(DIR)/libnanopolish.a $(HTS_LIB) $(MINIMAP2_LIB) $(H5_LIB) $(LIBS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(RUN_ALIGN)
