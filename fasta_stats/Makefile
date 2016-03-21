# Makefile for chrom-gc

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')

CXX  = g++
CC   = g++

ifeq ($(uname_S),Darwin)
	CXX  = llvm-g++
	CC   = llvm-g++
endif

CXXINCLUDEDIR =
CXXFLAGS =  $(CXXINCLUDEDIR) -D_FILE_OFFSET_BITS=64 -Wall -ggdb -g3 -O0 -fno-inline -fno-eliminate-unused-debug-types
LIBS = -lz
RM = rm -f

BIN  = fasta_stats

OBJ  = fasta_stats.o

HEADER = FastaStats.h

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

# $(HEADER): kseq.h  # no longer depends on kseq.h, it is included complete

clean: 
	$(RM) $(OBJ) $(BIN)

.PHONY: all all-before all-after clean clean-custom
