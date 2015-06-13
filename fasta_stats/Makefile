# Makefile for chrom-gc

CXX  = llvm-g++ # g++
CC   = llvm-g++ # g++

CXXINCLUDEDIR = 
CXXFLAGS = $(CXXINCLUDEDIR) -D_FILE_OFFSET_BITS=64 -Wall -ggdb -g3 -O0 -fno-inline -fno-eliminate-unused-debug-types
LIBS = -lz
RM = rm -f

BIN  = fasta_stats

OBJ  = fasta_stats.o

HEADER = FastaStats.h

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

$(HEADER): kseq.h

clean: 
	$(RM) $(OBJ) $(BIN)

.PHONY: all all-before all-after clean clean-custom
