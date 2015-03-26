# Makefile for chrom-gc

CPP  = g++
CC   = g++

CXXINCLUDEDIR = 
CXXFLAGS = $(CXXINCLUDEDIR) -D_FILE_OFFSET_BITS=64 -Wall -ggdb -g3 -fno-inline-small-functions -O0 -fno-inline -fno-eliminate-unused-debug-types
LIBS = -lz
RM = rm -f

OBJ  = fasta_stats.o

HEADER = kseq.h

BIN  = fasta_stats

all: $(BIN)


$(BIN): $(OBJ)
	$(CPP) $(CXXFLAGS) $(OBJ) -o $@ $(LIBS)

$(OBJ): $(HEADER)

clean: 
	$(RM) $(OBJ) $(BIN)

.PHONY: all all-before all-after clean clean-custom
