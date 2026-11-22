# ─────────────────────────────────────────────────────────────────────────────
# Makefile — CluMP Optimization Code
# Cluster-based Message-Passing (CluMP) Optimization in Rugged Energy Landscapes
# Rissone, Boettcher, Amendola, Sala, Ricci-Tersenghi (2025)
#
# Usage:
#   make          → build all binaries
#   make <target> → build a single binary (e.g. make CluMP)
#   make clean    → remove all compiled binaries
# ─────────────────────────────────────────────────────────────────────────────

# ── Compiler & flags ─────────────────────────────────────────────────────────
CC      = gcc
LIBS    = -lm

# Detect OS: use -std=c99 on macOS (Darwin), plain -O2 on Linux
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CFLAGS = -std=c99 -O2 -Wall
else
    CFLAGS = -O2 -Wall
endif

# ── Targets ──────────────────────────────────────────────────────────────────
TARGETS = GraphGenerator GetCluster CluMP RCluMP SimulatedAnnealing PopulationAnnealing ExtremalOpt

.PHONY: all clean

all: $(TARGETS)

GraphGenerator: GraphGenerator.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

GetCluster: GetCluster.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

CluMP: CluMP.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

RCluMP: RCluMP.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

SimulatedAnnealing: SimulatedAnnealing.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

PopulationAnnealing: PopulationAnnealing.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

ExtremalOpt: ExtremalOpt.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -f $(TARGETS)
