# Needleman–Wunsch Sequence Alignment (C++17, Multithreaded)

This project is a **C++17 implementation of the Needleman–Wunsch algorithm** for global sequence alignment.  
It demonstrates both the **dynamic programming (DP) matrix filling** and a **multithreaded traceback** to recover optimal alignments.

The program is designed for **educational and experimental use**.  
It shows how parallelism can be applied to traceback enumeration, but is **not optimized** for very large genomes (megabase scale).

---

## Features

- Global alignment via the **Needleman–Wunsch algorithm**.
- Configurable scoring:
  - Match reward +1
  - Mismatch penalty -1
  - Gap penalty -1
- DP matrix with **pointer bitmasks** to record multiple equally optimal paths. This is used if we needed multiple optimal alignments.
- **Multithreaded traceback** using `std::async`, exploring different paths in parallel (Only used in tracebacks).
- Support for reading sequences from files (FASTA-like: plain text, single line per sequence).
- Output of one or all optimal alignments.

---

## Building

Requires:
- A C++17 compiler (tested with `g++` 9+ or `clang++` 10+).
- POSIX threads (`-pthread`).

Compile with:

```bash
g++ -std=c++17 -O2 -pthread -o align align.cpp


```bash
#usage
./align <seq1_file> <seq2_file> [options]
```

## Options

--one-alignment
Return a single optimal deterministic alignment (default).

--score-only
Only give you the alignment score without the actual alignment

--max-alignment <int>
Controls how many alignments the program will actually return, default is 10000

--max-tasks
This controls how many asynchronous threads the program can spawn during traceback

--file1

--file2

### simply run ./align if you want to see the options

## Example usage (files in this repository)
```bash
./align --one-alignment ms2.fasta phiX174.fasta
```

## How It Works
1. DP Recurrence

The algorithm fills a DP matrix of size (len(seq1)+1) x (len(seq2)+1).
Each cell stores the best score from:

Diagonal (match/mismatch)

Up (gap in sequence 2)

Left (gap in sequence 1)


2. Pointer Bitmask

Instead of storing only the best score, we also store which direction(s) led to it:

Diagonal = 1

Up = 2

Left = 4

Just in case we needed all optimal alignments

3. Multithreaded Traceback

Traceback starts from the bottom-right of the matrix.

If multiple paths exist (bitmask > 1), each path may be explored in parallel with std::async.

This speeds up enumeration of many alignments, though in the worst case (many ties) it can still explode combinatorially.

4. Async Throttling

To prevent spawning unbounded threads:

A global std::atomic<int> counts active tasks.

If too many tasks are running, traceback continues synchronously in the current thread.

## Performance Notes

- Suitable for sequences up to a few thousand nucleotides.
- Not intended for whole-genome scale (megabases).
- Enumerating all alignments grows exponentially; in practice, use --one-alignment.
