# SeqRush

A dynamic, parallel, in-memory bioinformatics tool implementing biWFA+seqwish for pangenome graph construction.

## Features

- **Lockless Union-Find**: Uses UFRush for parallel graph construction
- **Dynamic Alignment**: Performs random pairwise alignments until convergence
- **Parallel Processing**: Multi-threaded alignment using Rayon
- **GFA Output**: Generates standard GFA format graphs
- **Configurable Parameters**: Alignment scoring, convergence criteria

## Usage

```bash
seqrush -s sequences.fasta -o output.gfa --threads 8 --verbose
```

## Algorithm

SeqRush implements the core concept from seqwish: using a union-find data structure to build pangenome graphs by iteratively aligning sequence pairs until the graph converges.

Key improvements:
- Lock-free parallel union-find operations
- Dynamic convergence detection
- Configurable alignment parameters
- Efficient GFA output with edge detection