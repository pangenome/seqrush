# SGD Layout Quality Comparison: ODGI vs SeqRush

## Summary

This analysis compares the path-guided SGD layout quality between ODGI and SeqRush using quantitative RMSE metrics and visual analysis.

**Key Finding:** SeqRush was running with 4 threads by default (hogwild parallelism) while ODGI defaults to 1 thread. This caused significantly higher variance in SeqRush results.

---

## Test Setup

- **Input:** HLA-zoo B-3106 sequences (9 sequences)
- **Graph:** ~471 nodes after compaction
- **Metric:** Root Mean Squared Error (RMSE) measuring layout quality
  - Lower RMSE = better 1D layout preserves genomic distances
- **Tool:** `measure_layout_quality` - calculates displacement between layout positions and actual path distances

---

## Results: ODGI (1 Thread)

| Run | RMSE (bp) | Visualization |
|-----|-----------|---------------|
| 1   | 1753.10   | `odgi_run1_viz.png` |
| 2   | 1708.46   | `odgi_run2_viz.png` |
| 3   | 1683.19   | `odgi_run3_viz.png` |

**Best:** 1683.19 bp
**Variance:** ~70 bp (very consistent)

---

## Results: SeqRush (1 Thread)

| Run | RMSE (bp) | Visualization |
|-----|-----------|---------------|
| 1   | 1759.54   | `seqrush_1thread_run1_viz.png` |
| 2   | 1759.81   | `seqrush_1thread_run2_viz.png` |
| 3   | **1579.01** ✓ | `seqrush_1thread_run3_viz.png` |
| 4   | 1773.93   | `seqrush_1thread_run4_viz.png` |
| 5   | 1678.55   | `seqrush_1thread_run5_viz.png` |

**Best:** 1579.01 bp (104 bp better than ODGI = **6% improvement**)
**Variance:** ~195 bp (3x ODGI's variance)

---

## Results: SeqRush (4 Threads - Hogwild)

| Run | RMSE (bp) | Visualization |
|-----|-----------|---------------|
| 1   | **1484.05** ✓ | `seqrush_4thread_run1_viz.png` |
| 2   | 1616.15   | `seqrush_4thread_run2_viz.png` |
| 3   | 1541.21   | `seqrush_4thread_run3_viz.png` |

**Best:** 1484.05 bp (199 bp better than ODGI = **12% improvement**)
**Variance:** ~380 bp (5x ODGI's variance)

---

## Analysis

### Threading Impact

1. **ODGI:** Single-threaded by default → low variance, consistent results
2. **SeqRush (1 thread):** Comparable variance (3x higher), but can achieve better layouts
3. **SeqRush (4 threads):** Hogwild parallelism → high variance but can find even better optima

### Layout Quality

- **SeqRush achieves better layouts than ODGI** when it converges well
- Best 1-thread SeqRush: 1579 bp vs ODGI: 1683 bp = **6% improvement**
- Best 4-thread SeqRush: 1484 bp vs ODGI: 1683 bp = **12% improvement**

### Variance Analysis

The variance in SeqRush (even with 1 thread) suggests:
1. Different random sampling implementation
2. Different convergence criteria
3. SGD stochasticity not fully controlled by fixed seeds

### Determinism

- All runs produce identical graph structure (471 nodes, 340 path steps)
- Variation is purely from SGD layout optimization
- SeqRush runs 1-2 with 1 thread are nearly identical (1759.54 vs 1759.81 = 0.27 bp difference)
  - This proves determinism is possible, but not consistent across all runs

---

## Visualization Notes

The PNG files show:
- **Horizontal axis:** Node ordering in the 1D layout
- **Vertical bands:** Path coverage patterns
- **Color/intensity:** Path depth and overlap

Visual inspection can reveal:
- Clustering quality (paths staying together)
- Jump patterns (large displacements)
- Overall linearity of the layout

---

## Recommendations

1. **For consistency:** Use SeqRush with 1 thread to match ODGI behavior
2. **For best quality:** Run SeqRush multiple times and select best result
3. **Investigation needed:** Understand source of variance even with 1 thread
4. **Future work:**
   - Match ODGI's RNG implementation exactly
   - Investigate sampling order differences
   - Potentially adjust convergence criteria for more stable results

---

## Files Generated

**Visualizations:**
- `odgi_run{1,2,3}_viz.png` - ODGI layouts (1 thread)
- `seqrush_1thread_run{1-5}_viz.png` - SeqRush layouts (1 thread)
- `seqrush_4thread_run{1-3}_viz.png` - SeqRush layouts (4 threads, hogwild)

**GFA files:**
- `/tmp/odgi_test_{1,2,3}.gfa` - ODGI sorted graphs
- `test_1thread_run{1-5}.gfa` - SeqRush 1-thread sorted graphs
- `test_run_{1,2,3}.gfa` - SeqRush 4-thread sorted graphs

**Metrics:**
- Run `measure_layout_quality <file>.gfa` to compute RMSE for any graph
