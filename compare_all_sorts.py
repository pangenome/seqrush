#!/usr/bin/env python3
"""
Compare edge distances across multiple sorting methods
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

def parse_gfa_edges(filename):
    """Parse L lines from GFA and extract edge distances."""
    edges = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('L'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    from_node = int(parts[1])
                    to_node = int(parts[3])
                    distance = abs(to_node - from_node)
                    edges.append(distance)
    return edges

def main():
    # Files to compare
    files = {
        'Unsorted': 'test_nosort.gfa',
        'Topological': 'test_normal.gfa',
        'Groom+Topo': 'test_groom_only.gfa',
        'SeqRush Ygs\n(SGD+Groom+Topo)': 'test_seqrush_ygs.gfa',
        'ODGI Ygs': 'test_odgi_ygs.gfa'
    }

    # Collect data
    data = {}
    for label, filename in files.items():
        try:
            distances = parse_gfa_edges(filename)
            if distances:
                data[label] = distances
                print(f"{label:15} - {len(distances)} edges")
        except FileNotFoundError:
            print(f"Warning: {filename} not found, skipping {label}")

    if not data:
        print("No data to plot")
        return

    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle('Edge Distance Comparison: SeqRush vs ODGI Sorting Methods', fontsize=16)

    # Box plot
    ax1 = axes[0, 0]
    bp = ax1.boxplot(list(data.values()), tick_labels=list(data.keys()), showfliers=False)
    ax1.set_ylabel('Edge Distance (Node ID Difference)')
    ax1.set_title('Edge Distance Distribution (outliers hidden)')
    ax1.grid(True, alpha=0.3)

    # Violin plot for better distribution visualization
    ax2 = axes[0, 1]
    positions = range(1, len(data) + 1)
    parts = ax2.violinplot(list(data.values()), positions=positions, showmeans=True)
    ax2.set_xticks(positions)
    ax2.set_xticklabels(list(data.keys()))
    ax2.set_ylabel('Edge Distance')
    ax2.set_title('Edge Distance Distribution (Violin Plot)')
    ax2.set_ylim(0, 50)  # Focus on local structure
    ax2.grid(True, alpha=0.3)

    # Statistics table
    ax3 = axes[1, 0]
    ax3.axis('tight')
    ax3.axis('off')

    # Calculate statistics
    stats_data = []
    for label in data.keys():
        distances = data[label]
        local = sum(1 for d in distances if d == 1) / len(distances) * 100
        nearby = sum(1 for d in distances if d <= 10) / len(distances) * 100
        stats_data.append([
            label,
            f"{np.mean(distances):.2f}",
            f"{np.median(distances):.1f}",
            f"{np.max(distances)}",
            f"{local:.1f}%",
            f"{nearby:.1f}%"
        ])

    table = ax3.table(cellText=stats_data,
                     colLabels=['Method', 'Mean', 'Median', 'Max', 'Local (d=1)', 'Nearby (d≤10)'],
                     cellLoc='center',
                     loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)

    # Color code the table
    for i in range(1, len(stats_data) + 1):
        # Color the mean column based on value
        mean_val = float(stats_data[i-1][1])
        if mean_val < 5:
            table[(i, 1)].set_facecolor('#90EE90')  # Light green for good
        elif mean_val < 10:
            table[(i, 1)].set_facecolor('#FFFFE0')  # Light yellow for okay
        else:
            table[(i, 1)].set_facecolor('#FFB6C1')  # Light red for poor

        # Color the local % column
        local_val = float(stats_data[i-1][4].rstrip('%'))
        if local_val > 92:
            table[(i, 4)].set_facecolor('#90EE90')
        elif local_val > 90:
            table[(i, 4)].set_facecolor('#FFFFE0')
        else:
            table[(i, 4)].set_facecolor('#FFB6C1')

    ax3.set_title('Statistics Summary', fontweight='bold')

    # Cumulative distribution
    ax4 = axes[1, 1]
    for label, distances in data.items():
        sorted_distances = np.sort(distances)
        cumulative = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances)
        ax4.plot(sorted_distances, cumulative, label=label, linewidth=2)

    ax4.set_xlabel('Edge Distance')
    ax4.set_ylabel('Cumulative Probability')
    ax4.set_title('Cumulative Distribution of Edge Distances')
    ax4.legend(loc='lower right')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim(0, 100)  # Focus on the interesting part

    plt.tight_layout()
    plt.savefig('all_methods_comparison.png', dpi=150)
    print("\nPlot saved to all_methods_comparison.png")

    # Print improvement summary
    print("\n=== IMPROVEMENT SUMMARY ===")
    unsorted_mean = np.mean(data['Unsorted']) if 'Unsorted' in data else 0
    for label in data.keys():
        if label != 'Unsorted' and unsorted_mean > 0:
            improvement = (unsorted_mean - np.mean(data[label])) / unsorted_mean * 100
            print(f"{label:15} - {improvement:.1f}% improvement over unsorted")

if __name__ == "__main__":
    main()