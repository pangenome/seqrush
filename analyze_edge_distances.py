#!/usr/bin/env python3
"""
Analyze edge distances in a GFA file to measure graph locality.
Edge distance = abs(from_node_id - to_node_id)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_gfa_edges(filename):
    """Parse L lines from GFA and extract edge distances."""
    edges = []
    node_degrees = defaultdict(int)

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('L'):
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    from_node = int(parts[1])
                    to_node = int(parts[3])

                    # Calculate distance in node ID space
                    distance = abs(to_node - from_node)
                    edges.append((from_node, to_node, distance))

                    # Track node degrees
                    node_degrees[from_node] += 1
                    node_degrees[to_node] += 1

    return edges, node_degrees

def analyze_distances(edges):
    """Compute statistics about edge distances."""
    if not edges:
        return None

    distances = [e[2] for e in edges]

    stats = {
        'count': len(distances),
        'mean': np.mean(distances),
        'median': np.median(distances),
        'std': np.std(distances),
        'min': np.min(distances),
        'max': np.max(distances),
        'p25': np.percentile(distances, 25),
        'p75': np.percentile(distances, 75),
        'p90': np.percentile(distances, 90),
        'p95': np.percentile(distances, 95),
        'p99': np.percentile(distances, 99),
    }

    # Count local edges (distance = 1)
    local_edges = sum(1 for d in distances if d == 1)
    stats['local_ratio'] = local_edges / len(distances) if distances else 0

    # Count "nearby" edges (distance <= 10)
    nearby_edges = sum(1 for d in distances if d <= 10)
    stats['nearby_ratio'] = nearby_edges / len(distances) if distances else 0

    return stats, distances

def plot_distance_distribution(distances, title, output_file=None):
    """Plot histogram of edge distances."""
    if not distances:
        print("No distances to plot")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Linear scale histogram
    ax1.hist(distances, bins=50, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Edge Distance (node ID difference)')
    ax1.set_ylabel('Count')
    ax1.set_title(f'{title} - Linear Scale')
    ax1.axvline(x=np.median(distances), color='red', linestyle='--', label=f'Median: {np.median(distances):.1f}')
    ax1.axvline(x=np.mean(distances), color='green', linestyle='--', label=f'Mean: {np.mean(distances):.1f}')
    ax1.legend()

    # Log scale histogram
    # Filter out zeros for log scale
    nonzero_distances = [d for d in distances if d > 0]
    if nonzero_distances:
        ax2.hist(nonzero_distances, bins=50, edgecolor='black', alpha=0.7)
        ax2.set_xlabel('Edge Distance (node ID difference)')
        ax2.set_ylabel('Count')
        ax2.set_title(f'{title} - Log Scale')
        ax2.set_yscale('log')
        ax2.axvline(x=np.median(distances), color='red', linestyle='--', label=f'Median: {np.median(distances):.1f}')
        ax2.axvline(x=np.mean(distances), color='green', linestyle='--', label=f'Mean: {np.mean(distances):.1f}')
        ax2.legend()

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=150)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

def print_stats(stats, title):
    """Print statistics in a formatted way."""
    print(f"\n=== {title} ===")
    print(f"Total edges: {stats['count']}")
    print(f"Mean distance: {stats['mean']:.2f}")
    print(f"Median distance: {stats['median']:.1f}")
    print(f"Std deviation: {stats['std']:.2f}")
    print(f"Min distance: {stats['min']}")
    print(f"Max distance: {stats['max']}")
    print(f"\nPercentiles:")
    print(f"  25th: {stats['p25']:.1f}")
    print(f"  50th (median): {stats['median']:.1f}")
    print(f"  75th: {stats['p75']:.1f}")
    print(f"  90th: {stats['p90']:.1f}")
    print(f"  95th: {stats['p95']:.1f}")
    print(f"  99th: {stats['p99']:.1f}")
    print(f"\nLocality metrics:")
    print(f"  Local edges (distance=1): {stats['local_ratio']*100:.1f}%")
    print(f"  Nearby edges (distance≤10): {stats['nearby_ratio']*100:.1f}%")

def compare_files(file1, file2):
    """Compare edge distances between two GFA files."""
    print(f"Comparing {file1} vs {file2}")

    # Analyze first file
    edges1, degrees1 = parse_gfa_edges(file1)
    stats1, distances1 = analyze_distances(edges1)

    # Analyze second file
    edges2, degrees2 = parse_gfa_edges(file2)
    stats2, distances2 = analyze_distances(edges2)

    # Print statistics
    print_stats(stats1, file1)
    print_stats(stats2, file2)

    # Print comparison
    print(f"\n=== Comparison ===")
    print(f"Mean distance reduction: {stats1['mean'] - stats2['mean']:.2f} ({(stats1['mean'] - stats2['mean'])/stats1['mean']*100:.1f}%)")
    print(f"Median distance reduction: {stats1['median'] - stats2['median']:.1f} ({(stats1['median'] - stats2['median'])/stats1['median']*100:.1f}%)")
    print(f"Max distance reduction: {stats1['max'] - stats2['max']:.0f}")
    print(f"Local edge improvement: {(stats2['local_ratio'] - stats1['local_ratio'])*100:.1f} percentage points")
    print(f"Nearby edge improvement: {(stats2['nearby_ratio'] - stats1['nearby_ratio'])*100:.1f} percentage points")

    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot both distributions
    axes[0,0].hist(distances1, bins=50, alpha=0.5, label=file1, edgecolor='black')
    axes[0,0].hist(distances2, bins=50, alpha=0.5, label=file2, edgecolor='black')
    axes[0,0].set_xlabel('Edge Distance')
    axes[0,0].set_ylabel('Count')
    axes[0,0].set_title('Edge Distance Distribution Comparison')
    axes[0,0].legend()

    # Log scale version
    axes[0,1].hist([d for d in distances1 if d > 0], bins=50, alpha=0.5, label=file1, edgecolor='black')
    axes[0,1].hist([d for d in distances2 if d > 0], bins=50, alpha=0.5, label=file2, edgecolor='black')
    axes[0,1].set_xlabel('Edge Distance')
    axes[0,1].set_ylabel('Count (log scale)')
    axes[0,1].set_yscale('log')
    axes[0,1].set_title('Edge Distance Distribution (Log Scale)')
    axes[0,1].legend()

    # Box plot comparison
    axes[1,0].boxplot([distances1, distances2], labels=[file1.split('/')[-1], file2.split('/')[-1]])
    axes[1,0].set_ylabel('Edge Distance')
    axes[1,0].set_title('Edge Distance Box Plot')

    # Cumulative distribution
    axes[1,1].hist(distances1, bins=100, cumulative=True, density=True, alpha=0.5, label=file1)
    axes[1,1].hist(distances2, bins=100, cumulative=True, density=True, alpha=0.5, label=file2)
    axes[1,1].set_xlabel('Edge Distance')
    axes[1,1].set_ylabel('Cumulative Probability')
    axes[1,1].set_title('Cumulative Distribution of Edge Distances')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('edge_distance_comparison.png', dpi=150)
    print("\nComparison plot saved to edge_distance_comparison.png")

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_edge_distances.py <gfa_file> [<gfa_file2>]")
        print("  Single file: analyze edge distances")
        print("  Two files: compare edge distances between files")
        sys.exit(1)

    file1 = sys.argv[1]

    if len(sys.argv) == 2:
        # Single file analysis
        edges, degrees = parse_gfa_edges(file1)
        stats, distances = analyze_distances(edges)

        if stats:
            print_stats(stats, file1)
            plot_distance_distribution(distances, file1, 'edge_distances.png')
        else:
            print("No edges found in file")

    else:
        # Compare two files
        file2 = sys.argv[2]
        compare_files(file1, file2)

if __name__ == "__main__":
    main()