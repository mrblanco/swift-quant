#!/usr/bin/env python3
"""
DASH Depletion Analysis Tool - Command Line Interface

A comprehensive suite for analyzing DASH depletion reactions and
evaluating gRNA design performance.
"""

import argparse
import sys
from pathlib import Path

from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed
from dash_visualizer import (DASHVisualizer, export_metrics_to_csv,
                             export_metrics_to_excel)
from dash_statistics import (DASHStatistics, export_comparison_results)


def analyze_command(args):
    """Execute the analyze subcommand."""
    print(f"Analyzing DASH depletion for targets in {args.targets}")
    print(f"Using BAM file: {args.bam}")

    # Load targets
    targets = load_grna_targets_from_bed(args.targets)
    print(f"Loaded {len(targets)} gRNA targets")

    # Analyze
    with DASHAnalyzer(args.bam, args.control) as analyzer:
        metrics_list = analyzer.analyze_multiple_targets(targets)

    # Export results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "depletion_metrics.csv"
    export_metrics_to_csv(metrics_list, str(csv_path))
    print(f"\nMetrics exported to: {csv_path}")

    if args.excel:
        excel_path = output_dir / "depletion_metrics.xlsx"
        export_metrics_to_excel(metrics_list, str(excel_path))
        print(f"Excel report exported to: {excel_path}")

    # Generate visualizations if requested
    if args.plot:
        visualizer = DASHVisualizer()

        # Metrics comparison plot
        comparison_plot = output_dir / "metrics_comparison.png"
        visualizer.plot_metrics_comparison(metrics_list, str(comparison_plot))
        print(f"Metrics comparison plot: {comparison_plot}")

        # Heatmap
        heatmap_plot = output_dir / "depletion_heatmap.png"
        visualizer.plot_depletion_heatmap(metrics_list, str(heatmap_plot))
        print(f"Depletion heatmap: {heatmap_plot}")

        # Individual coverage plots if requested
        if args.coverage_plots:
            with DASHAnalyzer(args.bam, args.control) as analyzer:
                coverage_dir = output_dir / "coverage_plots"
                coverage_dir.mkdir(exist_ok=True)

                for target in targets[:min(len(targets), 20)]:  # Limit to 20
                    plot_path = coverage_dir / f"{target.name}_coverage.png"
                    visualizer.plot_coverage_comparison(analyzer, target, str(plot_path))

                print(f"Coverage plots saved to: {coverage_dir}")

    # Generate statistics
    if args.stats:
        stats = DASHStatistics()

        # Summary statistics
        summary = stats.generate_summary_statistics(metrics_list)
        summary_path = output_dir / "summary_statistics.csv"
        summary.to_csv(summary_path)
        print(f"\nSummary statistics: {summary_path}")

        # Correlation analysis
        correlation = stats.correlation_analysis(metrics_list)
        corr_path = output_dir / "correlation_matrix.csv"
        correlation.to_csv(corr_path)
        print(f"Correlation matrix: {corr_path}")

        # Rank gRNAs
        rankings = stats.rank_grnas(metrics_list)
        rank_path = output_dir / "grna_rankings.csv"
        rankings.to_csv(rank_path, index=False)
        print(f"gRNA rankings: {rank_path}")

        # Identify outliers
        outliers = stats.identify_outliers(metrics_list)
        print(f"\nHigh performers: {', '.join(outliers['high_performers']) if outliers['high_performers'] else 'None'}")
        print(f"Low performers: {', '.join(outliers['low_performers']) if outliers['low_performers'] else 'None'}")

    print("\nAnalysis complete!")


def compare_command(args):
    """Execute the compare subcommand."""
    print(f"Comparing two groups of gRNA experiments")
    print(f"Group 1: {args.targets1} with {args.bam1}")
    print(f"Group 2: {args.targets2} with {args.bam2}")

    # Load targets
    targets1 = load_grna_targets_from_bed(args.targets1)
    targets2 = load_grna_targets_from_bed(args.targets2)

    # Analyze both groups
    with DASHAnalyzer(args.bam1, args.control1) as analyzer1:
        metrics1 = analyzer1.analyze_multiple_targets(targets1)

    with DASHAnalyzer(args.bam2, args.control2) as analyzer2:
        metrics2 = analyzer2.analyze_multiple_targets(targets2)

    # Statistical comparison
    stats = DASHStatistics()
    comparison_results = stats.compare_groups(metrics1, metrics2, alpha=args.alpha)

    # Export results
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    results_path = output_dir / "comparison_results.csv"
    export_comparison_results(comparison_results, str(results_path))
    print(f"\nComparison results saved to: {results_path}")

    # Print significant differences
    print("\nStatistical Comparison Results:")
    print("=" * 80)
    for result in comparison_results:
        sig_marker = "***" if result.significant else ""
        print(f"\n{result.metric_name} {sig_marker}")
        print(f"  Group 1: {result.group1_mean:.4f} ± {result.group1_std:.4f}")
        print(f"  Group 2: {result.group2_mean:.4f} ± {result.group2_std:.4f}")
        print(f"  t-statistic: {result.t_statistic:.4f}, p-value: {result.p_value:.4f}")
        print(f"  Effect size (Cohen's d): {result.effect_size:.4f}")

    print("\n*** indicates statistically significant difference")


def visualize_command(args):
    """Execute the visualize subcommand."""
    print(f"Creating visualizations for {args.targets}")

    # Load targets
    targets = load_grna_targets_from_bed(args.targets)

    # Create visualizer
    visualizer = DASHVisualizer()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate coverage plots
    with DASHAnalyzer(args.bam, args.control) as analyzer:
        if args.all_in_one:
            # Multiple targets in one figure
            plot_path = output_dir / "all_targets_coverage.png"
            visualizer.plot_multiple_targets_coverage(analyzer, targets, str(plot_path))
            print(f"Combined coverage plot: {plot_path}")
        else:
            # Individual plots
            for target in targets:
                plot_path = output_dir / f"{target.name}_coverage.png"
                visualizer.plot_coverage_comparison(analyzer, target, str(plot_path))

            print(f"Individual coverage plots saved to: {output_dir}")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="DASH Depletion Analysis Tool Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze gRNA performance
  %(prog)s analyze --bam treated.bam --control untreated.bam --targets grnas.bed --output results/

  # Compare two experimental conditions
  %(prog)s compare --bam1 condition1.bam --bam2 condition2.bam \\
                   --targets1 grnas1.bed --targets2 grnas2.bed --output comparison/

  # Generate coverage visualizations
  %(prog)s visualize --bam treated.bam --targets grnas.bed --output plots/
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Analyze subcommand
    analyze_parser = subparsers.add_parser('analyze',
                                           help='Analyze DASH depletion performance')
    analyze_parser.add_argument('--bam', required=True,
                               help='BAM file with treated sample')
    analyze_parser.add_argument('--control',
                               help='BAM file with control/untreated sample')
    analyze_parser.add_argument('--targets', required=True,
                               help='BED file with gRNA target regions')
    analyze_parser.add_argument('--output', '-o', required=True,
                               help='Output directory for results')
    analyze_parser.add_argument('--plot', action='store_true',
                               help='Generate visualization plots')
    analyze_parser.add_argument('--coverage-plots', action='store_true',
                               help='Generate individual coverage plots for each target')
    analyze_parser.add_argument('--stats', action='store_true',
                               help='Generate statistical summaries')
    analyze_parser.add_argument('--excel', action='store_true',
                               help='Export results to Excel format')

    # Compare subcommand
    compare_parser = subparsers.add_parser('compare',
                                          help='Compare two groups of experiments')
    compare_parser.add_argument('--bam1', required=True,
                               help='BAM file for group 1')
    compare_parser.add_argument('--bam2', required=True,
                               help='BAM file for group 2')
    compare_parser.add_argument('--control1',
                               help='Control BAM file for group 1')
    compare_parser.add_argument('--control2',
                               help='Control BAM file for group 2')
    compare_parser.add_argument('--targets1', required=True,
                               help='BED file with targets for group 1')
    compare_parser.add_argument('--targets2', required=True,
                               help='BED file with targets for group 2')
    compare_parser.add_argument('--output', '-o', required=True,
                               help='Output directory for results')
    compare_parser.add_argument('--alpha', type=float, default=0.05,
                               help='Significance level for statistical tests (default: 0.05)')

    # Visualize subcommand
    visualize_parser = subparsers.add_parser('visualize',
                                            help='Generate coverage visualizations')
    visualize_parser.add_argument('--bam', required=True,
                                 help='BAM file to visualize')
    visualize_parser.add_argument('--control',
                                 help='Control BAM file for comparison')
    visualize_parser.add_argument('--targets', required=True,
                                 help='BED file with target regions')
    visualize_parser.add_argument('--output', '-o', required=True,
                                 help='Output directory for plots')
    visualize_parser.add_argument('--all-in-one', action='store_true',
                                 help='Create one figure with all targets')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        if args.command == 'analyze':
            analyze_command(args)
        elif args.command == 'compare':
            compare_command(args)
        elif args.command == 'visualize':
            visualize_command(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
