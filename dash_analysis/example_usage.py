#!/usr/bin/env python3
"""
Example usage of the DASH Analysis Tool Suite

This script demonstrates how to use the Python API for analyzing
DASH depletion experiments.
"""

from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed, GRNATarget
from dash_visualizer import DASHVisualizer, export_metrics_to_csv, export_metrics_to_excel
from dash_statistics import DASHStatistics, export_comparison_results


def example_basic_analysis():
    """Basic analysis example."""
    print("=" * 60)
    print("Example 1: Basic gRNA Performance Analysis")
    print("=" * 60)

    # Load gRNA targets from BED file
    targets = load_grna_targets_from_bed('example_targets.bed')
    print(f"\nLoaded {len(targets)} gRNA targets")

    # Analyze with treated BAM file (and optional control)
    with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
        metrics_list = analyzer.analyze_multiple_targets(targets)

    # Print results
    print("\nResults:")
    print("-" * 60)
    for metrics in metrics_list:
        print(f"\n{metrics.grna_name}:")
        print(f"  Depletion Efficiency: {metrics.depletion_efficiency:.3f}")
        print(f"  Coverage Uniformity:  {metrics.coverage_uniformity:.3f}")
        print(f"  Zero Coverage:        {metrics.zero_coverage_fraction:.3f}")
        print(f"  Mean Coverage:        {metrics.mean_coverage:.1f}")

    # Export to CSV
    export_metrics_to_csv(metrics_list, 'results.csv')
    print("\nResults exported to results.csv")


def example_with_visualizations():
    """Analysis with visualizations."""
    print("\n" + "=" * 60)
    print("Example 2: Analysis with Visualizations")
    print("=" * 60)

    targets = load_grna_targets_from_bed('example_targets.bed')
    visualizer = DASHVisualizer()

    with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
        # Analyze metrics
        metrics_list = analyzer.analyze_multiple_targets(targets)

        # Create coverage plot for first target
        visualizer.plot_coverage_comparison(
            analyzer,
            targets[0],
            'coverage_plot.png'
        )
        print("\nCoverage plot saved to coverage_plot.png")

    # Create performance heatmap
    visualizer.plot_depletion_heatmap(metrics_list, 'heatmap.png')
    print("Performance heatmap saved to heatmap.png")

    # Create metrics comparison plots
    visualizer.plot_metrics_comparison(metrics_list, 'comparison.png')
    print("Metrics comparison saved to comparison.png")


def example_statistical_analysis():
    """Statistical analysis example."""
    print("\n" + "=" * 60)
    print("Example 3: Statistical Analysis and Ranking")
    print("=" * 60)

    targets = load_grna_targets_from_bed('example_targets.bed')

    with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
        metrics_list = analyzer.analyze_multiple_targets(targets)

    stats = DASHStatistics()

    # Rank gRNAs by performance
    rankings = stats.rank_grnas(metrics_list)
    print("\nTop 3 Performing gRNAs:")
    print(rankings.head(3).to_string())

    # Identify outliers
    outliers = stats.identify_outliers(metrics_list)
    print(f"\nHigh Performers: {outliers['high_performers']}")
    print(f"Low Performers:  {outliers['low_performers']}")

    # Generate summary statistics
    summary = stats.generate_summary_statistics(metrics_list)
    print("\nSummary Statistics:")
    print(summary.to_string())

    # Correlation analysis
    correlation = stats.correlation_analysis(metrics_list)
    print("\nMetric Correlations:")
    print(correlation.to_string())


def example_compare_conditions():
    """Compare two experimental conditions."""
    print("\n" + "=" * 60)
    print("Example 4: Comparing Two Conditions")
    print("=" * 60)

    # Load targets for both conditions
    targets1 = load_grna_targets_from_bed('condition1_targets.bed')
    targets2 = load_grna_targets_from_bed('condition2_targets.bed')

    # Analyze condition 1
    with DASHAnalyzer('condition1_treated.bam', 'condition1_control.bam') as analyzer1:
        metrics1 = analyzer1.analyze_multiple_targets(targets1)

    # Analyze condition 2
    with DASHAnalyzer('condition2_treated.bam', 'condition2_control.bam') as analyzer2:
        metrics2 = analyzer2.analyze_multiple_targets(targets2)

    # Statistical comparison
    stats = DASHStatistics()
    comparison_results = stats.compare_groups(metrics1, metrics2, alpha=0.05)

    print("\nStatistical Comparison:")
    print("-" * 60)
    for result in comparison_results:
        sig = "***" if result.significant else ""
        print(f"\n{result.metric_name} {sig}")
        print(f"  Condition 1: {result.group1_mean:.4f} ± {result.group1_std:.4f}")
        print(f"  Condition 2: {result.group2_mean:.4f} ± {result.group2_std:.4f}")
        print(f"  p-value:     {result.p_value:.4f}")

    # Export comparison results
    export_comparison_results(comparison_results, 'comparison_stats.csv')
    print("\nComparison results saved to comparison_stats.csv")


def example_custom_targets():
    """Create targets programmatically."""
    print("\n" + "=" * 60)
    print("Example 5: Creating Targets Programmatically")
    print("=" * 60)

    # Create GRNATarget objects manually
    targets = [
        GRNATarget(
            name="custom_gRNA_1",
            chrom="chr1",
            start=100000,
            end=100500,
            strand="+"
        ),
        GRNATarget(
            name="custom_gRNA_2",
            chrom="chr2",
            start=200000,
            end=200500,
            strand="-"
        )
    ]

    with DASHAnalyzer('treated.bam') as analyzer:
        for target in targets:
            metrics = analyzer.analyze_grna_target(target)
            print(f"\n{metrics.grna_name}:")
            print(f"  Depletion: {metrics.depletion_efficiency:.3f}")
            print(f"  Uniformity: {metrics.coverage_uniformity:.3f}")


def example_batch_processing():
    """Process multiple BAM files."""
    print("\n" + "=" * 60)
    print("Example 6: Batch Processing Multiple Samples")
    print("=" * 60)

    from pathlib import Path

    targets = load_grna_targets_from_bed('example_targets.bed')
    bam_files = ['sample1.bam', 'sample2.bam', 'sample3.bam']

    all_results = {}

    for bam_file in bam_files:
        print(f"\nProcessing {bam_file}...")

        with DASHAnalyzer(bam_file) as analyzer:
            metrics_list = analyzer.analyze_multiple_targets(targets)

        # Store results
        all_results[bam_file] = metrics_list

        # Export individual results
        output_file = f"results_{Path(bam_file).stem}.csv"
        export_metrics_to_csv(metrics_list, output_file)
        print(f"  Results saved to {output_file}")


if __name__ == '__main__':
    print("\nDASH Analysis Tool Suite - Example Usage\n")

    # Run examples (comment out examples that require specific files)

    # Example 1: Basic analysis
    # example_basic_analysis()

    # Example 2: With visualizations
    # example_with_visualizations()

    # Example 3: Statistical analysis
    # example_statistical_analysis()

    # Example 4: Compare conditions
    # example_compare_conditions()

    # Example 5: Custom targets
    # example_custom_targets()

    # Example 6: Batch processing
    # example_batch_processing()

    print("\nUncomment the examples you want to run in the script.")
    print("Make sure you have the appropriate BAM and BED files available.")
