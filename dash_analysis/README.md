# DASH Depletion Analysis Tool Suite

A comprehensive suite of tools for analyzing DASH (DNA-based Adaptive Screening with CRISPR/Cas) depletion reactions to evaluate and improve gRNA design.

## Features

- **Coverage Analysis**: Calculate per-base coverage across gRNA target regions
- **Depletion Metrics**: Quantify depletion efficiency, coverage uniformity, and completeness
- **Statistical Comparison**: Compare performance across different experimental conditions
- **Visualizations**: Generate plots, heatmaps, and coverage profiles
- **gRNA Ranking**: Identify best and worst performing gRNAs
- **Export Options**: CSV, Excel, and PDF reports

## Installation

```bash
pip install -r requirements.txt
```

## Quick Start

### 1. Analyze gRNA Performance

```bash
python dash_cli.py analyze \
    --bam treated.bam \
    --control untreated.bam \
    --targets grna_targets.bed \
    --output results/ \
    --plot \
    --stats \
    --excel
```

### 2. Compare Two Conditions

```bash
python dash_cli.py compare \
    --bam1 condition1.bam \
    --bam2 condition2.bam \
    --targets1 grnas_set1.bed \
    --targets2 grnas_set2.bed \
    --output comparison/ \
    --alpha 0.05
```

### 3. Generate Visualizations

```bash
python dash_cli.py visualize \
    --bam treated.bam \
    --control untreated.bam \
    --targets grnas.bed \
    --output plots/
```

## Input File Formats

### BED File Format (gRNA targets)

Tab-separated file with at least 4 columns:

```
chr1    1000    1500    gRNA_001    .    +
chr1    2000    2500    gRNA_002    .    +
chr2    3000    3500    gRNA_003    .    -
```

Columns:
1. Chromosome name
2. Start position (0-based)
3. End position (exclusive)
4. gRNA name
5. Score (optional, use '.' if not applicable)
6. Strand (optional, default: '+')

### BAM Files

Standard BAM format alignment files. Must be sorted and indexed.

```bash
samtools sort input.bam -o sorted.bam
samtools index sorted.bam
```

## Output Files

### Analysis Command

- `depletion_metrics.csv`: Core metrics for each gRNA
- `depletion_metrics.xlsx`: Excel format with formatted tables
- `metrics_comparison.png`: Bar plots comparing key metrics
- `depletion_heatmap.png`: Heatmap of all performance metrics
- `coverage_plots/`: Individual coverage plots for each target
- `summary_statistics.csv`: Descriptive statistics
- `correlation_matrix.csv`: Correlation between metrics
- `grna_rankings.csv`: gRNAs ranked by composite performance score

### Comparison Command

- `comparison_results.csv`: Statistical test results comparing two groups

## Metrics Explained

### Depletion Efficiency
- **Definition**: Reduction in coverage compared to control (0-1, higher = better)
- **Formula**: 1 - (treated_coverage / control_coverage)
- **Interpretation**: 0.8+ indicates excellent depletion, <0.5 may need optimization

### Coverage Uniformity
- **Definition**: How evenly coverage is distributed across the target (0-1, higher = better)
- **Formula**: 1 / (1 + coefficient_of_variation)
- **Interpretation**: >0.7 indicates uniform depletion, <0.5 suggests uneven depletion

### Zero Coverage Fraction
- **Definition**: Proportion of bases with zero coverage (0-1)
- **Interpretation**: Higher values indicate more complete depletion

### Coefficient of Variation
- **Definition**: Standard deviation / mean coverage
- **Interpretation**: Lower values indicate more consistent coverage

## Python API Usage

### Basic Analysis

```python
from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed

# Load targets
targets = load_grna_targets_from_bed('grnas.bed')

# Analyze
with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
    metrics_list = analyzer.analyze_multiple_targets(targets)

    for metrics in metrics_list:
        print(f"{metrics.grna_name}: {metrics.depletion_efficiency:.3f}")
```

### Visualization

```python
from dash_visualizer import DASHVisualizer, export_metrics_to_csv

visualizer = DASHVisualizer()

# Create plots
with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
    visualizer.plot_coverage_comparison(analyzer, targets[0], 'coverage.png')
    visualizer.plot_depletion_heatmap(metrics_list, 'heatmap.png')

# Export data
export_metrics_to_csv(metrics_list, 'results.csv')
```

### Statistical Analysis

```python
from dash_statistics import DASHStatistics

stats = DASHStatistics()

# Rank gRNAs
rankings = stats.rank_grnas(metrics_list)
print(rankings.head())

# Identify outliers
outliers = stats.identify_outliers(metrics_list)
print(f"High performers: {outliers['high_performers']}")
print(f"Low performers: {outliers['low_performers']}")

# Compare groups
comparison = stats.compare_groups(group1_metrics, group2_metrics)
for result in comparison:
    if result.significant:
        print(f"{result.metric_name}: p={result.p_value:.4f}")
```

## Interpreting Results

### Identifying Poor Performers

Look for gRNAs with:
- Low depletion efficiency (<0.5)
- Low coverage uniformity (<0.5)
- High coefficient of variation (>1.0)

### Design Recommendations

Based on analysis results:

1. **Low Depletion Efficiency**:
   - Check for off-target binding
   - Verify gRNA sequence matches target
   - Consider alternative gRNA locations

2. **Poor Uniformity**:
   - May indicate secondary structure issues
   - Check GC content (optimal: 40-60%)
   - Consider target accessibility

3. **Edge Effects**:
   - Coverage drops at target boundaries
   - Extend target region by 50-100bp
   - Use multiple overlapping gRNAs

## Advanced Usage

### Custom Scoring Weights

```python
from dash_statistics import DASHStatistics

# Customize metric weights
custom_weights = {
    'depletion': 0.5,      # Prioritize depletion efficiency
    'uniformity': 0.3,
    'zero_coverage': 0.15,
    'cv': 0.05
}

rankings = DASHStatistics.rank_grnas(metrics_list, weights=custom_weights)
```

### Batch Processing

```python
from pathlib import Path
from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed
from dash_visualizer import export_metrics_to_csv

bam_files = Path('data/').glob('*.bam')

for bam_file in bam_files:
    targets = load_grna_targets_from_bed('targets.bed')

    with DASHAnalyzer(str(bam_file)) as analyzer:
        metrics = analyzer.analyze_multiple_targets(targets)

    output_path = f"results/{bam_file.stem}_metrics.csv"
    export_metrics_to_csv(metrics, output_path)
```

## Troubleshooting

### Common Issues

**"No coverage found"**
- Verify BAM file is sorted and indexed
- Check chromosome naming matches between BED and BAM files
- Ensure coordinates are valid

**"Memory error with large BAM files"**
- Process targets in smaller batches
- Use `pysam.AlignmentFile(..., threads=4)` for parallel decompression

**"Plots not displaying"**
- Use `--plot` flag or specify output paths
- Check matplotlib backend configuration

## License

MIT License

## Contributing

Contributions welcome! Please open an issue or submit a pull request.
