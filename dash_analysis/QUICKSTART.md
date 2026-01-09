# DASH Analysis Quick Start Guide

This guide will help you get started with analyzing your DASH depletion experiments.

## Step 1: Installation

```bash
cd dash_analysis
pip install -r requirements.txt
```

## Step 2: Prepare Your Input Files

### Create a BED file with your gRNA target regions

Format (tab-separated):
```
chromosome  start   end     gRNA_name   score   strand
```

Example (`my_grnas.bed`):
```
chr1    10000   10500   gRNA_rRNA_1     .   +
chr1    15000   15500   gRNA_rRNA_2     .   +
chr12   6900000 6900500 gRNA_5S_1       .   +
```

### Prepare your BAM files

Your BAM files must be sorted and indexed:

```bash
# Sort BAM file
samtools sort your_file.bam -o treated.bam

# Index BAM file
samtools index treated.bam

# Do the same for control if you have one
samtools sort control_file.bam -o control.bam
samtools index control.bam
```

## Step 3: Run Your First Analysis

### Option A: Using the Command Line Interface

```bash
python dash_cli.py analyze \
    --bam treated.bam \
    --control control.bam \
    --targets my_grnas.bed \
    --output results/ \
    --plot \
    --stats \
    --excel
```

This will create:
- `results/depletion_metrics.csv` - Main results table
- `results/depletion_metrics.xlsx` - Excel formatted results
- `results/metrics_comparison.png` - Performance comparison plots
- `results/depletion_heatmap.png` - Heatmap of all metrics
- `results/summary_statistics.csv` - Statistical summary
- `results/grna_rankings.csv` - Ranked gRNAs

### Option B: Using Python

Create a script `analyze_my_data.py`:

```python
from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed
from dash_visualizer import DASHVisualizer, export_metrics_to_csv

# Load your gRNA targets
targets = load_grna_targets_from_bed('my_grnas.bed')

# Analyze
with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
    metrics = analyzer.analyze_multiple_targets(targets)

# Print results
for m in metrics:
    print(f"{m.grna_name}: Depletion={m.depletion_efficiency:.3f}")

# Save results
export_metrics_to_csv(metrics, 'my_results.csv')

# Create visualizations
viz = DASHVisualizer()
viz.plot_metrics_comparison(metrics, 'my_plots.png')
```

Run it:
```bash
python analyze_my_data.py
```

## Step 4: Interpret Your Results

### Understanding the Metrics

Open `results/depletion_metrics.csv` and look at:

1. **Depletion Efficiency** (0-1, higher is better)
   - 0.8+ : Excellent depletion
   - 0.5-0.8 : Good depletion
   - <0.5 : Poor depletion - consider redesigning gRNA

2. **Coverage Uniformity** (0-1, higher is better)
   - >0.7 : Uniform depletion across target
   - 0.5-0.7 : Moderate uniformity
   - <0.5 : Uneven depletion - check for secondary structures

3. **Zero Coverage Fraction** (0-1, higher is better)
   - >0.8 : Most of target region depleted
   - 0.5-0.8 : Partial depletion
   - <0.5 : Limited depletion

### Check the Rankings

Open `results/grna_rankings.csv` to see which gRNAs performed best/worst.

### Look at the Visualizations

1. **metrics_comparison.png** - Shows performance across all gRNAs
2. **depletion_heatmap.png** - Color-coded performance matrix

## Step 5: Identify Problem gRNAs

### Using Python:

```python
from dash_statistics import DASHStatistics

stats = DASHStatistics()
outliers = stats.identify_outliers(metrics)

print("Low performers:", outliers['low_performers'])
```

### What to do with poor performers:

1. **Low Depletion Efficiency**:
   - Check gRNA sequence matches reference genome
   - Look for off-target binding sites
   - Try alternative target location

2. **Poor Uniformity**:
   - Check GC content (aim for 40-60%)
   - Avoid strong secondary structures
   - Ensure target region is accessible

3. **Edge Effects** (coverage drops at edges):
   - Extend target region by 50-100bp
   - Use overlapping gRNAs

## Step 6: Compare Different Conditions (Optional)

If you want to compare two experimental conditions:

```bash
python dash_cli.py compare \
    --bam1 condition1.bam \
    --bam2 condition2.bam \
    --targets1 grnas1.bed \
    --targets2 grnas2.bed \
    --output comparison/
```

This will tell you which metrics are significantly different between conditions.

## Common Issues and Solutions

### "No coverage found"
- Check that chromosome names match between BAM and BED files
  - BAM uses "chr1" vs BED uses "1"? Modify one to match
- Verify coordinates are correct

### "BAM file not indexed"
```bash
samtools index your_file.bam
```

### "Memory error"
- Process fewer targets at a time
- Use smaller regions

### "Permission denied"
- Make sure output directory exists or script has write permissions

## Next Steps

1. Review the full [README.md](README.md) for advanced features
2. Look at [example_usage.py](example_usage.py) for more examples
3. Customize metric weights for ranking:

```python
from dash_statistics import DASHStatistics

custom_weights = {
    'depletion': 0.5,      # 50% weight
    'uniformity': 0.3,     # 30% weight
    'zero_coverage': 0.15, # 15% weight
    'cv': 0.05             # 5% weight
}

rankings = DASHStatistics.rank_grnas(metrics, weights=custom_weights)
```

## Getting Help

- Check the [README.md](README.md) for detailed documentation
- Look at [example_usage.py](example_usage.py) for code examples
- Review error messages carefully - they usually indicate the issue

## Quick Reference

### Common Commands

```bash
# Basic analysis
python dash_cli.py analyze --bam TREATED.bam --targets GRNAS.bed --output results/

# With control and plots
python dash_cli.py analyze --bam TREATED.bam --control CONTROL.bam --targets GRNAS.bed --output results/ --plot --stats

# Visualize only
python dash_cli.py visualize --bam TREATED.bam --targets GRNAS.bed --output plots/

# Compare conditions
python dash_cli.py compare --bam1 A.bam --bam2 B.bam --targets1 A.bed --targets2 B.bed --output comp/
```

Good luck with your DASH analysis!
