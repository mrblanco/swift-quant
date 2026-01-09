# DASH Analysis Tool Suite - Project Summary

## Overview

A complete Python-based analysis suite for DASH (DNA-based Adaptive Screening with CRISPR/Cas) depletion experiments. This toolkit helps identify areas for improvement in gRNA design by analyzing BAM files and quantifying depletion performance.

## What Was Built

### Core Modules

1. **[dash_analyzer.py](dash_analyzer.py)** - Core analysis engine
   - `DASHAnalyzer` class for BAM file processing
   - Per-base coverage calculation
   - Depletion efficiency metrics
   - Coverage uniformity quantification
   - `GRNATarget` dataclass for target regions
   - `DepletionMetrics` dataclass for results
   - BED file loader for gRNA targets

2. **[dash_visualizer.py](dash_visualizer.py)** - Visualization tools
   - `DASHVisualizer` class for creating plots
   - Coverage comparison plots (treated vs control)
   - Multi-target grid plots
   - Performance heatmaps
   - Metrics comparison bar charts
   - CSV/Excel export functions
   - PDF report generation

3. **[dash_statistics.py](dash_statistics.py)** - Statistical analysis
   - `DASHStatistics` class for statistical tests
   - Two-group comparison with t-tests
   - Cohen's d effect size calculation
   - gRNA ranking system with weighted scoring
   - Outlier detection (high/low performers)
   - Correlation analysis between metrics
   - Summary statistics generation

4. **[dash_design_helper.py](dash_design_helper.py)** - Design recommendations
   - `DASHDesignHelper` class for automated analysis
   - Issue detection and categorization
   - Severity-based recommendations (critical/warning/info)
   - Design improvement suggestions
   - Performance recognition
   - Comprehensive report generation
   - CSV export of recommendations

5. **[dash_cli.py](dash_cli.py)** - Command-line interface
   - Three main commands: `analyze`, `compare`, `visualize`
   - Full argument parsing and validation
   - Batch processing support
   - Multiple output formats
   - Progress reporting

### Documentation

1. **[README.md](README.md)** - Complete documentation
   - Feature overview
   - Installation instructions
   - Usage examples (CLI and Python API)
   - Metrics explanations
   - Troubleshooting guide
   - Advanced usage patterns

2. **[QUICKSTART.md](QUICKSTART.md)** - Getting started guide
   - Step-by-step setup
   - First analysis walkthrough
   - Results interpretation
   - Common issues and solutions
   - Quick reference commands

3. **[example_usage.py](example_usage.py)** - Code examples
   - 6 complete usage examples
   - Basic analysis
   - Visualization generation
   - Statistical analysis
   - Condition comparison
   - Custom target creation
   - Batch processing

4. **[example_targets.bed](example_targets.bed)** - Sample data
   - Example gRNA targets for rRNA depletion
   - Proper BED format demonstration

## Key Features

### Analysis Capabilities

- **Depletion Efficiency**: Quantifies reduction in coverage (treated vs control)
- **Coverage Uniformity**: Measures evenness of depletion across target regions
- **Zero Coverage Fraction**: Identifies completely depleted bases
- **Coefficient of Variation**: Detects inconsistent depletion patterns
- **Statistical Comparison**: Compare experimental conditions with significance testing
- **Performance Ranking**: Identify best and worst performing gRNAs

### Visualizations

- Coverage profiles (treated vs control overlay)
- Multi-target comparison grids
- Performance heatmaps
- Metric comparison bar charts
- PDF summary reports

### Design Recommendations

- Automated issue detection
- Severity-based prioritization
- Actionable improvement suggestions
- Performance benchmarking
- Design pattern recognition

## File Structure

```
dash_analysis/
├── dash_analyzer.py           # Core analysis engine
├── dash_visualizer.py         # Plotting and export
├── dash_statistics.py         # Statistical tools
├── dash_design_helper.py      # Design recommendations
├── dash_cli.py                # Command-line interface
├── example_usage.py           # Python API examples
├── example_targets.bed        # Sample gRNA targets
├── requirements.txt           # Python dependencies
├── README.md                  # Full documentation
├── QUICKSTART.md              # Getting started guide
└── PROJECT_SUMMARY.md         # This file
```

## Dependencies

- **pysam** >= 0.19.0 - BAM file handling
- **numpy** >= 1.21.0 - Numerical computations
- **pandas** >= 1.3.0 - Data manipulation
- **matplotlib** >= 3.4.0 - Plotting
- **seaborn** >= 0.11.0 - Statistical visualizations
- **scipy** >= 1.7.0 - Statistical tests
- **openpyxl** >= 3.0.0 - Excel export

## Usage Modes

### 1. Command-Line Interface

```bash
# Analyze with plots and statistics
python dash_cli.py analyze \
    --bam treated.bam \
    --control untreated.bam \
    --targets grnas.bed \
    --output results/ \
    --plot --stats --excel

# Compare two conditions
python dash_cli.py compare \
    --bam1 cond1.bam --bam2 cond2.bam \
    --targets1 grnas1.bed --targets2 grnas2.bed \
    --output comparison/

# Visualize coverage
python dash_cli.py visualize \
    --bam treated.bam \
    --targets grnas.bed \
    --output plots/
```

### 2. Python API

```python
from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed
from dash_visualizer import DASHVisualizer, export_metrics_to_csv
from dash_statistics import DASHStatistics
from dash_design_helper import DASHDesignHelper

# Load targets
targets = load_grna_targets_from_bed('grnas.bed')

# Analyze
with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
    metrics = analyzer.analyze_multiple_targets(targets)

# Visualize
viz = DASHVisualizer()
viz.plot_metrics_comparison(metrics, 'comparison.png')

# Get statistics
stats = DASHStatistics()
rankings = stats.rank_grnas(metrics)

# Get recommendations
helper = DASHDesignHelper()
report = helper.generate_report(metrics)
print(report)
```

## Output Files

### From `analyze` command:

- `depletion_metrics.csv` - Main results table
- `depletion_metrics.xlsx` - Excel formatted results
- `metrics_comparison.png` - Bar chart comparisons
- `depletion_heatmap.png` - Performance heatmap
- `coverage_plots/` - Individual coverage profiles
- `summary_statistics.csv` - Descriptive stats
- `correlation_matrix.csv` - Metric correlations
- `grna_rankings.csv` - Ranked performance

### From `compare` command:

- `comparison_results.csv` - Statistical test results

### From design helper:

- `design_recommendations.csv` - Improvement suggestions

## Metrics Explained

| Metric | Range | Interpretation |
|--------|-------|----------------|
| Depletion Efficiency | 0-1 | Higher = better depletion. 0.8+ excellent, <0.5 needs work |
| Coverage Uniformity | 0-1 | Higher = more even. >0.7 uniform, <0.5 uneven |
| Zero Coverage Fraction | 0-1 | Higher = more complete. >0.8 excellent coverage |
| Coefficient of Variation | 0+ | Lower = more consistent. <1.0 good, >2.0 problematic |

## Common Workflows

### 1. Initial Analysis
```bash
python dash_cli.py analyze --bam treated.bam --control control.bam \
    --targets grnas.bed --output results/ --plot --stats
```

### 2. Review Results
- Open `results/grna_rankings.csv` to see performance
- Check `results/depletion_heatmap.png` for visual overview
- Review `results/metrics_comparison.png` for detailed metrics

### 3. Get Design Recommendations
```python
from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed
from dash_design_helper import DASHDesignHelper

targets = load_grna_targets_from_bed('grnas.bed')
with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
    metrics = analyzer.analyze_multiple_targets(targets)

helper = DASHDesignHelper()
report = helper.generate_report(metrics)
print(report)
```

### 4. Redesign Poor Performers
- Identify low performers from rankings
- Read design recommendations
- Modify gRNA sequences/locations
- Re-run analysis

### 5. Compare Improvements
```bash
python dash_cli.py compare \
    --bam1 original.bam --bam2 improved.bam \
    --targets1 original_grnas.bed --targets2 improved_grnas.bed \
    --output comparison/
```

## Design Recommendations System

The tool automatically identifies:

- **Critical Issues** (Priority 5)
  - Very low depletion (<0.3)
  - Poor uniformity (<0.3)
  - Extreme variation (CV >2.0)

- **Warnings** (Priority 2-3)
  - Moderate depletion (0.3-0.6)
  - Moderate uniformity (0.3-0.5)
  - High variation (CV 1.0-2.0)

- **Info/Recognition**
  - Edge effects
  - Incomplete depletion
  - Excellent performance

Each recommendation includes:
- Specific issue description
- Root cause analysis
- Actionable improvement steps

## Next Steps

1. **Get Started**: Follow [QUICKSTART.md](QUICKSTART.md)
2. **Learn the API**: Review [example_usage.py](example_usage.py)
3. **Customize**: Modify metric weights, thresholds, and plotting styles
4. **Integrate**: Build into your analysis pipeline
5. **Extend**: Add new metrics or visualization types

## Future Enhancements (Ideas)

- Sequence-based features (GC content, secondary structure)
- Machine learning for gRNA performance prediction
- Integration with gRNA design tools
- Multi-sample batch comparison
- Interactive HTML reports
- Off-target prediction integration
- Automated gRNA redesign suggestions

## Support

- Review the [README.md](README.md) for detailed documentation
- Check [QUICKSTART.md](QUICKSTART.md) for getting started
- See [example_usage.py](example_usage.py) for code examples
- Examine error messages for debugging hints

## License

MIT License - Free to use and modify

---

**Built for analyzing DASH depletion reactions and improving gRNA design**
