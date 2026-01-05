# swift-quant

A Python tool for calculating the percentage of reads in BAM files that overlap with specific gene annotations.

## Overview

`swift-quant` provides efficient quantification of sequencing read overlap with genomic features. It supports both GTF and BED annotation formats and includes flexible filtering options for quality control.

## Features

- **Multiple annotation formats**: Support for GTF and BED files
- **Efficient overlap detection**: Optimized algorithm using sorted intervals
- **Quality filtering**: Filter by mapping quality and read pairing status
- **Feature-specific analysis**: Extract specific features from GTF files (exons, genes, etc.)
- **Progress tracking**: Real-time progress updates for large datasets

## Installation

### Requirements

- Python 3.6 or higher
- pysam library

### Install dependencies

```bash
pip install -r requirements.txt
```

Or install directly:

```bash
pip install pysam
```

## Usage

### Basic Usage

Calculate overlap with a GTF annotation file:

```bash
python calculate_overlap_percentage.py input.bam genes.gtf
```

### With BED Annotations

```bash
python calculate_overlap_percentage.py input.bam genes.bed --format bed
```

### Advanced Options

Apply quality filtering:

```bash
python calculate_overlap_percentage.py input.bam genes.gtf --min-mapq 30 --proper-pairs
```

Specify feature type from GTF:

```bash
python calculate_overlap_percentage.py input.bam genes.gtf --feature-type gene
```

### Command-Line Options

```
positional arguments:
  bam_file              Input BAM file (must be indexed)
  annotation_file       Annotation file (GTF or BED format)

optional arguments:
  -h, --help            Show this help message and exit
  --format {gtf,bed}    Annotation file format (default: gtf)
  --feature-type FEATURE_TYPE
                        Feature type to extract from GTF (default: exon)
  --min-mapq MIN_MAPQ   Minimum mapping quality (default: 0)
  --proper-pairs        Only count properly paired reads
```

## Output

The tool reports:
- Total number of reads analyzed
- Number of reads overlapping annotations
- Percentage of reads overlapping

Example output:

```
============================================================
RESULTS
============================================================
Total reads analyzed: 1,234,567
Reads overlapping annotations: 987,654
Percentage overlapping: 80.00%
============================================================
```

## Input Requirements

- **BAM file**: Must be coordinate-sorted and indexed (`.bai` file in same directory)
- **Annotation file**: GTF or BED format with gene/feature coordinates

## How It Works

1. Loads annotations from GTF/BED file into memory, organized by chromosome
2. Iterates through aligned reads in the BAM file
3. Filters reads based on specified criteria (unmapped, secondary, supplementary alignments are excluded)
4. Checks each read for overlap with annotations on the same chromosome
5. Reports statistics on total and overlapping reads

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Author

Created by [@mrblanco](https://github.com/mrblanco)
