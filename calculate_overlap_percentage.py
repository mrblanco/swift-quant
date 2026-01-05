#!/usr/bin/env python3
"""
Calculate the percentage of reads in a BAM file that overlap gene annotations.

Requirements:
    pip install pysam pybedtools

Usage:
    python calculate_overlap_percentage.py <bam_file> <annotation_file> [options]
"""

import argparse
import pysam
from collections import defaultdict
import sys


def load_annotations_from_gtf(gtf_file, feature_type='exon'):
    """
    Load gene annotations from a GTF file.

    Args:
        gtf_file: Path to GTF/GFF file
        feature_type: Feature type to extract (default: 'exon')

    Returns:
        Dictionary mapping chromosome to list of (start, end, gene_name) tuples
    """
    annotations = defaultdict(list)

    print(f"Loading annotations from {gtf_file}...", file=sys.stderr)

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
            end = int(fields[4])
            attributes = fields[8]

            if feature == feature_type:
                # Extract gene name from attributes
                gene_name = None
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                        break
                    elif attr.startswith('gene_id') and gene_name is None:
                        gene_name = attr.split('"')[1]

                annotations[chrom].append((start, end, gene_name))

    # Sort annotations by start position for efficient searching
    for chrom in annotations:
        annotations[chrom].sort()

    print(f"Loaded {sum(len(v) for v in annotations.values())} annotations", file=sys.stderr)
    return annotations


def load_annotations_from_bed(bed_file):
    """
    Load gene annotations from a BED file.

    Args:
        bed_file: Path to BED file

    Returns:
        Dictionary mapping chromosome to list of (start, end, name) tuples
    """
    annotations = defaultdict(list)

    print(f"Loading annotations from {bed_file}...", file=sys.stderr)

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else 'unknown'

            annotations[chrom].append((start, end, name))

    # Sort annotations by start position
    for chrom in annotations:
        annotations[chrom].sort()

    print(f"Loaded {sum(len(v) for v in annotations.values())} annotations", file=sys.stderr)
    return annotations


def check_overlap(read_start, read_end, annotations):
    """
    Check if a read overlaps with any annotation.

    Args:
        read_start: Read start position
        read_end: Read end position
        annotations: List of (start, end, name) tuples sorted by start

    Returns:
        True if overlap exists, False otherwise
    """
    for ann_start, ann_end, _ in annotations:
        # Since annotations are sorted, we can stop early
        if ann_start > read_end:
            break

        # Check for overlap
        if read_start < ann_end and read_end > ann_start:
            return True

    return False


def calculate_overlap_percentage(bam_file, annotations, min_mapq=0, require_proper_pair=False):
    """
    Calculate percentage of reads overlapping annotations.

    Args:
        bam_file: Path to BAM file
        annotations: Dictionary of annotations by chromosome
        min_mapq: Minimum mapping quality (default: 0)
        require_proper_pair: Only count properly paired reads (default: False)

    Returns:
        Tuple of (total_reads, overlapping_reads, percentage)
    """
    total_reads = 0
    overlapping_reads = 0

    print(f"Processing BAM file: {bam_file}...", file=sys.stderr)

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch():
            # Apply filters
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            if read.mapping_quality < min_mapq:
                continue

            if require_proper_pair and not read.is_proper_pair:
                continue

            total_reads += 1

            # Get read coordinates
            chrom = read.reference_name
            read_start = read.reference_start
            read_end = read.reference_end

            # Check for overlap
            if chrom in annotations:
                if check_overlap(read_start, read_end, annotations[chrom]):
                    overlapping_reads += 1

            # Progress indicator
            if total_reads % 100000 == 0:
                print(f"Processed {total_reads:,} reads...", file=sys.stderr)

    percentage = (overlapping_reads / total_reads * 100) if total_reads > 0 else 0

    return total_reads, overlapping_reads, percentage


def main():
    parser = argparse.ArgumentParser(
        description='Calculate percentage of BAM reads overlapping gene annotations'
    )
    parser.add_argument('bam_file', help='Input BAM file (must be indexed)')
    parser.add_argument('annotation_file', help='Annotation file (GTF or BED format)')
    parser.add_argument('--format', choices=['gtf', 'bed'], default='gtf',
                        help='Annotation file format (default: gtf)')
    parser.add_argument('--feature-type', default='exon',
                        help='Feature type to extract from GTF (default: exon)')
    parser.add_argument('--min-mapq', type=int, default=0,
                        help='Minimum mapping quality (default: 0)')
    parser.add_argument('--proper-pairs', action='store_true',
                        help='Only count properly paired reads')

    args = parser.parse_args()

    # Load annotations
    if args.format == 'gtf':
        annotations = load_annotations_from_gtf(args.annotation_file, args.feature_type)
    else:
        annotations = load_annotations_from_bed(args.annotation_file)

    # Calculate overlap
    total, overlapping, percentage = calculate_overlap_percentage(
        args.bam_file,
        annotations,
        min_mapq=args.min_mapq,
        require_proper_pair=args.proper_pairs
    )

    # Print results
    print("\n" + "="*60)
    print("RESULTS")
    print("="*60)
    print(f"Total reads analyzed: {total:,}")
    print(f"Reads overlapping annotations: {overlapping:,}")
    print(f"Percentage overlapping: {percentage:.2f}%")
    print("="*60)


if __name__ == '__main__':
    main()
