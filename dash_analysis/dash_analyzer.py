"""
DASH Depletion Analysis Core Module

This module provides core functionality for analyzing DASH depletion reactions,
including BAM file processing, coverage calculation, and depletion metrics.
"""

import pysam
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class GRNATarget:
    """Represents a gRNA target region."""
    name: str
    chrom: str
    start: int
    end: int
    sequence: str = ""
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start

    def __str__(self):
        return f"{self.name} ({self.chrom}:{self.start}-{self.end})"


@dataclass
class DepletionMetrics:
    """Metrics for gRNA depletion performance."""
    grna_name: str
    target_region: str
    mean_coverage: float
    median_coverage: float
    depletion_efficiency: float
    coverage_uniformity: float
    zero_coverage_fraction: float
    coefficient_of_variation: float

    def to_dict(self) -> Dict:
        """Convert metrics to dictionary."""
        return {
            'gRNA_name': self.grna_name,
            'target_region': self.target_region,
            'mean_coverage': round(self.mean_coverage, 2),
            'median_coverage': round(self.median_coverage, 2),
            'depletion_efficiency': round(self.depletion_efficiency, 4),
            'coverage_uniformity': round(self.coverage_uniformity, 4),
            'zero_coverage_fraction': round(self.zero_coverage_fraction, 4),
            'coefficient_of_variation': round(self.coefficient_of_variation, 4)
        }


class DASHAnalyzer:
    """Main analyzer class for DASH depletion experiments."""

    def __init__(self, bam_path: str, control_bam_path: Optional[str] = None):
        """
        Initialize DASH analyzer.

        Args:
            bam_path: Path to treated BAM file
            control_bam_path: Optional path to control/untreated BAM file
        """
        self.bam_path = Path(bam_path)
        self.control_bam_path = Path(control_bam_path) if control_bam_path else None
        self.bam = None
        self.control_bam = None

    def __enter__(self):
        """Context manager entry."""
        self.bam = pysam.AlignmentFile(str(self.bam_path), "rb")
        if self.control_bam_path:
            self.control_bam = pysam.AlignmentFile(str(self.control_bam_path), "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self.bam:
            self.bam.close()
        if self.control_bam:
            self.control_bam.close()

    def get_coverage(self, chrom: str, start: int, end: int,
                     bam_file: Optional[pysam.AlignmentFile] = None) -> np.ndarray:
        """
        Get per-base coverage for a genomic region.

        Args:
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position (exclusive)
            bam_file: Optional specific BAM file to use

        Returns:
            Array of coverage values for each position
        """
        bam = bam_file if bam_file else self.bam
        coverage = np.zeros(end - start, dtype=np.int32)

        for pileup_column in bam.pileup(chrom, start, end, truncate=True):
            if start <= pileup_column.pos < end:
                coverage[pileup_column.pos - start] = pileup_column.n

        return coverage

    def calculate_uniformity(self, coverage: np.ndarray) -> float:
        """
        Calculate coverage uniformity using coefficient of variation.

        Lower values indicate more uniform coverage.
        Returns 1.0 - CV to make higher values better.

        Args:
            coverage: Array of coverage values

        Returns:
            Uniformity score (0-1, higher is better)
        """
        if len(coverage) == 0 or np.mean(coverage) == 0:
            return 0.0

        cv = np.std(coverage) / np.mean(coverage)
        # Return inverse normalized CV (higher = more uniform)
        return 1.0 / (1.0 + cv)

    def calculate_depletion_efficiency(self, treated_coverage: np.ndarray,
                                      control_coverage: Optional[np.ndarray] = None) -> float:
        """
        Calculate depletion efficiency.

        If control is provided: 1 - (treated_mean / control_mean)
        If no control: based on fraction of bases with zero coverage

        Args:
            treated_coverage: Coverage in treated sample
            control_coverage: Coverage in control sample (optional)

        Returns:
            Depletion efficiency (0-1, higher is better)
        """
        if control_coverage is not None:
            treated_mean = np.mean(treated_coverage)
            control_mean = np.mean(control_coverage)

            if control_mean == 0:
                return 0.0

            efficiency = 1.0 - (treated_mean / control_mean)
            return max(0.0, min(1.0, efficiency))
        else:
            # Without control, use zero coverage fraction as proxy
            if len(treated_coverage) == 0:
                return 0.0
            return np.sum(treated_coverage == 0) / len(treated_coverage)

    def analyze_grna_target(self, target: GRNATarget) -> DepletionMetrics:
        """
        Analyze a single gRNA target region.

        Args:
            target: GRNATarget object

        Returns:
            DepletionMetrics object
        """
        # Get coverage for treated sample
        treated_cov = self.get_coverage(target.chrom, target.start, target.end)

        # Get coverage for control if available
        control_cov = None
        if self.control_bam:
            control_cov = self.get_coverage(target.chrom, target.start, target.end,
                                           self.control_bam)

        # Calculate metrics
        mean_cov = np.mean(treated_cov)
        median_cov = np.median(treated_cov)
        uniformity = self.calculate_uniformity(treated_cov)
        depletion_eff = self.calculate_depletion_efficiency(treated_cov, control_cov)
        zero_frac = np.sum(treated_cov == 0) / len(treated_cov) if len(treated_cov) > 0 else 0
        cv = np.std(treated_cov) / mean_cov if mean_cov > 0 else 0

        return DepletionMetrics(
            grna_name=target.name,
            target_region=f"{target.chrom}:{target.start}-{target.end}",
            mean_coverage=float(mean_cov),
            median_coverage=float(median_cov),
            depletion_efficiency=float(depletion_eff),
            coverage_uniformity=float(uniformity),
            zero_coverage_fraction=float(zero_frac),
            coefficient_of_variation=float(cv)
        )

    def analyze_multiple_targets(self, targets: List[GRNATarget]) -> List[DepletionMetrics]:
        """
        Analyze multiple gRNA targets.

        Args:
            targets: List of GRNATarget objects

        Returns:
            List of DepletionMetrics objects
        """
        results = []
        for target in targets:
            metrics = self.analyze_grna_target(target)
            results.append(metrics)
        return results


def load_grna_targets_from_bed(bed_path: str) -> List[GRNATarget]:
    """
    Load gRNA targets from a BED file.

    BED format: chrom start end name [score] [strand]

    Args:
        bed_path: Path to BED file

    Returns:
        List of GRNATarget objects
    """
    targets = []
    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            strand = fields[5] if len(fields) > 5 else "+"

            targets.append(GRNATarget(
                name=name,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand
            ))

    return targets
