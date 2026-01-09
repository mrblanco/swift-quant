"""
DASH Depletion Statistical Analysis Module

Provides statistical comparison and analysis tools for comparing
gRNA performance across different conditions or experiments.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from scipy import stats
from dataclasses import dataclass

from dash_analyzer import DepletionMetrics


@dataclass
class ComparisonResult:
    """Results from statistical comparison."""
    metric_name: str
    group1_mean: float
    group2_mean: float
    group1_std: float
    group2_std: float
    t_statistic: float
    p_value: float
    significant: bool
    effect_size: float

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'metric': self.metric_name,
            'group1_mean': round(self.group1_mean, 4),
            'group2_mean': round(self.group2_mean, 4),
            'group1_std': round(self.group1_std, 4),
            'group2_std': round(self.group2_std, 4),
            't_statistic': round(self.t_statistic, 4),
            'p_value': round(self.p_value, 4),
            'significant': self.significant,
            'effect_size_cohens_d': round(self.effect_size, 4)
        }


class DASHStatistics:
    """Statistical analysis tools for DASH experiments."""

    @staticmethod
    def calculate_cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
        """
        Calculate Cohen's d effect size.

        Args:
            group1: First group values
            group2: Second group values

        Returns:
            Cohen's d effect size
        """
        n1, n2 = len(group1), len(group2)
        var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)

        # Pooled standard deviation
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

        if pooled_std == 0:
            return 0.0

        return (np.mean(group1) - np.mean(group2)) / pooled_std

    @staticmethod
    def compare_groups(group1_metrics: List[DepletionMetrics],
                      group2_metrics: List[DepletionMetrics],
                      alpha: float = 0.05) -> List[ComparisonResult]:
        """
        Compare two groups of gRNA metrics using t-tests.

        Args:
            group1_metrics: Metrics from first group
            group2_metrics: Metrics from second group
            alpha: Significance level (default: 0.05)

        Returns:
            List of ComparisonResult objects
        """
        metrics_to_compare = [
            ('depletion_efficiency', 'Depletion Efficiency'),
            ('coverage_uniformity', 'Coverage Uniformity'),
            ('zero_coverage_fraction', 'Zero Coverage Fraction'),
            ('coefficient_of_variation', 'Coefficient of Variation')
        ]

        results = []

        for metric_attr, metric_name in metrics_to_compare:
            group1_values = np.array([getattr(m, metric_attr) for m in group1_metrics])
            group2_values = np.array([getattr(m, metric_attr) for m in group2_metrics])

            # Perform t-test
            t_stat, p_val = stats.ttest_ind(group1_values, group2_values)

            # Calculate effect size
            effect_size = DASHStatistics.calculate_cohens_d(group1_values, group2_values)

            results.append(ComparisonResult(
                metric_name=metric_name,
                group1_mean=float(np.mean(group1_values)),
                group2_mean=float(np.mean(group2_values)),
                group1_std=float(np.std(group1_values, ddof=1)),
                group2_std=float(np.std(group2_values, ddof=1)),
                t_statistic=float(t_stat),
                p_value=float(p_val),
                significant=p_val < alpha,
                effect_size=float(effect_size)
            ))

        return results

    @staticmethod
    def rank_grnas(metrics_list: List[DepletionMetrics],
                  weights: Optional[Dict[str, float]] = None) -> pd.DataFrame:
        """
        Rank gRNAs by overall performance using weighted scoring.

        Args:
            metrics_list: List of DepletionMetrics
            weights: Optional dict of weights for each metric
                    Keys: 'depletion', 'uniformity', 'zero_coverage', 'cv'

        Returns:
            DataFrame with gRNAs ranked by composite score
        """
        if weights is None:
            weights = {
                'depletion': 0.4,
                'uniformity': 0.3,
                'zero_coverage': 0.2,
                'cv': 0.1
            }

        scores = []

        for metrics in metrics_list:
            # Normalize metrics to 0-1 scale (higher is better)
            depletion_score = metrics.depletion_efficiency
            uniformity_score = metrics.coverage_uniformity
            zero_cov_score = metrics.zero_coverage_fraction  # Higher = more depletion = better
            cv_score = 1 / (1 + metrics.coefficient_of_variation)  # Lower CV = better

            # Calculate weighted composite score
            composite_score = (
                weights['depletion'] * depletion_score +
                weights['uniformity'] * uniformity_score +
                weights['zero_coverage'] * zero_cov_score +
                weights['cv'] * cv_score
            )

            scores.append({
                'gRNA_name': metrics.grna_name,
                'composite_score': composite_score,
                'depletion_efficiency': metrics.depletion_efficiency,
                'coverage_uniformity': metrics.coverage_uniformity,
                'zero_coverage_fraction': metrics.zero_coverage_fraction,
                'coefficient_of_variation': metrics.coefficient_of_variation
            })

        df = pd.DataFrame(scores)
        df = df.sort_values('composite_score', ascending=False).reset_index(drop=True)
        df['rank'] = range(1, len(df) + 1)

        return df[['rank', 'gRNA_name', 'composite_score', 'depletion_efficiency',
                  'coverage_uniformity', 'zero_coverage_fraction', 'coefficient_of_variation']]

    @staticmethod
    def identify_outliers(metrics_list: List[DepletionMetrics],
                         n_std: float = 2.0) -> Dict[str, List[str]]:
        """
        Identify gRNAs with outlier performance (both good and bad).

        Args:
            metrics_list: List of DepletionMetrics
            n_std: Number of standard deviations for outlier threshold

        Returns:
            Dict with 'high_performers' and 'low_performers' lists
        """
        depletion_values = np.array([m.depletion_efficiency for m in metrics_list])
        uniformity_values = np.array([m.coverage_uniformity for m in metrics_list])

        depletion_mean = np.mean(depletion_values)
        depletion_std = np.std(depletion_values)
        uniformity_mean = np.mean(uniformity_values)
        uniformity_std = np.std(uniformity_values)

        high_performers = []
        low_performers = []

        for metrics in metrics_list:
            depletion_z = (metrics.depletion_efficiency - depletion_mean) / depletion_std if depletion_std > 0 else 0
            uniformity_z = (metrics.coverage_uniformity - uniformity_mean) / uniformity_std if uniformity_std > 0 else 0

            # High performer: significantly better in both metrics
            if depletion_z > n_std and uniformity_z > n_std:
                high_performers.append(metrics.grna_name)

            # Low performer: significantly worse in either metric
            if depletion_z < -n_std or uniformity_z < -n_std:
                low_performers.append(metrics.grna_name)

        return {
            'high_performers': high_performers,
            'low_performers': low_performers
        }

    @staticmethod
    def correlation_analysis(metrics_list: List[DepletionMetrics]) -> pd.DataFrame:
        """
        Analyze correlations between different metrics.

        Args:
            metrics_list: List of DepletionMetrics

        Returns:
            Correlation matrix as DataFrame
        """
        data = {
            'Depletion Efficiency': [m.depletion_efficiency for m in metrics_list],
            'Coverage Uniformity': [m.coverage_uniformity for m in metrics_list],
            'Zero Coverage Fraction': [m.zero_coverage_fraction for m in metrics_list],
            'Coefficient of Variation': [m.coefficient_of_variation for m in metrics_list],
            'Mean Coverage': [m.mean_coverage for m in metrics_list]
        }

        df = pd.DataFrame(data)
        correlation_matrix = df.corr()

        return correlation_matrix

    @staticmethod
    def generate_summary_statistics(metrics_list: List[DepletionMetrics]) -> pd.DataFrame:
        """
        Generate summary statistics for all metrics.

        Args:
            metrics_list: List of DepletionMetrics

        Returns:
            DataFrame with summary statistics
        """
        metrics_dict = {
            'Depletion Efficiency': [m.depletion_efficiency for m in metrics_list],
            'Coverage Uniformity': [m.coverage_uniformity for m in metrics_list],
            'Zero Coverage Fraction': [m.zero_coverage_fraction for m in metrics_list],
            'Coefficient of Variation': [m.coefficient_of_variation for m in metrics_list],
            'Mean Coverage': [m.mean_coverage for m in metrics_list],
            'Median Coverage': [m.median_coverage for m in metrics_list]
        }

        summary = {}
        for metric_name, values in metrics_dict.items():
            values_array = np.array(values)
            summary[metric_name] = {
                'mean': np.mean(values_array),
                'median': np.median(values_array),
                'std': np.std(values_array),
                'min': np.min(values_array),
                'max': np.max(values_array),
                'q25': np.percentile(values_array, 25),
                'q75': np.percentile(values_array, 75)
            }

        return pd.DataFrame(summary).T


def export_comparison_results(results: List[ComparisonResult], output_path: str):
    """
    Export comparison results to CSV.

    Args:
        results: List of ComparisonResult objects
        output_path: Path to save CSV file
    """
    data = [r.to_dict() for r in results]
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
