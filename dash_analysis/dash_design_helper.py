"""
DASH gRNA Design Helper

Provides recommendations for improving gRNA design based on analysis results.
"""

import numpy as np
from typing import List, Dict, Tuple
from dataclasses import dataclass

from dash_analyzer import DepletionMetrics


@dataclass
class DesignRecommendation:
    """Recommendation for gRNA design improvement."""
    grna_name: str
    severity: str  # 'critical', 'warning', 'info'
    category: str
    issue: str
    recommendation: str
    priority: int  # 1-5, 5 being highest priority


class DASHDesignHelper:
    """Analyze metrics and provide gRNA design recommendations."""

    # Thresholds for identifying issues
    THRESHOLDS = {
        'depletion_critical': 0.3,
        'depletion_warning': 0.6,
        'uniformity_critical': 0.3,
        'uniformity_warning': 0.5,
        'cv_warning': 1.0,
        'cv_critical': 2.0,
        'zero_coverage_low': 0.5,
    }

    @staticmethod
    def analyze_grna_performance(metrics: DepletionMetrics) -> List[DesignRecommendation]:
        """
        Analyze a single gRNA and provide recommendations.

        Args:
            metrics: DepletionMetrics for the gRNA

        Returns:
            List of DesignRecommendation objects
        """
        recommendations = []

        # Check depletion efficiency
        if metrics.depletion_efficiency < DASHDesignHelper.THRESHOLDS['depletion_critical']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='critical',
                category='Depletion Efficiency',
                issue=f'Very low depletion efficiency ({metrics.depletion_efficiency:.3f})',
                recommendation='Consider complete redesign. Check: (1) gRNA sequence matches '
                              'reference genome, (2) no off-target binding, (3) target region '
                              'is accessible (not in heterochromatin)',
                priority=5
            ))
        elif metrics.depletion_efficiency < DASHDesignHelper.THRESHOLDS['depletion_warning']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='warning',
                category='Depletion Efficiency',
                issue=f'Moderate depletion efficiency ({metrics.depletion_efficiency:.3f})',
                recommendation='Consider optimization: (1) Try alternative target location, '
                              '(2) Increase gRNA concentration, (3) Extend incubation time',
                priority=3
            ))

        # Check coverage uniformity
        if metrics.coverage_uniformity < DASHDesignHelper.THRESHOLDS['uniformity_critical']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='critical',
                category='Coverage Uniformity',
                issue=f'Very poor coverage uniformity ({metrics.coverage_uniformity:.3f})',
                recommendation='Uneven depletion suggests: (1) Strong secondary structures '
                              'in target region - use RNA structure prediction tools, '
                              '(2) Repetitive sequences causing mapping issues, '
                              '(3) Consider using multiple shorter gRNAs instead',
                priority=4
            ))
        elif metrics.coverage_uniformity < DASHDesignHelper.THRESHOLDS['uniformity_warning']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='warning',
                category='Coverage Uniformity',
                issue=f'Moderate coverage uniformity ({metrics.coverage_uniformity:.3f})',
                recommendation='Check: (1) GC content is 40-60%, (2) Avoid strong hairpins, '
                              '(3) Target region is not in low complexity sequence',
                priority=3
            ))

        # Check coefficient of variation
        if metrics.coefficient_of_variation > DASHDesignHelper.THRESHOLDS['cv_critical']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='critical',
                category='Coverage Variation',
                issue=f'Very high coverage variation (CV={metrics.coefficient_of_variation:.3f})',
                recommendation='Highly inconsistent depletion. Possible causes: '
                              '(1) PCR amplification bias, (2) Secondary structures, '
                              '(3) Mapping artifacts. Consider redesigning target region',
                priority=4
            ))
        elif metrics.coefficient_of_variation > DASHDesignHelper.THRESHOLDS['cv_warning']:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='warning',
                category='Coverage Variation',
                issue=f'High coverage variation (CV={metrics.coefficient_of_variation:.3f})',
                recommendation='Moderate inconsistency in depletion. Consider: '
                              '(1) Checking for repetitive elements, '
                              '(2) Increasing sequencing depth',
                priority=2
            ))

        # Check zero coverage
        if metrics.zero_coverage_fraction < DASHDesignHelper.THRESHOLDS['zero_coverage_low']:
            if metrics.depletion_efficiency >= 0.5:  # Only if depletion is otherwise OK
                recommendations.append(DesignRecommendation(
                    grna_name=metrics.grna_name,
                    severity='info',
                    category='Coverage Completeness',
                    issue=f'Incomplete depletion (only {metrics.zero_coverage_fraction:.1%} '
                          'bases with zero coverage)',
                    recommendation='Consider: (1) Increasing DASH reaction time, '
                                  '(2) Increasing Cas9/gRNA concentration, '
                                  '(3) Using overlapping gRNAs for this region',
                    priority=2
                ))

        # Edge effects detection (this is a simplified heuristic)
        if metrics.coefficient_of_variation > 0.8 and metrics.coverage_uniformity < 0.6:
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='info',
                category='Edge Effects',
                issue='Possible edge effects in coverage',
                recommendation='Extend target region by 50-100bp on each side to ensure '
                              'complete coverage of intended region',
                priority=2
            ))

        # Good performance recognition
        if (metrics.depletion_efficiency > 0.8 and
            metrics.coverage_uniformity > 0.7):
            recommendations.append(DesignRecommendation(
                grna_name=metrics.grna_name,
                severity='info',
                category='Performance',
                issue='Excellent performance',
                recommendation='This gRNA shows excellent depletion efficiency and uniformity. '
                              'No changes needed. Consider using similar design principles '
                              'for other gRNAs.',
                priority=0
            ))

        return recommendations

    @staticmethod
    def generate_report(metrics_list: List[DepletionMetrics]) -> str:
        """
        Generate a comprehensive design recommendations report.

        Args:
            metrics_list: List of DepletionMetrics

        Returns:
            Formatted report string
        """
        all_recommendations = []

        for metrics in metrics_list:
            recs = DASHDesignHelper.analyze_grna_performance(metrics)
            all_recommendations.extend(recs)

        # Sort by priority
        all_recommendations.sort(key=lambda x: x.priority, reverse=True)

        # Generate report
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("DASH gRNA Design Recommendations Report")
        report_lines.append("=" * 80)
        report_lines.append("")

        # Summary statistics
        critical_count = sum(1 for r in all_recommendations if r.severity == 'critical')
        warning_count = sum(1 for r in all_recommendations if r.severity == 'warning')
        info_count = sum(1 for r in all_recommendations if r.severity == 'info')

        report_lines.append(f"Total gRNAs analyzed: {len(metrics_list)}")
        report_lines.append(f"Critical issues: {critical_count}")
        report_lines.append(f"Warnings: {warning_count}")
        report_lines.append(f"Info/Good performers: {info_count}")
        report_lines.append("")

        # Critical issues first
        if critical_count > 0:
            report_lines.append("=" * 80)
            report_lines.append("CRITICAL ISSUES (Immediate attention required)")
            report_lines.append("=" * 80)
            report_lines.append("")

            for rec in [r for r in all_recommendations if r.severity == 'critical']:
                report_lines.append(f"gRNA: {rec.grna_name}")
                report_lines.append(f"Category: {rec.category}")
                report_lines.append(f"Issue: {rec.issue}")
                report_lines.append(f"Recommendation: {rec.recommendation}")
                report_lines.append("-" * 80)
                report_lines.append("")

        # Warnings
        if warning_count > 0:
            report_lines.append("=" * 80)
            report_lines.append("WARNINGS (Should be addressed)")
            report_lines.append("=" * 80)
            report_lines.append("")

            for rec in [r for r in all_recommendations if r.severity == 'warning']:
                report_lines.append(f"gRNA: {rec.grna_name}")
                report_lines.append(f"Category: {rec.category}")
                report_lines.append(f"Issue: {rec.issue}")
                report_lines.append(f"Recommendation: {rec.recommendation}")
                report_lines.append("-" * 80)
                report_lines.append("")

        # Good performers
        good_performers = [r for r in all_recommendations
                          if r.severity == 'info' and r.priority == 0]
        if good_performers:
            report_lines.append("=" * 80)
            report_lines.append("HIGH PERFORMING gRNAs")
            report_lines.append("=" * 80)
            report_lines.append("")

            for rec in good_performers:
                report_lines.append(f"gRNA: {rec.grna_name}")
                report_lines.append(f"Recommendation: {rec.recommendation}")
                report_lines.append("")

        return "\n".join(report_lines)

    @staticmethod
    def export_recommendations_csv(metrics_list: List[DepletionMetrics],
                                   output_path: str):
        """
        Export recommendations to CSV file.

        Args:
            metrics_list: List of DepletionMetrics
            output_path: Path to save CSV
        """
        import csv

        all_recommendations = []
        for metrics in metrics_list:
            recs = DASHDesignHelper.analyze_grna_performance(metrics)
            all_recommendations.extend(recs)

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'grna_name', 'severity', 'category', 'issue',
                'recommendation', 'priority'
            ])
            writer.writeheader()

            for rec in sorted(all_recommendations, key=lambda x: x.priority, reverse=True):
                writer.writerow({
                    'grna_name': rec.grna_name,
                    'severity': rec.severity,
                    'category': rec.category,
                    'issue': rec.issue,
                    'recommendation': rec.recommendation,
                    'priority': rec.priority
                })


def main():
    """Example usage of the design helper."""
    from dash_analyzer import DASHAnalyzer, load_grna_targets_from_bed

    print("DASH Design Helper - Analyzing gRNA performance...\n")

    # Load and analyze
    targets = load_grna_targets_from_bed('example_targets.bed')

    with DASHAnalyzer('treated.bam', 'control.bam') as analyzer:
        metrics_list = analyzer.analyze_multiple_targets(targets)

    # Generate report
    helper = DASHDesignHelper()
    report = helper.generate_report(metrics_list)

    print(report)

    # Export to CSV
    helper.export_recommendations_csv(metrics_list, 'design_recommendations.csv')
    print("\nRecommendations exported to design_recommendations.csv")


if __name__ == '__main__':
    main()
