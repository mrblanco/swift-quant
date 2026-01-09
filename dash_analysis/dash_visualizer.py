"""
DASH Depletion Visualization Module

Provides visualization tools for DASH depletion analysis including
coverage plots, heatmaps, and performance metrics.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Optional, Dict
import pandas as pd

from dash_analyzer import DASHAnalyzer, GRNATarget, DepletionMetrics


class DASHVisualizer:
    """Visualization tools for DASH analysis."""

    def __init__(self, style: str = 'seaborn-v0_8-darkgrid'):
        """Initialize visualizer with plot style."""
        try:
            plt.style.use(style)
        except:
            plt.style.use('default')

        sns.set_palette("husl")

    def plot_coverage_comparison(self, analyzer: DASHAnalyzer, target: GRNATarget,
                                output_path: Optional[str] = None):
        """
        Plot coverage comparison between treated and control samples.

        Args:
            analyzer: DASHAnalyzer instance
            target: GRNATarget to visualize
            output_path: Optional path to save figure
        """
        # Get coverage data
        treated_cov = analyzer.get_coverage(target.chrom, target.start, target.end)
        positions = np.arange(target.start, target.end)

        fig, ax = plt.subplots(figsize=(12, 4))

        # Plot treated
        ax.fill_between(positions, treated_cov, alpha=0.6, label='Treated', color='#e74c3c')

        # Plot control if available
        if analyzer.control_bam:
            control_cov = analyzer.get_coverage(target.chrom, target.start, target.end,
                                               analyzer.control_bam)
            ax.fill_between(positions, control_cov, alpha=0.4, label='Control', color='#3498db')

        ax.set_xlabel(f'Position on {target.chrom}')
        ax.set_ylabel('Coverage')
        ax.set_title(f'Coverage Profile: {target.name}\n{target.chrom}:{target.start}-{target.end}')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_multiple_targets_coverage(self, analyzer: DASHAnalyzer,
                                       targets: List[GRNATarget],
                                       output_path: Optional[str] = None):
        """
        Plot coverage for multiple targets in a grid.

        Args:
            analyzer: DASHAnalyzer instance
            targets: List of GRNATarget objects
            output_path: Optional path to save figure
        """
        n_targets = len(targets)
        n_cols = min(3, n_targets)
        n_rows = (n_targets + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 3*n_rows))
        if n_targets == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for idx, target in enumerate(targets):
            ax = axes[idx]
            treated_cov = analyzer.get_coverage(target.chrom, target.start, target.end)
            positions = np.arange(target.start, target.end)

            ax.fill_between(positions, treated_cov, alpha=0.6, color='#e74c3c')

            if analyzer.control_bam:
                control_cov = analyzer.get_coverage(target.chrom, target.start, target.end,
                                                   analyzer.control_bam)
                ax.fill_between(positions, control_cov, alpha=0.4, color='#3498db')

            ax.set_xlabel(f'Position')
            ax.set_ylabel('Coverage')
            ax.set_title(f'{target.name}', fontsize=10)
            ax.grid(True, alpha=0.3)

        # Hide unused subplots
        for idx in range(n_targets, len(axes)):
            axes[idx].axis('off')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_depletion_heatmap(self, metrics_list: List[DepletionMetrics],
                              output_path: Optional[str] = None):
        """
        Create heatmap of depletion metrics across gRNAs.

        Args:
            metrics_list: List of DepletionMetrics
            output_path: Optional path to save figure
        """
        # Prepare data
        data = []
        grna_names = []

        for metrics in metrics_list:
            grna_names.append(metrics.grna_name)
            data.append([
                metrics.depletion_efficiency,
                metrics.coverage_uniformity,
                1 - metrics.zero_coverage_fraction,  # Invert so higher is better
                1 / (1 + metrics.coefficient_of_variation)  # Normalized CV
            ])

        df = pd.DataFrame(data,
                         columns=['Depletion\nEfficiency', 'Coverage\nUniformity',
                                 'Coverage\nCompleteness', 'Normalized\nVariation'],
                         index=grna_names)

        # Create heatmap
        fig, ax = plt.subplots(figsize=(8, max(6, len(grna_names) * 0.4)))
        sns.heatmap(df, annot=True, fmt='.3f', cmap='RdYlGn', center=0.5,
                   vmin=0, vmax=1, ax=ax, cbar_kws={'label': 'Score'})

        ax.set_title('gRNA Performance Metrics Heatmap\n(Higher values indicate better performance)')
        ax.set_ylabel('gRNA')
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def plot_metrics_comparison(self, metrics_list: List[DepletionMetrics],
                               output_path: Optional[str] = None):
        """
        Create bar plots comparing key metrics across gRNAs.

        Args:
            metrics_list: List of DepletionMetrics
            output_path: Optional path to save figure
        """
        grna_names = [m.grna_name for m in metrics_list]
        depletion_eff = [m.depletion_efficiency for m in metrics_list]
        uniformity = [m.coverage_uniformity for m in metrics_list]
        zero_frac = [m.zero_coverage_fraction for m in metrics_list]

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Depletion efficiency
        axes[0].bar(range(len(grna_names)), depletion_eff, color='#3498db', alpha=0.7)
        axes[0].set_xticks(range(len(grna_names)))
        axes[0].set_xticklabels(grna_names, rotation=45, ha='right')
        axes[0].set_ylabel('Depletion Efficiency')
        axes[0].set_title('Depletion Efficiency by gRNA')
        axes[0].set_ylim([0, 1])
        axes[0].axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='Target: 0.5')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3, axis='y')

        # Coverage uniformity
        axes[1].bar(range(len(grna_names)), uniformity, color='#2ecc71', alpha=0.7)
        axes[1].set_xticks(range(len(grna_names)))
        axes[1].set_xticklabels(grna_names, rotation=45, ha='right')
        axes[1].set_ylabel('Coverage Uniformity')
        axes[1].set_title('Coverage Uniformity by gRNA')
        axes[1].set_ylim([0, 1])
        axes[1].axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='Target: 0.5')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3, axis='y')

        # Zero coverage fraction
        axes[2].bar(range(len(grna_names)), zero_frac, color='#e74c3c', alpha=0.7)
        axes[2].set_xticks(range(len(grna_names)))
        axes[2].set_xticklabels(grna_names, rotation=45, ha='right')
        axes[2].set_ylabel('Zero Coverage Fraction')
        axes[2].set_title('Zero Coverage Fraction by gRNA')
        axes[2].set_ylim([0, 1])
        axes[2].grid(True, alpha=0.3, axis='y')

        plt.tight_layout()

        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def create_summary_report(self, metrics_list: List[DepletionMetrics],
                            output_path: str):
        """
        Create a comprehensive PDF report with all visualizations.

        Args:
            metrics_list: List of DepletionMetrics
            output_path: Path to save PDF report
        """
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(output_path) as pdf:
            # Page 1: Metrics comparison
            self.plot_metrics_comparison(metrics_list)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            # Page 2: Heatmap
            self.plot_depletion_heatmap(metrics_list)
            pdf.savefig(bbox_inches='tight')
            plt.close()

            # Add metadata
            d = pdf.infodict()
            d['Title'] = 'DASH Depletion Analysis Report'
            d['Author'] = 'DASH Analyzer'
            d['Subject'] = 'gRNA Performance Metrics'
            d['Keywords'] = 'DASH, depletion, gRNA, analysis'


def export_metrics_to_csv(metrics_list: List[DepletionMetrics], output_path: str):
    """
    Export metrics to CSV file.

    Args:
        metrics_list: List of DepletionMetrics
        output_path: Path to save CSV file
    """
    data = [m.to_dict() for m in metrics_list]
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)


def export_metrics_to_excel(metrics_list: List[DepletionMetrics], output_path: str):
    """
    Export metrics to Excel file with formatting.

    Args:
        metrics_list: List of DepletionMetrics
        output_path: Path to save Excel file
    """
    data = [m.to_dict() for m in metrics_list]
    df = pd.DataFrame(data)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name='Metrics', index=False)

        # Get workbook and worksheet
        workbook = writer.book
        worksheet = writer.sheets['Metrics']

        # Auto-adjust column widths
        for idx, col in enumerate(df.columns):
            max_length = max(df[col].astype(str).apply(len).max(), len(col)) + 2
            worksheet.column_dimensions[chr(65 + idx)].width = max_length
