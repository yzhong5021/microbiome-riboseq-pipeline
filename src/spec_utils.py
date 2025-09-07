"""
Helper functions for specificity analysis and visualization.

Functions:
    plot_alignment_heatmap: Generate alignment heatmaps with top value highlighting
    plot_corrected_heatmap: Plot corrected alignment heatmaps
    normalize: Normalize alignment matrices by pangenome sizes
    plot_alignment_comparison: Compare internal vs external alignments
    plot_corrected_alignment_subplots: Create subplot visualizations
    calculate_log2_fc: Calculate log2 fold changes with error propagation
    plot_log2_fc: Visualize log2 fold change results
"""

from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_alignment_heatmap(df, title, out_path, figsize=(11, 11), highlight_top=True):
    """plot alignment heatmap."""
    plt.figure(figsize=figsize)
    df_r = df.copy().T
    sns.heatmap(df_r, cmap='viridis', annot=True, annot_kws={"size": 8}, fmt=".3f", 
                square=True, cbar=False, cbar_kws={'shrink': 0.7})
    plt.yticks(fontsize=11)
    plt.xticks(fontsize=11, rotation=30, ha='right')

    if highlight_top:
        for y, row in enumerate(df_r.values):
            top2 = row.argsort()[-2:][::-1]
            plt.gca().add_patch(
                plt.Rectangle((top2[0], y), 1, 1, fill=False, edgecolor='violet', lw=3)
            )
            plt.gca().add_patch(
                plt.Rectangle((top2[1], y), 1, 1, fill=False, edgecolor='pink', lw=3)
            )

    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path)


def plot_corrected_heatmap(df, title, out_path, figsize=(12, 12)):
    """plot corrected alignment heatmap."""
    plt.figure(figsize=figsize)
    sns.heatmap(df, cmap='viridis', annot=True, annot_kws={"size": 8}, fmt=".3f", square=True)
    plt.xticks(fontsize=8, rotation=45, ha='right')
    plt.yticks(fontsize=8)
    plt.title(title)
    plt.tight_layout()

    for y, row in enumerate(df.values):
        top2 = row.argsort()[-2:][::-1]
        plt.gca().add_patch(
            plt.Rectangle((top2[0], y), 1, 1, fill=False, edgecolor='white', lw=3)
        )
        plt.gca().add_patch(
            plt.Rectangle((top2[1], y), 1, 1, fill=False, edgecolor='gray', lw=3)
        )

    plt.savefig(out_path)


def normalize(df, pg_size_mapping):
    """normalize alignment similarity matrix by pangenome sizes."""
    df_norm = df.copy()
    
    for ind in df_norm.index:
        real_name = ind.split('_')[0]
        size = pg_size_mapping[real_name]
        df_norm.loc[ind] = df_norm.loc[ind] / size
    
    return df_norm


def plot_alignment_comparison(ast_df, se_int, sem_ext, out_path, title='Internal vs. external alignment'):
    """plot internal vs external alignment comparison with error bars."""
    n_int = ast_df['int_unal_ct'] + ast_df['int_al1_ct'] + ast_df['int_al2_ct']
    p_int = ast_df['int_al']
    yerr_int = 1.96 * se_int
    yerr_ext = 1.96 * sem_ext

    fig, ax = plt.subplots(figsize=(12, 6))
    width = 0.35
    x = np.arange(len(ast_df.index))
    ax.bar(x - width/2, ast_df['int_al'], width=width, color='blue', label='int_al', yerr=yerr_int, capsize=4)
    ax.bar(x + width/2, ast_df['ext_al_mu'], width=width, color='red', label='ext_al_mu', yerr=yerr_ext, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(ast_df.index, rotation=45, ha='right')
    ax.set_ylabel('Proportion')
    ax.set_title(title)
    ax.legend()
    plt.savefig(out_path)
    plt.tight_layout()
    plt.show()


def plot_corrected_alignment_subplots(ast_df, yerr_int_cor, yerr_ext_cor, out_path):
    """plot corrected internal vs external alignment as individual subplots."""
    sorted_species = sorted(ast_df.index)
    n_items = len(sorted_species)
    n_cols = n_items
    n_rows = (n_items + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4))
    axes = axes.flatten()

    for i, species_name in enumerate(sorted_species):
        ax = axes[i]

        int_val = ast_df.loc[species_name, 'int_al_cor']
        ext_val = ast_df.loc[species_name, 'ext_al_cor_mu']
        int_err = yerr_int_cor.loc[species_name]
        ext_err = yerr_ext_cor.loc[species_name]

        ax.bar(
            ["int", "ext"],
            [int_val, ext_val],
            yerr=[int_err, ext_err],
            color=['blue', 'red'],
            capsize=5,
            align='center',
            zorder=2
        )
        
        ax.set_ylim(bottom=0)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_xlabel(species_name, fontsize=8, rotation=45, ha='right')

    for i in range(n_items, len(axes)):
        axes[i].set_visible(False)

    fig.suptitle('Corrected internal vs. external alignment', fontsize=16)
    legend_handles = [
        Patch(color='blue', label='Internal alignment'),
        Patch(color='red', label='External alignment')
    ]
    fig.legend(
        handles=legend_handles,
        loc='upper right',
        ncol=1,
        fontsize=8
    )

    plt.subplots_adjust(bottom=0.25, left=0.4, wspace=0.2, hspace=0.4)
    plt.savefig(out_path)
    plt.show()


def calculate_log2_fc(ast_df, se_values, sem_values, int_col='int_al', ext_col='ext_al_mu', epsilon=1e-9):
    """calculate log2 fold change with error bars"""
    se_log2fc = (1 / np.log(2)) * np.sqrt(
        (se_values**2 / (ast_df[int_col]**2 + epsilon)) + (sem_values**2 / (ast_df[ext_col]**2 + epsilon))
    )
    
    yerr_log2fc = 1.96 * se_log2fc
    
    log2_fold_change = np.log2((ast_df[int_col] + epsilon) / (ast_df[ext_col] + epsilon))
    
    return log2_fold_change, yerr_log2fc


def plot_log2_fc(ast_df, log2_fold_change, yerr_log2fc, title, out_path):
    """plot log2 fc bar plot."""
    fig, ax = plt.subplots(figsize=(12, 6))

    colors = ['red' if x > 0 else 'blue' for x in log2_fold_change]
    ax.bar(
        ast_df.index,
        log2_fold_change,
        yerr=yerr_log2fc,
        color=colors,
        capsize=4
    )

    ax.axhline(0, color='black', linewidth=0.8)
    ax.set_ylabel('Log2 FC (int/ext)')
    ax.set_title(title)
    ax.tick_params(axis='x', rotation=30)
    plt.setp(ax.get_xticklabels(), ha='right')
    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()