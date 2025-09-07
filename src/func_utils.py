"""
Helper functions for processing and visualizing functional analysis results.

Functions:
    _gsea_mat: process raw data directory and create gsea enrichment score matrix for functional terms, along with species names
    _prop_mat: process raw data directory and create raw proportion matrix for functional terms, along with species names
    _load_species: Loads and aggregates all individual species count files from a directory into a single Pandas df.
    _stacked_bp: Creates a horizontal stacked bar plot on an existing plot.
    plot_heatmap: Calculates and plots pairwise cosine distances between species based on functional term enrichments.
    pcoa_plot: Performs PCoA analysis on a matrix of functional term enrichments.
    plot_rawct_enrichment: Plots a simple horizontal bar chart of top N enriched pathways among all aggregated species.
    plot_rawct_enrichment_contr: Identifies the top N enriched pathways and plots a stacked bar chart showing the abundance contribution of each species to those pathways.
    analyze_and_plot_gsea: Calculates an aggregate GSEA score for all terms and plots top_n terms for KEGG pathways ('kp') or KOs ('ko').
    plot_top_pathway_abundance: Plots a stacked bar chart of the top N most abundant pathways by summing counts from individual species files.
    translate_ko_ids: Translates KO IDs to KO names.
    translate_pathway_ids: Translates pathway IDs to pathway names.
    translate_ko: Translates KO IDs to KO names.
    translate_pathway: Translates pathway IDs to pathway names.
    translate_merged_ko: Translates merged KO IDs to KO names.
    translate_merged_pathway: Translates merged pathway IDs to pathway names.
    merge_kocounts: Merges KO counts from individual species files into a single file.
    merge_pathwaycounts: Merges pathway counts from individual species files into a single file.
"""

import os
import random
import string
import subprocess
import textwrap

from adjustText import adjust_text
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from skbio.stats.ordination import pcoa
from sklearn.metrics import pairwise_distances

from src import utils

def _gsea_mat(dir_path, term_col):
    """process raw data directory and create gsea enrichment score matrix for functional terms, along with species names"""
    cats = []
    raw_enr = []

    for file in os.listdir(dir_path):
        rfile = os.path.join(dir_path, file)
        name = utils.extract_species_name(rfile)
        cats.append(name)
        bac_df = pd.read_csv(rfile, sep='\t')
        path = bac_df[term_col]
        enr = bac_df['enrichment']

        # gsea results may contain -inf values; replace these with 0 and create a random unique string
        mask = bac_df['enrichment'] == -np.inf
        if mask.any():
            rand_str = ''.join(random.choices(string.ascii_letters, k=5))
            bac_df.loc[mask, 'enrichment'] = 0
            bac_df.loc[mask, term_col] = rand_str

        raw_enr.append(dict(zip(path, enr)))
    
    return cats, raw_enr


def _prop_mat(dir_path, term_col):
    """process raw data directory and create raw proportion matrix for functional terms, along with species names"""
    cats = []
    raw_prop = []

    for file in os.listdir(dir_path):
        rfile = os.path.join(dir_path, file)
        name = utils.extract_species_name(rfile)
        cats.append(name)
        df = pd.read_csv(rfile, sep='\t')   
        total_count = df['count'].sum()
        proportions = dict(zip(df[term_col], df['count']/total_count))
        raw_prop.append(proportions)
        
    return cats, raw_prop


def _load_species(data_path):
    """Loads and aggregates all individual species count files from a directory into a single Pandas df."""

    all_species_counts = []

    for filename in os.listdir(data_path):
        try:
            file_path = os.path.join(data_path, filename)
            name = utils.extract_species_name(filename)
            counts_df = pd.read_csv(file_path, sep='\t', on_bad_lines='skip')
            
            if 'Pathway_Name' in counts_df.columns:
                term_col = 'Pathway_Name'
            elif 'KO_Name' in counts_df.columns:
                term_col = 'KO_Name'
            else:
                print(f"Missing either 'Pathway_Name' or 'KO_Name' in '{filename}'.")
                continue

            if 'count' in counts_df.columns:
                counts_df[term_col] = counts_df[term_col].str.split(' - ').str[0]
                counts_df['species'] = name
                grouped_df = counts_df.groupby([term_col, 'species'], as_index=False)['count'].sum()
                all_species_counts.append(grouped_df)

        except Exception as e:
            print(f"Error with '{filename}': {e}.")
            continue

        
    return pd.concat(all_species_counts, ignore_index=True)


def _stacked_bp(ax, plot_data, title, xlabel, fontsize):
    """Creates a horizontal stacked bar plot on an existing plot."""

    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel('')

    pathways = [textwrap.fill(name, width=60, max_lines=2, placeholder='...') for name in plot_data.index]
    ax.set_yticklabels(pathways, fontsize=fontsize, va='center')
    
    ax.tick_params(axis='x', which='major', labelsize=8, rotation=30)
    ax.tick_params(axis='y', which='major', labelsize=fontsize, length=0)

    # remove borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.legend(
        bbox_to_anchor=(0.5, -0.1),
        ncols=4,
        loc='upper center'
    )


def plot_heatmap(dir_path, title, term_col, dtype, colormap):
    """calculates and plots pairwise cosine distances between species based on functional term enrichments"""
    
    if 'proportions' in dtype:
        cats, raw_mat = _prop_mat(dir_path, term_col)
    elif 'gsea' in dtype:
        cats, raw_mat = _gsea_mat(dir_path, term_col)

    raw_mat_df = pd.DataFrame(raw_mat).T
    raw_mat_df.columns = cats
    raw_mat_df.fillna(0, inplace=True)

    # calculate pairwise cosine distances
    enr_dist = pairwise_distances(raw_mat_df.T, metric='cosine')
    enr_dist_df = pd.DataFrame(enr_dist, index=cats, columns=cats)
    plt.figure(figsize=(10, 8))
    sns.heatmap(enr_dist_df, cmap=colormap, annot=True, fmt='.2f')
    plt.xticks(rotation=30, ha='right')
    plt.title(title)
    plt.savefig(f'figures/functional_analysis/dist_heatmap_{title.lower().replace(' ', '_')}.png')
    plt.show()

    return raw_mat_df, enr_dist_df


def pcoa_plot(mat, title):
    """performs PCoA analysis on a matrix of functional term enrichments"""

    coords = pcoa(mat, dimensions=2)

    coords = coords.samples.values

    plt.figure(figsize=(8, 6))
    plt.scatter(coords[:, 0], coords[:, 1], s=80)

    texts = []
    for i, label in enumerate(mat.index):
        texts.append(plt.text(coords[i, 0], coords[i, 1], label, fontsize=10))

    adjust_text(texts, arrowprops=dict(arrowstyle='fancy', color='black', lw=0.5))

    plt.xlabel('PCoA1')
    plt.ylabel('PCoA2')
    plt.title(f'PCoA of species, cosine distance of {" ".join(title.split(" ")[0:2])} raw proportions')
    plt.savefig(f'figures/functional_analysis/pcoa_{title.lower().replace(' ', '_')}.png')
    plt.tight_layout()
    plt.show()


def plot_rawct_enrichment(df, top_n, title, color, fontsize=8):
    """plots a simple horizontal bar chart of top N enriched pathways among all aggregated species"""

    df = df.sort_values(by='enrichment', ascending=False)
    
    # get top n rows and reverse
    top_data = df.head(top_n)

    fig, ax = plt.subplots(figsize=(8, 8))

    agg_data_reversed = top_data.iloc[::-1]

    agg_data_reversed['pathway'] = agg_data_reversed['pathway'].astype(str)
    pathways = [textwrap.fill(name, width=50) for name in agg_data_reversed['pathway']]

    colors = [color for x in agg_data_reversed['enrichment']]

    ax.barh(y=pathways, width=agg_data_reversed['enrichment'], color=colors, height=0.6)
    ax.set_yticks(np.arange(len(pathways)))
    ax.set_yticklabels(pathways)

    ax.axvline(0, linewidth=0.5)

    ax.tick_params(axis='x', which='major', labelsize=8, rotation=90)
    ax.tick_params(axis='y', which='major', labelsize=fontsize)

    ax.set_title(f'T{top_n} aggregate enriched {" ".join(title.split("_"))}')
    ax.set_xlabel('Abundance')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    plt.savefig(rf'figures\functional_analysis\{title.lower()}_t{top_n}_rawct.png')
    plt.show()


def plot_rawct_enrichment_contr(df, top_n, title, t, fontsize=8):
    """
    Identifies the top N enriched pathways and plots a stacked bar chart showing the
    abundance contribution of each species to those pathways.
    """

    # remove redundant sub-categories of KEGG pathways
    df['pathway_clean'] = df['pathway'].str.split(' - ').str[0]

    # aggregate terms by summing their raw counts; sort by enrichment; get top n functional terms
    aggregated_df = df.groupby('pathway_clean')['enrichment'].sum().reset_index()
    df_sorted = aggregated_df.sort_values(by='enrichment', ascending=False)
    top_pathways_df = df_sorted.head(top_n).copy()
    top_pathways = top_pathways_df['pathway_clean'].tolist()
    
    raw_spcts_path = rf'raw_data\go_analysis\raw_counts\raw_spcts_{t}'
    
    combined_df = _load_species(raw_spcts_path)
    
    col_name = 'Pathway_Name' if 'Pathway_Name' in combined_df.columns else 'KO_Name'
    filtered_df = combined_df[combined_df[col_name].isin(top_pathways)]

    # pivot df to have species as columns and functional terms as rows
    pivot_df = filtered_df.pivot_table(index=col_name, columns='species', values='count', aggfunc='sum').fillna(0)
    species_list = pivot_df.sum().sort_values(ascending=False).index
    pivot_df = pivot_df[species_list]
    pivot_df = pivot_df.reindex(top_pathways).dropna(how='all').fillna(0)
    plot_data = pivot_df.iloc[::-1]

    fig, ax = plt.subplots(figsize=(12, 10))
    plot_data.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', width=0.9)
    
    plot_title = f'Top {top_n} Aggregate Enriched Pathways, \n{" ".join(title.split("_"))}'
    _stacked_bp(ax, plot_data, plot_title, 'Raw count', fontsize)

    fig.tight_layout()
    output_filename = rf'figures\functional_analysis\{title.lower()}_t{top_n}_rawct_stacked.png'
    plt.savefig(output_filename)
    print(f"Saved plot to {output_filename}")
    plt.show()


def analyze_and_plot_gsea(analysis_type, top_n, fontsize=8):
    '''
    Calculates an aggregate GSEA score for all terms and plots top_n terms for KEGG pathways ('kp') or KOs ('ko').
    The aggregate score is a weighted combination of every term's rank for:
        - median GSEA score
        - prevalence (existence of significant positive score) 
    across species.
    '''
    if analysis_type == 'kp':
        data_path = os.path.join('raw_data', 'go_analysis', 'ssgsea_res', 'pathway_gsea')
        file_suffix = '_pathway_gsea.tsv'
        term_col_name = 'pathway'
        title_entity = 'KEGG Pathways'
    elif analysis_type == 'ko':
        data_path = os.path.join('raw_data', 'go_analysis', 'ssgsea_res', 'ko_gsea')
        file_suffix = '_ko_gsea.tsv'
        term_col_name = 'ko_term'
        title_entity = 'KO Terms'

    all_data = []

    # extract GSEA results
    for filename in os.listdir(data_path):
        if filename.endswith(file_suffix):
            file_path = os.path.join(data_path, filename)
            species_name = utils.extract_species_name(filename)

            with open(file_path, 'r', encoding='utf-8') as f:
                next(f)
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.rsplit(None, 1)
                    if len(parts) == 2:
                        term = parts[0].strip()
                        score = float(parts[1])
                        all_data.append([term, score, species_name])
                    else:
                        continue
        
    
    # transform df, get median scores/prevalence scores
    full_df = pd.DataFrame(all_data, columns=[term_col_name, 'score', 'species'])
    pivoted_df = full_df.pivot_table(index=term_col_name, columns='species', values='score')
    imputed_df = pivoted_df.apply(lambda row: row.fillna(row.min()), axis=1)
    med_scores = imputed_df.quantile(0.5, axis=1)
    prevalence_scores = (pivoted_df > 0).sum(axis=1)

    print('\nTop 10 median scores:')
    sorted_med_scores = med_scores.sort_values(ascending=False)
    print(sorted_med_scores.head(10))
    print('\nTop 10 prevalence scores:')
    sorted_prevalence_scores = prevalence_scores.sort_values(ascending=False)
    print(sorted_prevalence_scores.head(10))

    scores_df = pd.DataFrame({
        'med_score': med_scores,
        'prevalence_score': prevalence_scores
    })

    # weighted combination of median/prevalence ranks to get final score
    scores_df['med_rank'] = scores_df['med_score'].rank()
    scores_df['prevalence_rank'] = scores_df['prevalence_score'].rank()
    scores_df['final_score'] = (0.7 * scores_df['prevalence_rank']) + (0.3 * scores_df['med_rank'])
    
    # sort top terms, scale data
    top_terms_df = scores_df.sort_values(by='final_score', ascending=False).head(top_n)
    top_terms_list = top_terms_df.index
    plot_data = imputed_df.loc[top_terms_list].copy()
    plot_data_shifted = plot_data.apply(lambda row: row - row.min(), axis=1)
    row_sums = plot_data_shifted.sum(axis=1)
    row_sums[row_sums == 0] = 1.0
    normalized_df = plot_data_shifted.divide(row_sums, axis=0)
    species_order = normalized_df.sum().sort_values(ascending=False).index
    plot_data_final = normalized_df[species_order]
    plot_data_final = plot_data_final.iloc[::-1]

    fig, ax = plt.subplots(figsize=(12, 10))
    plot_data_final.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', width=0.9)

    plot_title = f'T{top_n} Enriched {title_entity} by GSEA Combined Rank Score'
    xlabel = 'Normalized enrichment contribution by species'
    _stacked_bp(ax, plot_data_final, plot_title, xlabel, fontsize)

    fig.tight_layout()
    output_dir = os.path.join('figures', 'functional_analysis')
    os.makedirs(output_dir, exist_ok=True)
    output_filename = os.path.join(output_dir, f'{analysis_type}_gsea_top{top_n}_stacked.png')
    
    plt.savefig(output_filename, bbox_inches='tight', dpi=300)
    plt.show()


def plot_top_pathway_abundance(top_n, title, t, fontsize=8):
    """
    Plots a stacked bar chart of the top N most abundant pathways by summing
    counts from individual species files.
    """
    raw_spcts_path = rf'raw_data\go_analysis\raw_counts\raw_spcts_{t}'
    
    combined_df = _load_species(raw_spcts_path)

    if combined_df.empty:
        print("No valid species data was loaded. Cannot generate plot.")
        return

    # aggregate counts by term
    col_name = 'Pathway_Name' if 'Pathway_Name' in combined_df.columns else 'KO_Name'
    pathway_totals = combined_df.groupby(col_name)['count'].sum()

    # get top n pathways, pivot to make species into columns
    top_pathway_names = pathway_totals.sort_values(ascending=False).head(top_n).index.tolist()
    filtered_df = combined_df[combined_df[col_name].isin(top_pathway_names)]
    pivot_df = filtered_df.pivot_table(index=col_name, columns='species', values='count', aggfunc='sum').fillna(0)
    species_order = pivot_df.sum().sort_values(ascending=False).index
    pivot_df = pivot_df[species_order]
    pivot_df = pivot_df.reindex(top_pathway_names)
    plot_data = pivot_df.iloc[::-1]

    if plot_data.empty:
        print("Pivoted data is empty. Cannot generate plot.")
        return

    fig, ax = plt.subplots(figsize=(12, 10))
    plot_data.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', width=0.9)
    
    plot_title = f'T{top_n} Aggregate Enriched Pathways, {" ".join(title.split("_"))}'
    _stacked_bp(ax, plot_data, plot_title, 'Total Raw Count', fontsize)

    fig.tight_layout()
    output_filename = rf'figures\functional_analysis\{title.lower()}_t{top_n}_abundance_stacked.png'
    print(f"Saving plot to: {output_filename}")
    plt.savefig(output_filename, bbox_inches='tight', dpi=300)
    plt.show()

def translate_ko_ids(ko_counts_file):
    """
    Translates the ko_id column in a species-specific counts file to descriptive names,
    then overwrites the file.
    """
    try:
        kn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/ko_to_name.tsv",
            sep='\t',
            header=None,
            names=['ko', 'actual_name']
        )
        kn_mapping = pd.Series(kn_df.actual_name.values, index=kn_df.ko).to_dict()

        counts_df = pd.read_csv(ko_counts_file, sep='\t')

        if 'ko_id' not in counts_df.columns:
            print(f"Warning: 'ko_id' column not found in {ko_counts_file}. Skipping translation.")
            return

        cleaned_ko_ids = counts_df['ko_id'].str.replace('ko:', '', regex=False)
        counts_df['KO_Name'] = cleaned_ko_ids.map(kn_mapping)
        counts_df = counts_df[['ko_id', 'KO_Name', 'count']]

        print("Translated KO IDs! Sample:")
        print(counts_df.head())

        counts_df.to_csv(ko_counts_file, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: KO mapping file not found. Skipping translation of {ko_counts_file}")
    except Exception as e:
        print(f"An error occurred during KO translation for {ko_counts_file}: {e}")

def translate_pathway_ids(pathway_counts_file):
    """
    Translates the pathway_id column in a species-specific counts file to descriptive names,
    then overwrites the file.
    """
    try:
        pn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/pathway_to_name.tsv",
            sep='\t',
            header=None,
            names=['pathway', 'actual_name']
        )
        pn_mapping = pd.Series(pn_df.actual_name.values, index=pn_df.pathway).to_dict()

        counts_df = pd.read_csv(pathway_counts_file, sep='\t')

        if 'pathway_id' not in counts_df.columns:
            print(f"Warning: 'pathway_id' column not found in {pathway_counts_file}. Skipping translation.")
            return

        counts_df['pathway'] = counts_df['pathway_id'].map(pn_mapping)
        counts_df = counts_df[['pathway_id', 'pathway', 'count']]

        print("Translated pathway IDs! Sample:")
        print(counts_df.head())

        counts_df.to_csv(pathway_counts_file, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: Pathway mapping file not found. Skipping translation of {pathway_counts_file}")
    except Exception as e:
        print(f"An error occurred during pathway translation for {pathway_counts_file}: {e}")

def translate_ko(ko_path):
    """
    ko_path: path to the ko set.
    performed before ko gsea analysis.
    """
    try:
        kn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/ko_to_name.tsv",
            sep='\t',
            header=None,
            names=['ko', 'actual_name']
        )
        kn_mapping = pd.Series(kn_df.actual_name.values, index=kn_df.ko).to_dict()

        raw_df = pd.read_csv(ko_path, sep='\t')
        raw_df['ko'] = raw_df['ko'].str.replace('ko:', '', regex=False)
        raw_df['ko_name'] = raw_df['ko'].map(kn_mapping)
        raw_df = raw_df[['ko', 'ko_name', 'enrichment']]

        print("translated KEGG KO names! sample list:")
        print(raw_df.head(n = 5))

        raw_df.to_csv(ko_path, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: KO mapping file not found. Skipping translation of {ko_path}")
    except Exception as e:
        print(f"An error occurred during KO translation: {e}")

def translate_pathway(pw_path):
    """
    pw_path: path to the pathway set.
    performed before the pathway gsea analysis.
    """
    try:
        pn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/pathway_to_name.tsv",
            sep='\t',
            header=None,
            names=['pathway', 'actual_name']
        )
        pn_mapping = pd.Series(pn_df.actual_name.values, index=pn_df.pathway).to_dict()

        raw_df = pd.read_csv(pw_path, sep='\t')
        raw_df['pathway'] = raw_df['pathway'].map(pn_mapping)
        
        print("translated KEGG pathway names! sample list:")
        print(raw_df.head(n = 5))

        raw_df.to_csv(pw_path, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: Pathway mapping file not found. Skipping translation of {pw_path}")
    except Exception as e:
        print(f"An error occurred during pathway translation: {e}")

def translate_merged_ko(merged_ko_path):
    """
    Translates the KO_ID column in the merged ko_counts.tsv file to descriptive names.
    Overwrites the input file with the translated results.
    """
    try:
        kn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/ko_to_name.tsv",
            sep='\t',
            header=None,
            names=['ko', 'actual_name']
        )
        kn_mapping = pd.Series(kn_df.actual_name.values, index=kn_df.ko).to_dict()

        counts_df = pd.read_csv(merged_ko_path, sep='\t')
        
        if 'KO_ID' not in counts_df.columns:
            print(f"Warning: 'KO_ID' column not found in {merged_ko_path}. Skipping translation.")
            return merged_ko_path
        
        cleaned_ko_ids = counts_df['KO_ID'].str.replace('ko:', '', regex=False)
        counts_df['KO'] = cleaned_ko_ids.map(kn_mapping)
        counts_df = counts_df[['KO', 'count']]

        print("Translated merged KO counts! Sample:")
        print(counts_df.head(n=5))

        counts_df.to_csv(merged_ko_path, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: KO mapping file not found. Skipping translation of {merged_ko_path}")
    except Exception as e:
        print(f"An error occurred during merged KO translation: {e}")
    return merged_ko_path

def translate_merged_pathway(merged_pathway_path):
    """
    Translates the PATHWAY_ID column in the merged pathway_counts.tsv file to descriptive names.
    Overwrites the input file with the translated results.
    """
    try:
        pn_df = pd.read_csv(
            r"/home/yochen/microbiome/go_analysis/ko/pathway_to_name.tsv",
            sep='\t',
            header=None,
            names=['pathway', 'actual_name']
        )
        pn_mapping = pd.Series(pn_df.actual_name.values, index=pn_df.pathway).to_dict()

        counts_df = pd.read_csv(merged_pathway_path, sep='\t')
        
        if 'PATHWAY_ID' not in counts_df.columns:
            print(f"Warning: 'PATHWAY_ID' column not found in {merged_pathway_path}. Skipping translation.")
            return merged_pathway_path

        counts_df['pathway'] = counts_df['PATHWAY_ID'].map(pn_mapping)
        counts_df = counts_df[['pathway', 'count']]

        print("Translated merged pathway counts! Sample:")
        print(counts_df.head())

        counts_df.to_csv(merged_pathway_path, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Warning: Pathway mapping file not found. Skipping translation of {merged_pathway_path}")
    except Exception as e:
        print(f"An error occurred during merged pathway translation: {e}")
    return merged_pathway_path

def merge_kocounts(counts_dir, ko_sets_dir, out_path, ko_delim='[;,|]', type='KO_ID'):
    """
    1. Concatenate every tsv in counts_dir (genes, coutns) into combined/counts_agg.tsv
    2. Concatenate every tsv in ko_sets_dir (QUERY_NAME, KEGG_KO) into combined/ko_mapping.tsv
       Expected cols: QUERY_NAME, KEGG_KO
    3. Sum raw counts for each KO term (pooled across genes) into ko_counts.tsv

    ko_delim: describes, in regex, how KO IDs are separated
    """
    combined_dir = utils.create_dir(out_path, "combined")
    
    counts_files = [f for f in os.listdir(counts_dir) if f.endswith('.tsv')]
    ko_files = [f for f in os.listdir(ko_sets_dir) if f.endswith('.tsv')]

    counts_dfs = []
    for file in counts_files:
        df = pd.read_csv(os.path.join(counts_dir, file), sep='\t')
        counts_dfs.append(df)
    counts_df = pd.concat(counts_dfs, ignore_index=True)

    ko_dfs = []
    for file in ko_files:
        df = pd.read_csv(os.path.join(ko_sets_dir, file), sep='\t')
        ko_dfs.append(df)
    ko_df = pd.concat(ko_dfs, ignore_index=True)

    counts_agg_path = os.path.join(combined_dir, "counts_agg.tsv")
    ko_mapping_path = os.path.join(combined_dir, "ko_mapping.tsv")
    
    counts_df.to_csv(counts_agg_path, sep='\t', index=False)
    ko_df.to_csv(ko_mapping_path, sep='\t', index=False)

    merged = counts_df.merge(
        ko_df,
        left_on="gene",
        right_on="QUERY_NAME",
        how="inner"
    )

    merged[type] = merged['name'].astype(str).str.split(ko_delim, regex=True)
    exploded = merged.explode(type)
    exploded[type] = exploded[type].str.strip()

    exploded = exploded[
        exploded[type].str.lower().notna() & 
        (exploded[type].str.lower() != 'unknown') & 
        (exploded[type] != '')
    ]

    ko_counts = (
        exploded.groupby(type, as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )

    ko_agg_path = os.path.join(combined_dir, "ko_counts.tsv")
    ko_counts.to_csv(ko_agg_path, sep='\t', index=False)

    print(f"Merged {len(ko_counts)} KOs from {len(counts_files)} species")
    return ko_agg_path

def merge_pathwaycounts(counts_dir, path_sets_dir, out_path, path_delim='[;,|]'):
    """
    1. Concatenate every *.tsv in *counts_dir*
    2. Concatenate every *.tsv in *path_sets_dir*
    3. Sum raw counts for each Pathway ID (pooled across genes)
       → combined/pathway_counts.tsv (cols: PATHWAY_ID, count)

    *path_delim* is a regex describing how pathway IDs are separated.
    """
    combined_dir = utils.create_dir(out_path, "combined")
    
    counts_files = [f for f in os.listdir(counts_dir) if f.endswith('.tsv')]
    path_files = [f for f in os.listdir(path_sets_dir) if f.endswith('.tsv')]

    counts_dfs = []
    for file in counts_files:
        df = pd.read_csv(os.path.join(counts_dir, file), sep='\t')
        counts_dfs.append(df)
    counts_df = pd.concat(counts_dfs, ignore_index=True)

    path_dfs = []
    for file in path_files:
        df = pd.read_csv(os.path.join(path_sets_dir, file), sep='\t')
        path_dfs.append(df)
    path_df = pd.concat(path_dfs, ignore_index=True)

    merged = counts_df.merge(
        path_df,
        left_on="gene",
        right_on="QUERY_NAME",
        how="inner"
    )

    merged['PATHWAY_ID'] = merged['name'].astype(str).str.split(path_delim, regex=True)
    exploded = merged.explode('PATHWAY_ID')
    exploded['PATHWAY_ID'] = exploded['PATHWAY_ID'].str.strip()

    exploded = exploded[
        exploded['PATHWAY_ID'].str.lower().notna() & 
        (exploded['PATHWAY_ID'].str.lower() != 'unknown') & 
        (exploded['PATHWAY_ID'] != '')
    ]

    exploded = exploded[~exploded['PATHWAY_ID'].str.contains('ko', case=False, na=False)]

    pathway_counts = (
        exploded.groupby('PATHWAY_ID', as_index=False)["count"]
        .sum()
        .sort_values("count", ascending=False)
    )

    output_file = os.path.join(combined_dir, "pathway_counts.tsv")
    pathway_counts.to_csv(output_file, sep='\t', index=False)

    print(f"Merged {len(pathway_counts)} pathways from {len(counts_files)} species")
    return output_file

def gsea(tpm_tsv, set_tsv, type, out_path, species):
    """runs the following command:
    Rscript ~/microbiome/go_analysis/scripts/ssgsea.r -t {tpm_tsv} -g {ko_set_tsv} -o {out_path/ssgsea_res} -s {species}"""
    utils.create_dir(out_path, "")
    
    r_script_path = os.path.expanduser("~/microbiome/go_analysis/scripts/ko_ssgsea.r")

    if not os.path.exists(r_script_path):
        print(f"Error: R script not found at {r_script_path}")
        return

    cmd = [
        "Rscript", r_script_path,
        "-t", tpm_tsv,
        "-g", set_tsv,
        "-o", out_path,
        "-s", species
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running GSEA for {species}: {e}")
        return

    return out_path + f"/{species}_{type}_gsea.tsv"