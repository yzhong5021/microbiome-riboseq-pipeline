"""
Helper functions for taxonomic analysis and visualization.

Functions:
    plot_genus_distributions: Plot genus-level taxonomic distributions
    plot_species_distributions: Plot species-level taxonomic distributions
    plot_individual_distributions: Create individual classifier distribution plots
    create_stacked_bar_plot: Generate stacked bar plots with "Other" categories
    plot_other_category: Visualize breakdown of "Other" category
    extract_genus_from_species: Extract genus information from species names
    calculate_distance_matrix: Compute Aitchison or Bray-Curtis distances
    plot_sankey: Generate Sankey diagrams for classifier overlap
    barplot_hits: Plot hit counts per classifier
"""

import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go


def plot_genus_distributions(df_list, ind_to_df_mapping, cutoff=0.02, out_path=r'figures\barplots\all_genus_02.png'):
    """Plot genus-level taxonomic distributions"""
    
    n = len(df_list)
    dicts = []
    all_categories = set()

    for df in df_list:
        d = {}
        other = 0.0
        for sp, row in df.iterrows():
            prop = row['proportion']
            if prop >= cutoff:
                d[sp] = prop
                all_categories.add(sp)
            else:
                other += prop

        d['Other'] = other
        all_categories.add('Other')
        dicts.append(d)

    categories_for_legend = sorted([cat for cat in all_categories if cat != 'Other']) + ['Other']
    if 'Other' not in all_categories:
        categories_for_legend.remove('Other')
    color_map = plt.colormaps.get_cmap('tab20')
    colors = {cat: color_map(i) for i, cat in enumerate(categories_for_legend)}
    if 'Other' in colors:
        colors['Other'] = (0, 0, 0, 1)

    x = np.arange(n)
    plt.figure(figsize=(12, 8))

    for i, d in enumerate(dicts):
        bottom = 0
        other_val = d.pop('Other', 0)
        sorted_items = sorted(d.items(), key=lambda item: item[1], reverse=True)
        for category_name, proportion in sorted_items:
            plt.bar(i, proportion, bottom=bottom, color=colors.get(category_name, 'grey'), width=0.8)
            bottom += proportion
        if other_val > 0:
            plt.bar(i, other_val, bottom=bottom, color=colors['Other'], width=0.8)

    handles = [plt.Rectangle((0,0),1,1, color=colors[cat]) for cat in categories_for_legend]

    plt.title(f'Genus Distributions (other cutoff = {cutoff*100}%)')
    plt.xticks(x, [ind_to_df_mapping[i] for i in range(n)])

    plt.legend(handles, categories_for_legend, title='Genus', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()


def plot_species_distributions(df_list, ind_to_df_mapping, cutoff=0.03, out_path=r'figures\barplots\all_species_03.png'):
    """Plot species-level taxonomic distributions"""

    n = len(df_list)
    dicts = []
    all_categories = set()

    for df in df_list:
        d = {}
        other = 0.0
        for sp, row in df.iterrows():
            prop = row['proportion']
            sp = sp.strip()
            if prop >= cutoff:
                d[sp] = prop
                all_categories.add(sp)
            else:
                other += prop
        d['Other'] = other
        all_categories.add('Other')
        dicts.append(d)

    categories_for_legend = sorted([cat for cat in all_categories if cat != 'Other']) + ['Other']
    if 'Other' not in all_categories:
        categories_for_legend.remove('Other')
    color_map = plt.colormaps.get_cmap('tab20')
    colors = {cat: color_map(i) for i, cat in enumerate(categories_for_legend)}
    if 'Other' in colors:
        colors['Other'] = (0, 0, 0, 1)

    x = np.arange(n)
    plt.figure(figsize=(12, 8))

    for i, d in enumerate(dicts):
        bottom = 0
        other_val = d.pop('Other', 0)
        sorted_items = sorted(d.items(), key=lambda item: item[1], reverse=True)
        for category_name, proportion in sorted_items:
            plt.bar(i, proportion, bottom=bottom, color=colors.get(category_name, 'grey'), width=0.8)
            bottom += proportion
        if other_val > 0:
            plt.bar(i, other_val, bottom=bottom, color=colors['Other'], width=0.8)

    handles = [plt.Rectangle((0,0),1,1, color=colors[cat]) for cat in categories_for_legend]

    plt.title(f'Species Distributions (other cutoff = {int(cutoff*100)}%)')
    plt.xticks(x, [ind_to_df_mapping[i] for i in range(n)])

    plt.legend(handles, categories_for_legend, title='Species', bbox_to_anchor=(1.04, 1))
    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()


def barplot_hits(counts, classifier_names, out_path=r'figures\classifier_comp\total_hits.png'):
    """Plot hit counts per classifier"""

    plt.figure(figsize=(5, 5))
    plt.bar(classifier_names, counts, color="blue")
    plt.xlabel('')
    for i, v in enumerate(counts):
        plt.text(i, v + max(counts) * 0.01, f'{v}', ha='center', va='bottom', fontsize=8)

    plt.ylabel('')
    plt.title('Total species hits')
    plt.xticks(fontsize=8)

    max_v = max(counts)
    y = int(np.floor(np.log10(max_v)))
    base = 10 ** y

    step = 5 * base
    if max_v / step < 5:
        step = base

    num_ticks = int(np.ceil(max_v / step))

    yticks = np.arange(0, num_ticks * step, step)

    plt.yticks(yticks)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

    plt.tight_layout()
    plt.savefig(out_path)
    plt.show()


def plot_sankey(id_sets, out_path=r'figures\classifier_comp\sankey_counts.html'):
    """Plot Sankey diagram of classifier overlap"""

    # map ID's to classifiers
    id_to_classifier = {}
    for classifier, ids in id_sets.items():
        for _id in ids:
            id_to_classifier.setdefault(_id, set()).add(classifier)

    # count unique classifiers
    group_counts = {}
    for classifier in id_to_classifier.values():
        key = frozenset(classifier)
        group_counts[key] = group_counts.get(key, 0) + 1

    # labels for category combinations
    combo_labels = []
    for classifier, count in group_counts.items():
        if len(classifier) > 1:
            label = " & ".join(sorted(classifier))
        else:
            label = f"{next(iter(classifier))} only"
        combo_labels.append(f"{label} ({count})")

    # individual node labels
    categories = sorted(id_sets)
    cat_labels = [f'{c} ({len(id_sets[c])})' for c in categories]

    # combine node labels
    node_labels = combo_labels + cat_labels
    idx = {lab: i for i, lab in enumerate(node_labels)}

    # create sankey links
    sources, targets, values = [], [], []
    for (classifiers, n), g_label in zip(group_counts.items(), combo_labels):
        for classifier in classifiers:
            sources.append(idx[g_label])
            targets.append(idx[[l for l in cat_labels if l.startswith(classifier)][0]])
            values.append(n)

    # plotting
    fig = go.Figure(go.Sankey(
        node=dict(label=node_labels, pad=10, thickness=30, 
                line=dict(width=0, color='black')),
        link=dict(source=sources, target=targets, value=values,
                    hoverinfo='all'),
    ))

    fig.update_layout(
        height=800, width=1200,
    )

    fig.update_layout(title_text='classified reads', font=dict(size=10, color='black'))
    fig.write_html(out_path)
    fig.show()

    return id_to_classifier


def plot_individual_distributions(df_list, ind_to_df_mapping, ind_to_count_mapping, cutoff=0.01, out_dir=r'figures\barplots\individual'):
    """Plot individual classifier distribution plots"""

    n = len(df_list)
    dicts = []
    all_categories = set()

    for df in df_list:
        d = {}
        other = 0.0
        for sp, row in df.iterrows():
            prop = row['proportion']
            if prop >= cutoff:
                d[sp  + ' (' + str(round(row['proportion'], 3)) + ')'] = prop
                all_categories.add(sp)
            else:
                other += prop
        d['Other'] = other
        all_categories.add('Other')
        dicts.append(d)

    for i, d in enumerate(dicts):
        plt.figure(figsize=(8, 8))
        sample_name = ind_to_df_mapping[i]
        bottom = 0
        other_val = d.pop('Other', 0)
        sorted_items = sorted(d.items(), key=lambda item: item[1], reverse=True)
        current_labels = [item[0] for item in sorted_items]
        if other_val > 0:
            current_labels.append('Other')
        color_map = plt.colormaps.get_cmap('tab20')
        colors = {cat: color_map(i) for i, cat in enumerate(current_labels)}
        if 'Other' in colors:
            colors['Other'] = (0, 0, 0, 1)
        for category_name, proportion in sorted_items:
            plt.bar(0, proportion, bottom=bottom, color=colors.get(category_name, 'grey'), width=0.5)
            bottom += proportion
        if other_val > 0:
            plt.bar(0, other_val, bottom=bottom, color=colors['Other'], width=0.5)
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[cat]) for cat in current_labels]
        plt.title(f'Species Distribution, {sample_name}\n(other cutoff = {int(cutoff*100)}%, n = {ind_to_count_mapping[i]})')
        plt.xticks([])
        plt.ylabel('Proportion')
        plt.ylim(0, 1)
        plt.legend(handles, current_labels, title='Species', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(rf'{out_dir}\{sample_name}_species_{str(cutoff).replace("0.", "")}.png')
    plt.show()


def create_stacked_bar_plot(df, title, cutoff=0.005, out_dir=r'figures\barplots\final'):
    """Create stacked bar plot"""
    
    df = df.set_index('name') if 'name' in df.columns else df
    
    total_count = df['count'].sum()
    
    df_valid = df[df['proportion'] >= cutoff].copy()
    df_other = df[df['proportion'] < cutoff].copy()

    df_valid = df_valid.sort_values(by='proportion', ascending=False)

    if not df_other.empty:
        other_count = df_other['count'].sum()
        other_proportion = other_count / total_count
        df_other_row = pd.DataFrame({'count': other_count, 'proportion': other_proportion}, index=['Other'])
    else:
        df_other_row = pd.DataFrame(columns=['count', 'proportion'])

    final_df = pd.concat([df_valid, df_other_row])

    color_map = plt.colormaps.get_cmap('tab20')
    colors = [color_map(i) for i in range(len(final_df))]
    final_index = list(final_df.index)

    for i, species in enumerate(final_index):
        if species == "Other":
            colors[i] = (0, 0, 0, 1)

    x = [0]
    bottom = 0

    plt.figure(figsize=(10, 8))
    legend_handles = []
    legend_labels = []

    for i, (species, row) in enumerate(final_df.iterrows()):
        height = row['proportion']
        bar = plt.bar(x, height, bottom=bottom, color=colors[i], width=1.0)
        legend_handles.append(bar[0])
        legend_labels.append(f"{species} ({int(row['count'])}, {round(row['proportion'], 3)})")
        bottom += height

    plt.xticks([])
    plt.ylabel("Proportion")
    plt.title(f"{title} (other cutoff = {cutoff*100}%, n = {total_count})")
    plt.legend(legend_handles, legend_labels, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    final_path = os.path.join(out_dir, title.replace(' ', '_') + ".png")
    plt.savefig(final_path)
    plt.show()


def extract_genus_from_species(df_list):
    """Extract genus information from species names"""
    
    df_list_genera = []
    for df in df_list:
        df_copy = df.copy()
        def extract_genus(name):
            parts = name.split()
            if parts[0] == 'KR' or parts[0] == 'TSA:':
                return name.replace(':', '')
            return parts[0]
        df_copy['genus'] = df_copy.index.map(extract_genus)
        df_grouped = df_copy.groupby('genus').sum()
        df_list_genera.append(df_grouped)
    return df_list_genera


def calculate_distance_matrix(df_list_genera, ind_to_df_mapping, distance_type='aitchison', top_n=200):
    """Compute Aitchison or Bray-Curtis distances"""
    
    df_indices_to_compare = range(len(df_list_genera))
    df_labels = [ind_to_df_mapping[i] for i in df_indices_to_compare]
    distance_matrix = pd.DataFrame(0.0, index=df_labels, columns=df_labels)

    if distance_type == 'aitchison':
        def clr(proportions):
            logs = np.log(proportions + 1e-9)
            return logs - np.mean(logs)
        
        t_top_genera = [df.sort_values(by='proportion', ascending=False).iloc[0:top_n] for df in df_list_genera]

        for i in df_indices_to_compare:
            for j in range(i + 1, len(df_list_genera)):
                name1 = ind_to_df_mapping[i]
                name2 = ind_to_df_mapping[j]

                s1 = t_top_genera[i]['proportion']
                s2 = t_top_genera[j]['proportion']

                all_species = s1.index.union(s2.index)

                s1_full = s1.reindex(all_species, fill_value=0.0)
                s2_full = s2.reindex(all_species, fill_value=0.0)

                clr1 = clr(s1_full)
                clr2 = clr(s2_full)

                distance = np.linalg.norm(clr1 - clr2)

                distance_matrix.loc[name1, name2] = distance
                distance_matrix.loc[name2, name1] = distance
    
    elif distance_type == 'bray_curtis':
        for i in df_indices_to_compare:
            for j in range(i + 1, len(df_list_genera)):
                name1 = ind_to_df_mapping[i]
                name2 = ind_to_df_mapping[j]

                s1 = df_list_genera[i]['proportion']
                s2 = df_list_genera[j]['proportion']
                
                bc_dissimilarity = 1 - np.sum(np.minimum(s1, s2))
                
                distance_matrix.loc[name1, name2] = bc_dissimilarity
                distance_matrix.loc[name2, name1] = bc_dissimilarity
    
    return distance_matrix


def plot_other_category(original_df, title_prefix, cutoff=0.005, other_cutoff=0.01, out_dir=r'figures\barplots\final'):
    """Plot 'Other' category distribution"""
    
    df_other = original_df[original_df['proportion'] < cutoff].copy()

    other_total_count = df_other['count'].sum()
    df_other['proportion'] = df_other['count'] / other_total_count
    
    df_other_valid = df_other[df_other['proportion'] >= other_cutoff].copy()
    df_other_other = df_other[df_other['proportion'] < other_cutoff].copy()

    plot_title = f"'Other' distribution, {title_prefix} (n={other_total_count})"

    if not df_other_other.empty:
        other_other_count = df_other_other['count'].sum()
        other_other_proportion = df_other_other['proportion'].sum()
        plot_title = f"'Other' distribution, {title_prefix} (cutoff = {cutoff*100}/{other_cutoff*100}%, n={other_total_count})"
        
        other_other_row = pd.DataFrame({'count': other_other_count, 'proportion': other_other_proportion}, index=['Other'])
        other_other_row.index.name = 'name'
        final_df = pd.concat([df_other_valid, other_other_row])
    else:
        final_df = df_other_valid
    
    final_df = final_df.sort_values(by='proportion', ascending=False)
    if 'name' in final_df.columns:
         final_df = final_df.set_index('name')

    plt.figure(figsize=(10, 8))

    heights = final_df['proportion']
    bottoms = heights.cumsum().shift(1).fillna(0)
    
    color_map = plt.colormaps.get_cmap('tab20')
    colors = [color_map(i) for i in range(len(final_df))]
    
    bars = plt.bar([0] * len(final_df), height=heights, bottom=bottoms, color=colors, width=1.0)

    legend_labels = [f"{name} ({int(count)}, {prop:.3f})" 
                     for name, count, prop in zip(final_df.index, final_df['count'], final_df['proportion'])]

    if 'Other' in final_df.index:
        other_index = final_df.index.get_loc('Other')
        bars[other_index].set_color('black')

    plt.xticks([])
    plt.ylabel("Proportion of 'Other'")
    plt.title(plot_title)
    plt.legend(bars, legend_labels, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    final_path = os.path.join(out_dir, title_prefix.replace(' ', '_') + "_other.png")
    plt.savefig(final_path)
    plt.show()