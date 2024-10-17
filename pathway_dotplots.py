#!/usr/bin/env python

import argparse
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import glob
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import linkage, fcluster

def hierarchical_clustering(data, max_d=50, cluster_columns=None):
    '''
    Perform hierarchical clustering of data
    :param data:
    :param max_d:  distance threshold
    :param cluster_columns:
    :return:
    '''
    dfs = list(data.values())
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='pathway', suffixes=('', '_y'))
    merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]

    if cluster_columns is not None:
        merged_df = merged_df[['pathway'] + cluster_columns]
    features = merged_df.columns.drop('pathway')
    Z = linkage(merged_df[features], 'ward')
    clusters = fcluster(Z, max_d, criterion='distance')
    merged_df['cluster'] = clusters
    sorted_df = merged_df.sort_values(by='cluster',ascending=True)
    return sorted_df, clusters

def dotplot(data, list_of_samples, pvalue_column, pvalue_threshold, plot_type, set_of_path_ways):
    '''
    Create dotplot for pathway analysis
    :param data:
    :param list_of_samples: It should be list of samples
    :param pvalue_column:
    :param pvalue_threshold:
    :param plot_type:
    :param set_of_path_ways:
    :return:
    '''
    fig = plt.figure(figsize=[11, 12])
    min_color = 1.0
    max_color = 0
    min_size = 1.0
    max_size = 0
    enrich_column = 'NES'
    if plot_type == "invsize":
        size_column = enrich_column
        color_column = pvalue_column
    else:
        size_column = pvalue_column
        color_column = enrich_column
    for smp in list_of_samples:
        min_size_current = data[smp][size_column].min()
        max_size_current = data[smp][size_column].max()
        print(f"Min {pvalue_column} for {smp} is {min_size_current}, max {pvalue_column} is {max_size_current}")
        if min_size > min_size_current:
            min_size = min_size_current
        if max_size < max_size_current:
            max_size = max_size_current
        min_color_current = data[smp][color_column].min()
        max_color_current = data[smp][color_column].max()
        if min_color > min_color_current:
            min_color = min_color_current
        if max_color < max_color_current:
            max_color = max_color_current
    if size_column == pvalue_column:
        if max_size > 1.0:
            raise RuntimeError(f"Some p-values in column pval greater than 1.0 ")
        if min_size < 0.0:
            raise RuntimeError(f"Some p-values in column pval less than 0.0 ")
        print(f"min_size:{min_size}, max_size:{max_size}, min_color:{min_color}, max_color:{max_color}")
    i = 0
    # create list of pathways where all values greater threshold:
    skip_pathways = []
    accept_pathways = []
    for pw in set_of_path_ways:
        empty_num = 0
        for smp in list_of_samples:
            below_threshold = data[smp][data[smp]['pathway'] == pw][pvalue_column] < pvalue_threshold
            empty_num += int(below_threshold)
        if empty_num == 0:
            skip_pathways.append(pw)
        else:
            accept_pathways.append(pw)
    set_of_path_ways = accept_pathways
    for smp in list_of_samples:
        j = 0
        df = data[smp]
        for pw in set_of_path_ways:
            if plot_type == "invsize":
                msize = df[df['pathway'] == pw][enrich_column]
            else:
                msize = df[df['pathway'] == pw][pvalue_column]
            if len(msize) == 0:
                continue
            if len(msize) > 1:
                print(f" {smp} have more than one pathway {pw}. Will be skipped")
                continue
            if plot_type == "size":
                msize_inv = -np.log10(msize.iloc[0])
                mcolor = df[df['pathway'] == pw][enrich_column].iloc[0]
                bc, rc = (abs(mcolor) / abs(min_color), 0) if mcolor < 0 else (0, abs(mcolor) / abs(max_color))
                plt.scatter(i, j, s=50 * (msize_inv) / 1.6, color=(rc, 0.0, bc), alpha=1)
            elif plot_type == "invsize":
                msize_inv = msize.iloc[0] / max_size * 25
                mcolor =  -np.log10(df[df['pathway'] == pw][pvalue_column].iloc[0]) / -np.log10(min_color)
                bc, rc = (abs(mcolor) , 0) if mcolor < 0 else (0, abs(mcolor))
                plt.scatter(i, j, s=25 * (msize_inv) , color=(rc, 0.0, bc), alpha=1)
            elif plot_type == "gray":
                msize_inv = 5  # -np.log10(msize.iloc[0])
                mcolor = df[df['pathway'] == pw][enrich_column].iloc[0]
                bc, rc = (abs(mcolor) / abs(min_color), 0) if mcolor < 0 else (0, abs(mcolor) / abs(max_color))
                if msize.iloc[0] <= pvalue_threshold:
                    plt.scatter(i, j, s=25 * (msize_inv), color=(rc, 0.0, bc), alpha=1)
                else:
                    plt.scatter(i, j, s=25 * (msize_inv), color=(0.5, 0.5, 0.5), alpha=1)
            elif plot_type == "white":
                msize_inv = 5  # -np.log10(msize.iloc[0])
                mcolor = df[df['pathway'] == pw][enrich_column].iloc[0]
                if min_color >= 0:
                    bc, rc = (0, abs(mcolor) / abs(max_color))
                elif max_color <= 0:
                    bc, rc = (abs(mcolor) / abs(min_color), 0)
                else:
                    bc, rc = (abs(mcolor) / abs(min_color), 0) if mcolor < 0 else (0, abs(mcolor) / abs(max_color))
                if msize.iloc[0] <= pvalue_threshold:
                    plt.scatter(i, j, s=25 * (msize_inv), color=(rc, 0.0, bc), alpha=1)
            else:
                raise RuntimeError(f"Unknown plot type:{plot_type}")

            j += 1
        i += 1
    ax = plt.gca()
    fig = plt.gcf()
    plt.rc('axes', axisbelow=True)
    plt.xticks(np.linspace(0, i - 1, i))
    ax.set_xticklabels(list_of_samples, fontsize=10, rotation=90, ha='right')
    plt.yticks(np.linspace(0, j - 1, j))
    ax.set_yticklabels(set_of_path_ways, fontsize=8, rotation=0)
    # handles, labels = plt.gca().get_legend_handles_labels()

    plt.title("Gene set enrichment analysis of Hallmark pathways")
    if min_color < 0 and max_color > 0:
        colors = ['#0000FF', '#000088', '#000000', '#880000',
              '#FF0000']
        labels_text = [min_color, min_color / 2, 0, max_color / 2, max_color]
    else:
        labels_text = [min_color, (max_color-min_color)/4,(max_color-min_color)/4*2, (max_color-min_color)/4*3, max_color]
        colors = ['#000000', '#330000', '#660000', '#990000',
              '#FF0000']
    # Create patches for the custom legend
    if plot_type == "size":
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'{labels_text[i]:.2f}',
                                  markersize=10, markerfacecolor=color) for i, color in enumerate(colors)]
        title_legend_elements = [Line2D([0], [0], marker='', color='w', label='Enrichment')] + legend_elements

        labels2_text = [min_size, (min_size + max_size) * 10 ** (np.log10(min_size) / 4 * 3),
                        (min_size + max_size) * 10 ** (np.log10(min_size) / 2),
                        (min_size + max_size) * 10 ** (np.log10(min_size) / 4), max_size]
        title_legend_elements2 = [Line2D([0], [0], marker='', color='w', label='P-values')] + [Line2D([0], [0],
                                                                                                      marker='o',
                                                                                                      color='w',
                                                                                                      label=f'10^{np.log10(labels2_text[i]):.2f}',
                                                                                                      markersize=1.3 * -np.log10(
                                                                                                          labels2_text[
                                                                                                              i]) / 1.3,
                                                                                                      markerfacecolor='#FF0000')
                                                                                               for i, color in
                                                                                               enumerate(colors)]
        all_legend = title_legend_elements + title_legend_elements2
        ax.legend(handles=all_legend, title='', loc='upper right', bbox_to_anchor=(1.45, 0.75))
    elif plot_type == "invsize":
        #reverse list colors
        colors = colors[::-1]
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'{labels_text[i]:.2f}',
                                  markersize=10, markerfacecolor=color) for i, color in enumerate(colors)]
        title_legend_elements = [Line2D([0], [0], marker='', color='w', label='P-value')] + legend_elements
        labels2_text = [min_size,
                        (max_size-min_size)/4,
                        (max_size-min_size)/4*2,
                        (max_size-min_size)/4*3,
                        max_size]
        title_legend_elements2 = [Line2D([0], [0], marker='', color='w', label='Enrichment')] + [Line2D([0], [0],
                                                                                                      marker='o',
                                                                                                      color='w',
                                                                                                      label=f'{labels2_text[i]:.2f}',
                                                                                                      markersize= labels2_text[i] / max_size * 25,
                                                                                                      markerfacecolor=color)
                                                                                               for i, color in
                                                                                               enumerate(colors)]
        all_legend = title_legend_elements + title_legend_elements2
        ax.legend(handles=all_legend, title='', loc='upper right', bbox_to_anchor=(1.45, 0.75))

    elif plot_type == "gray":
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=f'{labels_text[i]:.2f}',
                                  markersize=10, markerfacecolor=color) for i, color in enumerate(colors)]
        title_legend_elements = [Line2D([0], [0], marker='', color='w', label='P-value')] + legend_elements
        colors = ['#0000FF', '#FF0000', '#888888']
        labels_text = [f'p_value =< {pvalue_threshold}', f'p_value =< {pvalue_threshold}', f'p_value >'
                                                                                           f' {pvalue_threshold}']
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=labels_text[i],
                                  markersize=10, markerfacecolor=color) for i, color in enumerate(colors)]
        title_legend_elements2 = [Line2D([0], [0], marker='', color='w', label='P-values')] + legend_elements
        all_legend = title_legend_elements + title_legend_elements2
        ax.legend(handles=all_legend, title='', loc='upper right', bbox_to_anchor=(1.45, 0.75))
    else:
        from matplotlib import colors as mcolors
        from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        title_legend_elements2 = []
        norm = mcolors.Normalize(min_color, max_color)
        N = 256
        blue = np.ones((N, 4))
        blue[:, 0] = np.zeros(N)
        blue[:, 1] = np.zeros(N)
        blue[:, 2] = np.linspace(1, 0, N)
        blue_cmp = ListedColormap(blue)
        red = np.ones((N, 4))
        red[:, 0] = np.linspace(1, 0, N)
        red[:, 1] = np.zeros(N)
        red[:, 2] = np.zeros(N)
        red_cmp = ListedColormap(red)
        if min_color >= 0:
            blue_red_cmp = ListedColormap(red_cmp(np.linspace(1, 0, 256)), name="blue_red_cmp")
        elif max_color <= 0:
            blue_red_cmp = ListedColormap(blue_cmp(np.linspace(0, 1, 256)), name="blue_red_cmp")
        else:
            blue_red_cmp = ListedColormap(np.vstack((blue_cmp(np.linspace(0, 1, 128)),
                                                     red_cmp(np.linspace(1, 0, 128)))), name="blue_red_cmp")
        cmap = blue_red_cmp
        cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Enrichment', shrink=0.4)
        cbar.ax.locator_params(nbins=20)

    plt.tight_layout()
    plt.grid(color='gray', linestyle='dashed')
    return fig
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="pathway-dotplots - Tool for ploting dot-plot for pathway analysis",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-indir", help="Directory with input *tsv files", type=str,
                        required=True)
    parser.add_argument("--input_file_mask", help="Input file mask", type=str,
                        default='*.csv')

    parser.add_argument("-outdir", help="Directory where pdf with plots will be saved", type=str,
                        required=True)
    parser.add_argument("--plot_type", help="Type of plot. Options:"
                                            "'size' - size of dot will be -log10(pvalue), "
                                            "'invsize' - size of dot will be proportional to enrichment, "
                                            "'gray' - dots with pvalue more than threshold (default 0.05) will be ploted in gray color"
                                            "'white - dots with pvalue more than threshold (default 0.05) will not be ploted",
                                            type=str,
                        default='size')
    parser.add_argument("--pathway_sort", help="Type of pathways sort. Supported options 'cluster',"
                                              "'first_enrich','enrich_threshold'", type=str,
                        default='cluster')
    parser.add_argument("--pvalue_threshold", help="pvalue threshold for 'gray' plot option ", type=float,
                        default=0.1)
    parser.add_argument("--show", help="If set, plots will be shown", default=False,
                        action='store_true')
    parser.add_argument("--imgname", help="filename for output plot", type=str,
                        default='dot_plot')

    args = parser.parse_args()
    indir = args.indir
    outdir = args.outdir
    plot_type = args.plot_type
    imgname = args.imgname
    show = args.show
    pvalue_column = 'padj'
    pvalue_threshold = args.pvalue_threshold
    files = glob.glob(f"{indir}/{args.input_file_mask}")
    mandatory_columns = ['pathway','NES',pvalue_column]
    data = {}
    list_of_pathways = []
    list_of_samples = []
    script_dir = os.path.dirname(os.path.realpath(__file__))
    theme_bw = os.path.join(script_dir,"theme_bw.mplstyle")
    plt.style.use(theme_bw)
    for fi in files:
        df = pd.read_csv(fi,delimiter=',')
        for mc in mandatory_columns:
            if mc not in df.columns:
                print(f"Mandatory column {mc} not found in {os.path.basename(fi)}. File will be skipped")
        list_of_samples.append(os.path.basename(fi).split('.csv')[0])
        data[list_of_samples[-1]] = df
        list_of_pathways += df['pathway'].to_list()

    def reverse_sort(s):
        return s[::-1]
    list_of_samples = sorted(list_of_samples,reverse=True)

    if (args.pathway_sort == 'cluster'):
        cldfs,clusters = hierarchical_clustering(data,cluster_columns=['NES'],max_d=10)
        set_of_path_ways = cldfs['pathway'].to_list()
    elif (args.pathway_sort == 'first_enrich'):
        sort_order_data = data[list_of_samples[0]]
        set_of_path_ways = sort_order_data.sort_values(by=['NES'],ascending=True)['pathway'].to_list()
    elif (args.pathway_sort == 'enrich_threshold'):
        sort_order_data = data[list_of_samples[0]]
        sort_order_data['pval_thres'] = sort_order_data[pvalue_column] <= pvalue_threshold
        set_of_path_ways = sort_order_data.sort_values(by=['pval_thres','NES'],ascending=True)['pathway'].to_list()
    else:
        raise RuntimeError(f"Unknown pathway sort type:{args.pathway_sort}")
    fig = dotplot(data, list_of_samples, pvalue_column, pvalue_threshold, plot_type, set_of_path_ways)
    if show:
        plt.show()
    else:
        plt.savefig(os.path.join(outdir,f'{imgname}_{plot_type}.pdf'))