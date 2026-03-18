import igraph as ig
import pandas as pd
import argparse
import math
import random

import igraph as ig
import plotly.graph_objects as go
import matplotlib.cm as cm

import os
from ctypes import *
import ctypes
import numpy as np


def read_graph_from_gml(gml_file_w_path):
    g = ig.Graph.Read_GML(gml_file_w_path)
    return g

def community_voronoi(g, weights):
    """
    Perform community detection using the Voronoi method.
    """
    if weights is None or len(weights) == 0:
        weights = g.es['weight']
            
    # Get absolute path to the .so file
    so_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'my_functions.so')
            # Convert numpy array to C compatible format

    lib = ctypes.CDLL(so_file)

    adjacency_matrix = g.get_adjacency(attribute='weight').data
    size = len(adjacency_matrix)
    
    # Create array of pointers to rows
    ptr_type = ctypes.POINTER(ctypes.c_double)
    row_pointers = (ptr_type * size)()
    
    # Create C arrays for each row
    rows = []
    for i in range(size):
        row = (ctypes.c_double * size)(*adjacency_matrix[i])
        rows.append(row)
        row_pointers[i] = row

    # Set function argument types and return type
    lib.community_voronoi.argtypes = [ctypes.POINTER(ptr_type), ctypes.c_int]
    lib.community_voronoi.restype = ctypes.POINTER(ctypes.c_int)

    # Call C function
    result_ptr = lib.community_voronoi(row_pointers, size)
    if not result_ptr:
        return None

    # Convert result to numpy array
    result = np.array([result_ptr[i] for i in range(size)])
    
    # Free C allocated memory
    libc = ctypes.CDLL(None)
    libc.free(result_ptr)
    result = abs(result)
    new_membership = pd.factorize(result)[0].tolist()
    return ig.VertexClustering(g, membership=new_membership)

def voronoi_clustering(g):
    return g.community_voronoi(weights=g.es['weight'])

def infomap_clustering(g):
    return g.community_infomap(edge_weights=g.es['weight'])

def calculate_betweenness_centrality(g):
    return g.betweenness(directed=True, weights=g.es['weight'])

def calculate_closeness(g):
    return g.closeness(mode='all', weights=g.es['weight'])

def calculate_energy_function(g, clustering, gene_p_values, key='mean_pvalue_all', energy_type='betweenness'):
    energies = {}
    if energy_type == 'betweenness':
        centrality = calculate_betweenness_centrality(g)
    elif energy_type == 'closeness':
        centrality = calculate_closeness(g)
    for idx, cluster in enumerate(clustering.membership):
        gene_in_cluster = g.vs[idx]['name']
        p_val = gene_p_values[gene_p_values['Gene_0'].isin([gene_in_cluster])][key]
        p_val = p_val.min() if not p_val.empty else 0
        energies[gene_in_cluster] = centrality[cluster] * - math.log10(p_val + 1e-300)
    return energies

def read_gene_p_values(file_w_path):
    df = pd.read_csv(file_w_path)
    return df

def get_p_value_stats_per_cluster(g, gene_p_values, clustering_p, key='mean_pvalue_all'):
    cluster_avg_pvalues = []
    energies_betweenness = calculate_energy_function(g, clustering_p, gene_p_values, key, energy_type='betweenness')
    energies_closeness = calculate_energy_function(g, clustering_p, gene_p_values, key, energy_type='closeness')
    print(clustering_p)
    for cluster in clustering_p:
        # if len(cluster) == 0:
        #     continue
        print(">>> ", cluster)
        genes_in_cluster = [g.vs[vid]['name'] for vid in cluster]
        print(genes_in_cluster)
        p_values = gene_p_values[gene_p_values['Gene_0'].isin(genes_in_cluster)][key]

        mean_p_value = p_values.mean() if not p_values.empty else 0
        median_p_value = p_values.median() if not p_values.empty else 0
        min_p_value = p_values.min() if not p_values.empty else 0
        energy_betweenness = sum(energies_betweenness[gene] for gene in genes_in_cluster if gene in energies_betweenness)
        energy_closeness = sum(energies_closeness[gene] for gene in genes_in_cluster if gene in energies_closeness)
        cluster_avg_pvalues.append({
            'cluster_size': len(genes_in_cluster),
            'mean_p_value': -math.log10(mean_p_value + 1e-300),
            'median_p_value': -math.log10(median_p_value + 1e-300),
            'min_p_value': -math.log10(min_p_value + 1e-300),
            'energy_betweenness': energy_betweenness,
            'energy_closeness': energy_closeness,
            'genes': genes_in_cluster
        })
    return cluster_avg_pvalues

def calculate_stats_main(graph, clustering, gene_p_values_df, key='mean_pvalue_all', clustering_name="Voronoi clustering"):
    clustering_stats = get_p_value_stats_per_cluster(graph, gene_p_values_df, clustering, key)
    print("--------------------------")
    print(f"{clustering_name} stats:")
    print("--------------------------")
    for i, stats in enumerate(clustering_stats):
        print([graph.vs[k]['name'] for k in clustering[i]])
        print(f"Cluster {i}: size={stats['cluster_size']}, avg_p={stats['mean_p_value']:.30f}, "
              f"median_p={stats['median_p_value']:.30f}, min_p={stats['min_p_value']:.30f}")
    return clustering_stats

def plot_based_on_pvals(g, cluster_stats, clusters, key, clustering_name, highlight_genes=None):
    if highlight_genes is None:
        # DEMENTIA
        highlight_genes = ["NECTIN2", "TOMM40", "APOE", "APOC1"]
        # highlight_genes = ["NECTIN2", "LNCOB1", "TOMM40", "APOE", "APOC1"]

        # POLIOMYELITIS
        # highlight_genes = ["BOLA3",  "MTHFD2", "TBC1D12", "ARPIN-AP3S2", "ARPIN", "ZNF600", "ZNF702P"] 

    # Map mean_p_values to colors
    mean_p_values = [s[key] for s in cluster_stats]
    
    # Ensure vertex order matches g.vs
    vertex_color_value = [0]*len(g.vs)
    vertex_cluster_value = [0]*len(g.vs)
    for cluster_idx, cluster in enumerate(clusters):
        for v in cluster:
            vertex_color_value[v] = mean_p_values[cluster_idx]
            vertex_cluster_value[v] = cluster_idx

    # Layout coordinates
    random.seed(42)
    layout = g.layout("fr")
    Xn = [coords[0] for coords in layout.coords]
    Yn = [coords[1] for coords in layout.coords]

    # Edges
    Xe = []
    Ye = []
    for e in g.es:
        source, target = e.tuple
        Xe += [layout.coords[source][0], layout.coords[target][0], None]
        Ye += [layout.coords[source][1], layout.coords[target][1], None]

    # Edge trace
    edge_trace = go.Scatter(
        x=Xe,
        y=Ye,
        mode='lines',
        line=dict(color='gray', width=1),
        hoverinfo='none'
    )

    # Generate distinct border colors for clusters
    num_clusters = len(clusters)
    cluster_border_colors = [f"rgb{tuple(int(255*c) for c in cm.tab20(i % 20)[:3])}" for i in range(num_clusters)]

    # Node borders
    line_colors = []
    line_widths = []
    vertex_size = []
    for v in range(len(g.vs)):
        name = g.vs[v]['name']
        cluster_idx = vertex_cluster_value[v]
        if name in highlight_genes:
            line_colors.append('red')
            line_widths.append(5)
            vertex_size.append(30)
        else:
            line_colors.append(cluster_border_colors[cluster_idx])
            line_widths.append(3)
            vertex_size.append(20) 

    # Node trace with colorbar
    node_trace = go.Scatter(
        x=Xn,
        y=Yn,
        mode='markers',
        marker=dict(
            size=vertex_size,
            color=vertex_color_value,       # continuous color based on mean_p_value
            colorscale='Viridis',           # colormap
            colorbar=dict(title=key),
            line=dict(color=line_colors, width=line_widths)
        ),
        text=[f"Node: {g.vs[v]['name']}<br>Cluster: {vertex_cluster_value[v]}<br>{key}: {vertex_color_value[v]:.30f}" 
              for v in range(len(g.vs))],
        hoverinfo='text'
    )

    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title=dict(
            text=f"Graph with clusters colored by {key} ({clustering_name})",
            font=dict(size=28)
        ),
        showlegend=False,
        
        hovermode='closest',
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
    )

    fig.show()

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Cluster a gene graph and calculate p-value statistics per cluster.")
    parser.add_argument("graph", type=str, help="Path to the GML graph file")
    parser.add_argument("pvalues", type=str, help="Path to the CSV file containing gene p-values")
    args = parser.parse_args()

    # Read graph
    graph = read_graph_from_gml(args.graph)
    print("Graph loaded. Number of vertices:", graph.vcount(), "Number of edges:", graph.ecount())

    # Read gene p-values
    gene_p_values_df = read_gene_p_values(args.pvalues)
    print("Gene p-values loaded. Number of genes:", len(gene_p_values_df))

    voronoi_clusters = community_voronoi(graph, weights=graph.es['weight'])
    print(voronoi_clusters)
    print("////////////////////////////////")
    # Run Voronoi clustering
    # voronoi_clusters = voronoi_clustering(graph)
    # print("Voronoi clustering done. Number of clusters:", len(voronoi_clusters))
    # Run Infomap clustering
    infomap_clusters = infomap_clustering(graph)
    print(infomap_clusters)
    # exit()

    # print("Infomap clustering done. Number of clusters:", len(infomap_clusters))

    # print("===============================================")
    # print("============== mean pvalue ====================")
    # print("===============================================")
    # # Calculate cluster p-value stats for Voronoi
    # mean_pv = calculate_stats_main(graph, voronoi_clusters, gene_p_values_df, "mean_pvalue_all", clustering_name="Voronoi clustering")
    # print()
    # print()
    # print()
    # mean_pi = calculate_stats_main(graph, infomap_clusters, gene_p_values_df, "mean_pvalue_all", clustering_name="Infomap clustering")
    

    # print("===============================================")
    # print("============== min pvalue ====================")
    # print("===============================================")
    # # Calculate cluster p-value stats for Voronoi
    # min_pv = calculate_stats_main(graph, voronoi_clusters, gene_p_values_df, "min_pvalue_all", clustering_name="Voronoi clustering")
    # print()
    # print()
    # print()
    # min_pi = calculate_stats_main(graph, infomap_clusters, gene_p_values_df, "min_pvalue_all", clustering_name="Infomap clustering")
    
    print("===============================================")
    print("============== median pvalue ====================")
    print("===============================================")
    # Calculate cluster p-value stats for Voronoi
    median_pv = calculate_stats_main(graph, voronoi_clusters, gene_p_values_df, "median_pvalue_all", clustering_name="Voronoi clustering")
    print()
    print()
    print()
    median_pi = calculate_stats_main(graph, infomap_clusters, gene_p_values_df, "median_pvalue_all", clustering_name="Infomap clustering")
    

    energy_type = 'closeness'
    # plot_based_on_pvals(graph, mean_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (mean_pv)")
    # plot_based_on_pvals(graph, mean_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (mean_pi)")

    # plot_based_on_pvals(graph, min_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (min_pv)")
    # plot_based_on_pvals(graph, min_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (min_pi)")

    plot_based_on_pvals(graph, median_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (median_pv)")
    plot_based_on_pvals(graph, median_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (median_pi)")

    # energy_type = 'betweenness'
    # plot_based_on_pvals(graph, mean_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (mean_pv)")
    # plot_based_on_pvals(graph, mean_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (mean_pi)")

    # # plot_based_on_pvals(graph, min_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (min_pv)")
    # # plot_based_on_pvals(graph, min_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (min_pi)")

    # plot_based_on_pvals(graph, median_pv, voronoi_clusters, f"energy_{energy_type}", "Voronoi Clustering (median_pv)")
    # plot_based_on_pvals(graph, median_pi, infomap_clusters, f"energy_{energy_type}", "Infomap Clustering (median_pi)")





    # plot_based_on_pvals(graph, mean_pv, voronoi_clusters, "mean_p_value", "Voronoi Clustering (mean_pv)")
    # plot_based_on_pvals(graph, mean_pi, infomap_clusters, "mean_p_value", "Infomap Clustering (mean_pi)")

    # plot_based_on_pvals(graph, mean_pv, voronoi_clusters, "min_p_value", "Voronoi Clustering (mean_pv)")
    # plot_based_on_pvals(graph, mean_pi, infomap_clusters, "min_p_value", "Infomap Clustering (mean_pi)")

    # plot_based_on_pvals(graph, mean_pv, voronoi_clusters, "median_p_value", "Voronoi Clustering (mean_pv)")
    # plot_based_on_pvals(graph, mean_pi, infomap_clusters, "median_p_value", "Infomap Clustering (mean_pi)")

    # plot_based_on_pvals(graph, min_pv, voronoi_clusters, "mean_p_value", "Voronoi Clustering (min_pv)")
    # plot_based_on_pvals(graph, min_pi, infomap_clusters, "mean_p_value", "Infomap Clustering (min_pi)")

    # plot_based_on_pvals(graph, min_pv, voronoi_clusters, "min_p_value", "Voronoi Clustering (min_pv)")
    # plot_based_on_pvals(graph, min_pi, infomap_clusters, "min_p_value", "Infomap Clustering (min_pi)")

    # plot_based_on_pvals(graph, min_pv, voronoi_clusters, "median_p_value", "Voronoi Clustering (min_pv)")
    # plot_based_on_pvals(graph, min_pi, infomap_clusters, "median_p_value", "Infomap Clustering (min_pi)")

    # plot_based_on_pvals(graph, median_pv, voronoi_clusters, "mean_p_value", "Voronoi Clustering (median_pv)")
    # plot_based_on_pvals(graph, median_pi, infomap_clusters, "mean_p_value", "Infomap Clustering (median_pi)")

    # plot_based_on_pvals(graph, median_pv, voronoi_clusters, "min_p_value", "Voronoi Clustering (median_pv)")
    # plot_based_on_pvals(graph, median_pi, infomap_clusters, "min_p_value", "Infomap Clustering (median_pi)")

    # plot_based_on_pvals(graph, median_pv, voronoi_clusters, "median_p_value", "Voronoi Clustering (median_pv)")
    # plot_based_on_pvals(graph, median_pi, infomap_clusters, "median_p_value", "Infomap Clustering (median_pi)")
