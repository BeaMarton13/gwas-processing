import pandas as pd
import argparse
import os
import igraph as ig

# def get_ig_graph_from_tsv(filenane, weight_col):


def argument_handler():
    parser = argparse.ArgumentParser(
        description="Creates an igraph graph from TSV file."
    )
    parser.add_argument(
        "file_w_path",
        type=str,
        help="filename of the TSV file containing STRING interactions with path (e.g., gene_networks/dementia001.tsv)"
    )
    args = parser.parse_args()
    return args

def read_tsv_as_df(filename):# Load the TSV file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_w_path = os.path.abspath(os.path.join(script_dir, f'../../{filename}'))
    return pd.read_csv(file_w_path, sep='\t')

def create_igraph_from_df(df, weight_col):
    # Create an igraph Graph from the DataFrame
    g = ig.Graph.DataFrame(df[['node1', 'node2', weight_col]], directed=True, use_vids=False)
    
    # Assign weights to edges based on the specified weight column
    if weight_col in df.columns:
        g.es['weight'] = df[weight_col].tolist()
    else:
        print(f"Warning: Weight column '{weight_col}' not found in DataFrame columns.")
    
    return g

def analyze_graph(g):
    print("\n--- Graph Attributes Check ---")
    print(f"Total Vertices (Proteins): {g.vcount()}")
    print(f"Total Edges (Interactions): {g.ecount()}")

    # Check if weights were successfully assigned
    if 'weight' in g.edge_attributes():
        print(f"Edge weights assigned successfully using column '{weight_columns[0]}'.")
        print(f"Min Weight: {min(g.es['weight']):.3f}")
        print(f"Max Weight: {max(g.es['weight']):.3f}")
    else:
        print("Warning: Edge weights were not found. Check the column name.")

    # 2. Calculate and display the degree of the top 5 most connected proteins
    # The degree() function returns a list of degrees for all nodes
    degrees = g.degree()

    # Add degree as a vertex attribute
    g.vs['degree'] = degrees

    # Create a list of tuples (protein_name, degree) and sort it
    protein_degrees = [(v['name'], v['degree']) for v in g.vs]
    protein_degrees.sort(key=lambda x: x[1], reverse=True)

    print("\n--- Top 5 Hub Proteins (by Degree) ---")
    for name, degree in protein_degrees[:5]:
        print(f"Protein: {name:<10} | Degree: {degree}")


def export_graph_to_gml(g, filename):
    gml_file = filename.replace('.tsv', '.gml')
    g.write_gml(gml_file)
    print(f"\nGraph exported to GML format at: {gml_file}")



if __name__ == '__main__':
    args = argument_handler()
    filename = args.file_w_path.strip()

    df = read_tsv_as_df(filename)

    # Clean the column names by removing any leading '#'
    df.columns = df.columns.str.replace('#', '')

    # Define the columns that represent the network weights
    # weight_columns = [
    #     'neighborhood_on_chromosome', 'gene_fusion', 'phylogenetic_cooccurrence',
    #     'homology', 'coexpression', 'experimentally_determined_interaction',
    #     'database_annotated', 'automated_textmining', 'combined_score'
    # ]
    weight_columns = [
        'combined_score'
    ]
    df = df[['node1', 'node2'] + weight_columns]

    g = create_igraph_from_df(df, weight_columns[0])

    export_graph_to_gml(g, filename)

    print(g.community_infomap(edge_weights=g.es['weight']))

    