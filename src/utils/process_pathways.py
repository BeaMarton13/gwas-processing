import pandas as pd
import os
import sys
import collections
from functools import reduce

def _get_intersection_and_counts(list_of_lists):
    """
    Calculates the intersection of a list of lists and the total frequency 
    of each intersecting element across all original lists.
    """
    if not list_of_lists:
        return set(), {}

    # 1. Get the intersection of all lists (treating them as sets)
    # Convert each sublist to a set, then use reduce to find the common elements.
    sets = [set(lst) for lst in list_of_lists]
    intersection_set = reduce(set.intersection, sets)

    # 2. Flatten the list of lists to count all elements
    all_elements = [item for sublist in list_of_lists for item in sublist]
    full_counts = collections.Counter(all_elements)

    # 3. Filter the counts to only include elements in the intersection set
    intersection_counts = {
        element: full_counts[element] 
        for element in intersection_set
    }
    
    return intersection_set, intersection_counts


def read_data(file_w_path):
    df = pd.read_csv(file_w_path)
    return df

def cut_by_column(df, col, min_val):
    return df[df[col] > min_val]

def get_pathways(df, col):
    return df[col]

def _iterate_through_directory(directory_path):
    csv_files = []
    for file in os.listdir(directory_path):
        if file.endswith('.csv'):
            csv_files.append(os.path.join(directory_path, file))

    return csv_files

def calculate_column(df, target_value, target_col='Pathway', calc_col='Fold Enrichment'):
    filtered_df = df[df[target_col] == target_value]

    if filtered_df.empty:
        return None
    
    return filtered_df[calc_col].values[0]

def process_pathways(background_type):
    directory_path = os.path.join(sys.path[0], f'../../data/shinygo/{background_type}')
    csv_files = _iterate_through_directory(directory_path)
    list_of_sets = []
    list_of_dfs = []
    for file in csv_files:
        df = read_data(file)
        # NOTE Do some cutting as needed
        df = cut_by_column(df, 'Enrichment FDR', 0.05)
        # df = cut_by_column(df, 'Fold Enrichment', 0.5)
        list_of_dfs.append(df)
        list_of_sets.append(set(get_pathways(df, 'Pathway')))
        # list_of_sets.append(list(get_pathways(df, 'Pathway')))

    common_elements = reduce(set.intersection, list_of_sets)
    # common_elements, common_count = _get_intersection_and_counts(list_of_sets)
    # print('===================================================')
    # print('===================================================')
    # print(common_count)
    # print('===================================================')
    # print('===================================================')

    avg_fdrs = {}
    for elem in common_elements:
        avg_fdrs[elem] = sum(calculate_column(df, elem, target_col='Pathway', calc_col='Enrichment FDR') for df in list_of_dfs) / len(list_of_dfs)
        # avg_fdrs[elem] = sum(calculate_column(df, elem, target_col='Pathway', calc_col='Fold Enrichment') for df in list_of_dfs) / len(list_of_dfs)
        
        # print(elem)
        # print('Avg FDR: ', avg_fdrs[elem])

    sorted_avg_fdrs = sorted(avg_fdrs.items(), key=lambda item: item[1], reverse=True)
    for pathway, avg_fdr in sorted_avg_fdrs:
        print(f'Average FDR: {avg_fdr}\t\tPathway: {pathway}')


if __name__ == '__main__':
    background_type = ['biological_process', 'cellular_component', 'molecular_function']
    for bt in background_type:
        print('====================================================')
        print('====================================================')
        print('====================================================')
        print(f'============ {bt.upper()} ============')
        print('====================================================')
        process_pathways(bt)

    

    
