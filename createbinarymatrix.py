import pandas as pd
import math
import csv



header = ['(0,0)', '(0,1)', '(0,2)', '(0,3)', '(0,4)', '(1,0)', '(1,1)', '(1,2)', '(1,3)', '(1,4)', '(2,0)', '(2,1)',
          '(2,2)', '(2,3)', '(2,4)', '(3,0)', '(3,1)', '(3,2)', '(3,3)', '(3,4)', '(4,0)', '(4,1)', '(4,2)', '(4,3)',
          '(4,4)']
smallest_dict = {}
def create_binary_matrix(matrix_file, output_file):
    matrix_df = pd.read_csv(matrix_file)
    for index, row in matrix_df.iterrows():
        which_sample = row['sample']
        if which_sample == 'blank':
            continue
        smallest = 2
        smallest_cluster = ""
        for cluster in header:
            curr_value = row[cluster]
            if (curr_value < smallest):
                smallest = curr_value
                smallest_cluster = cluster
        smallest_dict[which_sample] = smallest_cluster
    adjust_values(matrix_df, output_file)

def adjust_values(matrix_df, output_file):
    for index, row in matrix_df.iterrows():
        which_sample = row['sample']
        if which_sample in smallest_dict:
            smallest_cluster = smallest_dict[which_sample]
            for cluster in header:
                if cluster == smallest_cluster:
                    matrix_df.at[index, cluster] = 1
                else:
                    matrix_df.at[index, cluster] = 0
    matrix_df.to_csv(output_file, index=False)


def create_clustered_file():
    samples_array = []
    clusters = []
    for sample in smallest_dict:
        samples_array.append(sample)
        index = header.index(smallest_dict[sample])
        clusters.append(index)

    result_df = pd.DataFrame({
        'Sample': samples_array,
        'clust': clusters,
    })

    filename = "sample_closest_cluster.csv"
    result_df.to_csv(filename, index=False)






def main():
    matrix_file = "matrix.csv"
    output_file = "binary_matrix.csv"
    create_binary_matrix(matrix_file, output_file)
    create_clustered_file()

main()