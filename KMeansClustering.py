import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
sample_list = []

# this code uses Kmeans clustering to cluster all the samples on the non-binary matrix of euclidian distances

NUM_CLUSTERS = 6
# global variable to define how many clusters to group samples into


def create_aggregated_grayscales(matrix_file):
    matrix_df = pd.read_csv(matrix_file)
    data_list = []
    for index, row in matrix_df.iterrows():
        if row['sample'] == 'blank':
            continue
        sample_list.append(row['sample'])
        gray_scale_img = []  # 5x5 array of grayscale image values
        for i in range(5):
            gray_scale_img.append([])  # add an array for a row
            for j in range(5):
                gray_scale_img[i].append(row[f"({i},{j})"])  # unsure about format
        data_list.append(gray_scale_img)
    return data_list


def perform_clustering(grayscale_images_list):
    grayscale_images = np.array(grayscale_images_list)
    flattened_images = grayscale_images.reshape(grayscale_images.shape[0], -1)
    n_init_value = 50
    kmeans = KMeans(n_clusters=NUM_CLUSTERS, n_init=n_init_value, random_state=42)
    cluster_labels = kmeans.fit_predict(flattened_images)
    return cluster_labels

def produce_bar_plot(flattened_images):
    n_init_value = 50
    kmeans = KMeans(n_clusters=NUM_CLUSTERS, n_init=n_init_value, random_state=42)
    cluster_labels = kmeans.fit_predict(flattened_images)
    # Count the number of samples in each cluster
    cluster_counts = np.bincount(cluster_labels)

    # Plot the histogram or bar plot of cluster assignments
    plt.bar(range(NUM_CLUSTERS), cluster_counts, tick_label=[f"{i}" for i in range(NUM_CLUSTERS)])
    plt.xlabel("Cluster")
    plt.ylabel("Number of Samples")
    plt.title("Sample Distribution across Clusters")
    plt.show()


def elbow_method(flattened_images):
    # Perform K-means for different numbers of clusters and calculate WCSS
    wcss = []
    max_clusters = 10
    for NUM_CLUSTERS in range(1, max_clusters + 1):
        kmeans = KMeans(n_clusters=NUM_CLUSTERS, random_state=42)
        kmeans.fit(flattened_images)
        wcss.append(kmeans.inertia_)

    # Plot the Elbow curve
    plt.plot(range(1, max_clusters + 1), wcss)
    plt.xlabel('Number of Clusters')
    plt.ylabel('WCSS')
    plt.title('Elbow Method to Find Optimal Clusters')
    plt.show()


def score(flattened_images):
    # Calculate Silhouette scores for different numbers of clusters
    max_clusters = 10
    silhouette_scores = []
    for NUM_CLUSTERS in range(2, max_clusters + 1):
        kmeans = KMeans(n_clusters=NUM_CLUSTERS, random_state=42)
        cluster_labels = kmeans.fit_predict(flattened_images)
        silhouette_scores.append(silhouette_score(flattened_images, cluster_labels))

    # Plot the Silhouette scores
    plt.plot(range(2, max_clusters + 1), silhouette_scores)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Silhouette Score')
    plt.title('Silhouette Score to Find Optimal Clusters')
    plt.show()

def return_file(output_file, cluster_labels):
    data_frame = pd.DataFrame({
        'sample': sample_list,
        'clust': cluster_labels,
    })
    data_frame.to_csv(output_file, index=False)
def main():
    matrix_file_name = "matrix.csv"
    output_file_name = "KMeansClusters.csv"
    grayscale_images_list = create_aggregated_grayscales(matrix_file_name)
    cluster_labels = perform_clustering(grayscale_images_list)
    return_file(output_file_name, cluster_labels)


main()
