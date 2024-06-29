import argparse

import numpy as np
import sklearn.cluster
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", help="path to vectors file", required=True)
    parser.add_argument("-y", help="path to labels file")
    parser.add_argument("-o", help="output path")
    parser.add_argument(
        "-a", help="clustering algorithm", choices=["agg", "kmeans"], default="agg"
    )
    parser.add_argument("-c", help="number of clusters", type=int)
    return parser.parse_args()


def cluster(vectors, num_clusters, algorithm):
    if algorithm == "agg" and not num_clusters:
        return sklearn.cluster.AgglomerativeClustering(
            distance_threshold=0,
            n_clusters=None,
        ).fit_predict(vectors)
    elif algorithm == "agg" and num_clusters:
        model = sklearn.cluster.AgglomerativeClustering(num_clusters)
        return model.fit_predict(vectors)
    elif algorithm == "kmeans" and num_clusters:
        return sklearn.cluster.KMeans(num_clusters).fit_predict(vectors)
    else:
        raise ValueError("unsupported algorithm/num_clusters")


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    vectors, labels, _ = utils.load_data(args.x, args.y, None)
    print("done")
    print("clustering...", end=" ", flush=True)
    clusters = cluster(vectors, args.c, args.a)
    print("done")
    if labels:
        label_index = utils.index_labels(labels)
        y_true = list(map(lambda x: label_index[x], labels))
        print("ARI = ", sklearn.metrics.adjusted_rand_score(y_true, clusters))
        print("AMI = ", sklearn.metrics.adjusted_mutual_info_score(y_true, clusters))
    print("saving results...", end=" ", flush=True)
    np.save(args.o, clusters)
    print("done")


if __name__ == "__main__":
    main()
