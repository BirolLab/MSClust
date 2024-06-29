import argparse
import collections

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
import sklearn.neighbors
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", help="path to distance matrix file", required=True)
    parser.add_argument("-y", help="path to cluster labels file", required=True)
    parser.add_argument("-k", help="number of nearest neighbors", type=int, default=3)
    parser.add_argument(
        "-j", help="number of cpus, use -1 for all", type=int, default=1
    )
    parser.add_argument("-o", help="output prefix", required=True)
    return parser.parse_args()


def get_model(distances, num_neighbors, num_jobs):
    knn = sklearn.neighbors.NearestNeighbors(
        n_neighbors=num_neighbors,
        n_jobs=num_jobs,
        metric="precomputed",
    )
    knn.fit(distances)
    return knn


def plot_graph(model, labels, out_path):
    cmap = sns.husl_palette(20)
    nx_graph = nx.from_scipy_sparse_array(model.kneighbors_graph(mode="distance"))
    _, ax = plt.subplots(figsize=(10, 10), dpi=350)
    ax.axis("off")
    labels_index = dict(zip(sorted(list(set(labels))), range(len(labels))))
    node_colors = [cmap[labels_index[label]] for label in labels]
    edge_weights = [nx_graph.edges[u, v]["weight"] for u, v in nx_graph.edges()]
    nx.draw_networkx(
        nx_graph,
        ax=ax,
        node_color=node_colors,
        with_labels=False,
        node_size=50,
        edge_color=edge_weights,
        edge_cmap=plt.cm.Greys_r,
        edge_vmin=min(edge_weights),
        edge_vmax=max(edge_weights),
    )
    handles = []
    legend_labels = set()
    for label in labels:
        color = cmap[labels_index[label]]
        if label not in legend_labels:
            handles.append(mpatches.Patch(color=color, label=label))
        legend_labels.add(label)
    plt.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(legend_labels) // 5,
    )
    plt.tight_layout()
    plt.savefig(out_path + ".png")


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    x, labels, _ = utils.load_data(args.x, args.y, None)
    labels = list(map(utils.shorten_label, labels))
    print("done")
    print("constructing knn graph...", end=" ", flush=True)
    model = get_model(x, args.k, args.j)
    print("done")
    print("plotting graph...", end=" ", flush=True)
    plot_graph(model, labels, args.o)
    print("done")


if __name__ == "__main__":
    main()
