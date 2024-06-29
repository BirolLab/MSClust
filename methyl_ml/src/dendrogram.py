import argparse

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy
import scipy.spatial.distance
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-x", help="path to vectors/distance matrix file", required=True
    )
    parser.add_argument("-y", help="path to labels file")
    parser.add_argument("-o", help="output path")
    return parser.parse_args()


def plot(linkage_data, labels, out_path):
    _, ax = plt.subplots(figsize=(20, 10), dpi=350)
    scipy.cluster.hierarchy.dendrogram(
        linkage_data,
        labels=labels,
        leaf_rotation=90,
    )
    labels_index = dict(zip(set(labels), range(len(labels))))
    handles = []
    legend_labels = set()
    for lbl in ax.get_xmajorticklabels():
        label = lbl.get_text()
        color = plt.cm.tab20.colors[labels_index[label]]
        if label not in legend_labels:
            handles.append(mpatches.Patch(color=color, label=label))
        lbl.set_color(color)
        legend_labels.add(label)
    plt.tight_layout()
    ax.spines[["right", "top"]].set_visible(False)
    plt.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(legend_labels) // 4,
    )
    plt.savefig(out_path)


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    x, labels, _ = utils.load_data(args.x, args.y, None)
    labels = list(map(utils.shorten_label, labels))
    print("done")
    if x.shape[0] == x.shape[1]:
        print("converting distance matrix...", end=" ", flush=True)
        np.fill_diagonal(x, 0)
        x = scipy.spatial.distance.squareform(x)
        print("done")
    print("getting linkage matrix...", end=" ", flush=True)
    x = scipy.cluster.hierarchy.linkage(
        x,
        method="complete",
        metric="cosine",
    )
    print("done")
    print("plotting results...", end=" ", flush=True)
    plot(x, labels, args.o)
    print("done")


if __name__ == "__main__":
    main()
