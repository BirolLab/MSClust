import argparse

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sklearn.metrics
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", help="path to vectors file", required=True)
    parser.add_argument("-y", help="path to labels file", required=True)
    parser.add_argument("-o", help="output path")
    parser.add_argument(
        "-j", help="number of cpus, use -1 for all", type=int, default=1
    )
    return parser.parse_args()


def plot(dist, labels, out_path):
    mask = np.zeros_like(dist, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    plt.figure(figsize=(8, 8), dpi=350)
    cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
    mid_labels = []
    for i, label in enumerate(labels):
        if i % 20 == 10:
            mid_labels.append(label)
        else:
            mid_labels.append("")
    sns.heatmap(
        dist,
        mask=mask,
        xticklabels=mid_labels,
        yticklabels=mid_labels,
        cmap=cmap,
        cbar_kws={"label": "Distance"},
    )
    plt.title("Cosine Distance Matrix")
    plt.tight_layout()
    plt.savefig(out_path)


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    vectors, labels, _ = utils.load_data(args.x, args.y, None)
    labels = list(map(utils.shorten_label, labels))
    print("done")
    if vectors.shape[0] != vectors.shape[1]:
        print("constructing matrix...", end=" ", flush=True)
        vectors = sklearn.metrics.pairwise_distances(
            vectors,
            metric="cosine",
            n_jobs=args.j,
        )
        np.save(args.o + ".npy", vectors)
        print("done")
    print("plotting results...", end=" ", flush=True)
    plot(vectors, labels, args.o)
    print("done")


if __name__ == "__main__":
    main()
