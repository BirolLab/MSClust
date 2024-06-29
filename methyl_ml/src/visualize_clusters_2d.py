import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.decomposition
import sklearn.manifold
import umap
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", help="path to vectors file", required=True)
    parser.add_argument("-y", help="path to labels file", required=True)
    parser.add_argument("-c", help="path to clustering results")
    parser.add_argument(
        "-a",
        help="dimensionality reduction algorithm",
        choices=["pca", "tsne", "umap", "mds", "isomap"],
        default="pca",
    )
    parser.add_argument("-o", help="output path")
    parser.add_argument(
        "-j", help="number of cpus, use -1 for all", type=int, default=1
    )
    return parser.parse_args()


def get_2d_vectors(vectors, algorithm, num_jobs):
    if algorithm == "pca":
        return sklearn.decomposition.PCA(2).fit_transform(vectors)
    elif algorithm == "tsne":
        return sklearn.manifold.TSNE(2, n_jobs=num_jobs).fit_transform(vectors)
    elif algorithm == "umap":
        model = umap.UMAP(n_components=2, metric="cosine", n_jobs=num_jobs)
        return model.fit_transform(vectors)
    elif algorithm == "mds":
        return sklearn.manifold.MDS(2, n_jobs=num_jobs).fit_transform(vectors)
    elif algorithm == "isomap":
        model = sklearn.manifold.Isomap(n_jobs=num_jobs, metric="cosine")
        return model.fit_transform(vectors)
    else:
        raise ValueError("unsupported algorithm")


def plot(vectors_2d, labels, clusters, title, out_path):
    plt.style.use("seaborn-v0_8")
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=plt.cm.tab20.colors)
    df = pd.DataFrame(dict(x=vectors_2d[:, 0], y=vectors_2d[:, 1], label=labels))
    groups = df.groupby("label")
    plt.figure(figsize=(12, 10), dpi=350)
    plt.title(title)
    for name, group in groups:
        plt.scatter(group.x, group.y, label=name, edgecolor="k", s=100)
    if clusters is not None:
        for i, txt in enumerate(clusters):
            plt.annotate(txt, (vectors_2d[i, 0], vectors_2d[i, 1]))
    plt.legend(
        frameon=True,
        framealpha=0.9,
        facecolor="inherit",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
    )
    plt.tight_layout()
    plt.savefig(out_path)


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    vectors, labels, clusters = utils.load_data(args.x, args.y, args.c)
    labels = list(map(utils.shorten_label, labels))
    print("done")
    if vectors.shape[1] != 2:
        print("reducing dimensionality...", end=" ", flush=True)
        vectors = get_2d_vectors(vectors, args.a, args.j)
        np.save(args.o + ".npy", vectors)
        print("done")
    print("plotting results...", end=" ", flush=True)
    plot(vectors, labels, clusters, args.a.upper(), args.o)
    print("done")


if __name__ == "__main__":
    main()
