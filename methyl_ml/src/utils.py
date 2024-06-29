import re

import numpy as np


def load_data(
    vectors_path: str,
    labels_path: str | None = None,
    clusters_path: str | None = None,
):
    vectors = np.load(vectors_path)
    if labels_path:
        with open(labels_path) as f:
            labels = f.read().split()
    else:
        labels = None
    if clusters_path:
        clusters = np.load(clusters_path)
    else:
        clusters = None
    return vectors, labels, clusters


def shorten_label(label: str):
    return re.search(r"GSM\d+_(.+?)-Z0+", label).group(1)


def index_labels(labels: list[str]):
    return dict(zip(set(labels), range(len(labels))))
