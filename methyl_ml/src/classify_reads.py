import argparse

import numpy as np
import sklearn.metrics
import sklearn.neighbors
import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-xtrain", help="path to distance matrix file", required=True)
    parser.add_argument("-ytrain", help="path to cluster labels file", required=True)
    parser.add_argument("-xtest", help="path to query distances file", required=True)
    parser.add_argument("-ytest", help="path to query ground truth file", required=True)
    parser.add_argument("-k", help="number of nearest neighbors", type=int, default=3)
    parser.add_argument(
        "-j", help="number of cpus, use -1 for all", type=int, default=1
    )
    return parser.parse_args()


def get_model(distances, clusters, num_neighbors, num_jobs):
    knn = sklearn.neighbors.KNeighborsClassifier(
        n_neighbors=num_neighbors,
        weights="distance",
        metric="precomputed",
        n_jobs=num_jobs,
    )
    knn.fit(distances, clusters)
    return knn


def predict(model, x_test):
    unassigned = np.argwhere(np.isnan(x_test).all(axis=1))
    x_test[unassigned, :] = 0
    np.nan_to_num(x_test, copy=False, nan=1)
    labels = model.predict(x_test)
    return labels, unassigned


def report_results(y_pred, y_true):
    print(sklearn.metrics.classification_report(y_true, y_pred))


def main():
    args = get_args()
    print("loading data...", end=" ", flush=True)
    x_train, y_train, _ = utils.load_data(args.xtrain, args.ytrain, None)
    x_test, y_test, _ = utils.load_data(args.xtest, args.ytest, None)
    x_test[x_test < 0] = 1
    y_train = list(map(utils.shorten_label, y_train))
    print("done")
    print("constructing knn graph...", end=" ", flush=True)
    model = get_model(x_train, y_train, args.k, args.j)
    print("done")
    print("querying reads...", end=" ", flush=True)
    y_pred, unassigned = predict(model, x_test)
    y_pred[unassigned] = "unassigned"
    print("done")
    print("unassigned: ", *unassigned.reshape(-1))
    print("results:", end=" ", flush=True)
    report_results(y_pred, y_test)


if __name__ == "__main__":
    main()
