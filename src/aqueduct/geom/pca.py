import numpy as np
import scipy.linalg


class Center(object):
    def __init__(self, X):
        self.mean = np.mean(X, 0)

    def __call__(self, X):
        return X - self.mean

    def undo(self, X):
        return X + self.mean


class Normalize(object):
    def __init__(self, X):
        self.std = np.std(X, 0)

    def __call__(self, X):
        return X / self.std

    def undo(self, X):
        return X * self.std


class Standartize(Center, Normalize):
    def __init__(self, X):
        Center.__init__(self, X)
        Normalize.__init__(self, X)

    def __call__(self, X):
        return Normalize.__call__(self, Center.__call__(self, X))

    def undo(self, X):
        return Center.undo(self, Normalize.undo(self.X))


class PCA(object):
    def __init__(self, X, prepro=None):
        self.preprocess_method = prepro
        # SVD
        self.U, self.d, self.Pt = scipy.linalg.svd(self.preprocess(X), full_matrices=False)
        assert np.all(self.d[:-1] >= self.d[1:]), "SVD error."  # sorted
        # PCA
        self.T = self.U * self.d
        self.eigen = self.d ** 2
        self.sumvariance = np.cumsum(self.eigen)
        self.sumvariance /= self.sumvariance[-1]

    @property
    def P(self):
        return self.Pt.T

    def preprocess(self, X):
        if self.preprocess_method is None:
            return X
        return self.preprocess_method(X)

    def preprocess_undo(self, X):
        if self.preprocess_method is None:
            return X
        return self.preprocess_method.undo(X)

    def __call__(self, X):
        return np.dot(self.preprocess(X), self.P)

    def undo(self, T):
        return self.preprocess_undo(np.dot(T, self.Pt))
