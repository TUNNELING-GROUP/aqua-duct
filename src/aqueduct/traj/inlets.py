# -*- coding: utf-8 -*-
import numpy as np
from collections import namedtuple
from scipy.spatial.distance import pdist
import copy

from aqueduct.utils.helpers import is_iterable, listify, lind


class ProtoInletTypeCodes:
    surface = 'surface'
    internal = 'internal'

    incoming = 'inin'
    outgoing = 'inout'


class InletTypeCodes(ProtoInletTypeCodes):
    # TODO: write it in a more smart way

    all_surface = [(ProtoInletTypeCodes.surface, itype) for itype in
                   (ProtoInletTypeCodes.incoming, ProtoInletTypeCodes.outgoing)]
    all_internal = [(ProtoInletTypeCodes.internal, itype) for itype in
                    (ProtoInletTypeCodes.incoming, ProtoInletTypeCodes.outgoing)]
    all_incoming = [(itype, ProtoInletTypeCodes.incoming) for itype in
                    (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.internal)]
    all_outgoing = [(itype, ProtoInletTypeCodes.outgoing) for itype in
                    (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.internal)]

    surface_incoming = (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.incoming)
    internal_incoming = (ProtoInletTypeCodes.internal, ProtoInletTypeCodes.incoming)
    internal_outgoing = (ProtoInletTypeCodes.internal, ProtoInletTypeCodes.outgoing)
    surface_outgoing = (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.outgoing)


# clusers can be:
# None  - no such cluster
# nr - 0 means outliers


class InletClusterGenericType(object):
    def __init__(self, inp, out):
        self.clusters = [inp, out]

    @property
    def input(self):
        return self.clusters[0]

    @property
    def output(self):
        return self.clusters[-1]

    def cluster2str(self, cl):
        if cl is None:
            return 'N'
        return '%d' % int(cl)

    def __getitem__(self, item):
        return self.clusters[item]

    def __len__(self):
        return len(self.clusters)

    def __str__(self):
        return ':'.join(map(self.cluster2str, list(self)[::2] + list(self)[::-2]))

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, str(self))

    def make_val(self, base):
        val = 0.
        for nr, e in enumerate(tuple(self)[::-1]):
            if e is None:
                e = float(base) ** (nr + 1)
            else:
                e += 1
            for dummy_iter in range(nr + 1):
                e /= float(base)
            val += e
        return val

    def __cmp__(self, other):
        if other is None:
            return 1
        if not isinstance(other, self.__class__):
            return 1

        result = 0
        base = max(max(self), max(other), len(self), len(other)) + 2.

        if self.make_val(base) - other.make_val(base) > 0:
            return 1
        elif self.make_val(base) - other.make_val(base) < 0:
            return -1
        return 0

    def __hash__(self):
        return hash(str(self))


class InletClusterExtendedType(InletClusterGenericType):
    def __init__(self, surfin, interin, interout, surfout):
        InletClusterGenericType.__init__(self, surfin, surfout)
        self.clusters.extend([interin, interout])

    @property
    def generic(self):
        return InletClusterGenericType(*self.clusters[:2])


Inlet = namedtuple('Inlet', 'coords type reference')


class Inlets(object):
    # class for list of inlets
    def __init__(self, spaths, onlytype=InletTypeCodes.all_surface):

        self.onlytype = onlytype
        self.inlets_list = []
        self.inlets_ids = []
        self.clusters = []
        self.number_of_clustered_inlets = None

        for spath in spaths:
            self.extend_inlets(spath)

    def extend_inlets(self, spath, onlytype=None):

        if onlytype is None:
            onlytype = self.onlytype

        nr = len(self.inlets_list)
        for inlet in spath.get_inlets():
            if onlytype is not None:
                if inlet.type not in onlytype:
                    continue
            self.inlets_list.append(inlet)
            self.inlets_ids.append(nr)
            nr += 1

    def add_cluster_annotations(self, clusters):
        assert len(clusters) == len(self.inlets_list)
        self.clusters = clusters

    # basic properites

    @property
    def size(self):
        return len(self.inlets_list)

    @property
    def coords(self):
        return [inlet.coords.tolist() for inlet in self.inlets_list]

    @property
    def types(self):
        return [inlet.type for inlet in self.inlets_list]

    @property
    def refs(self):
        return [inlet.reference for inlet in self.inlets_list]

    def perform_clustering(self, method):
        # this do clean clustering, all previous clusters are discarded
        # 0 means outliers
        self.add_cluster_annotations(method(np.array(self.coords)))
        self.number_of_clustered_inlets = len(self.clusters)
        # renumber clusters
        self.renumber_clusters()

    def perform_reclustering(self, method):
        # this do reclusterization of all clusters, if no cluster exists perform_clustering is called
        if len(self.clusters) == 0:
            return self.perform_clustering(method)
        for cluster in self.clusters_list:
            self.recluster_cluster(method,cluster)
        # number of cluster
        self.number_of_clustered_inlets = len(self.clusters)
        # renumber clusters
        self.renumber_clusters()

    def recluster_cluster(self,method,cluster):
        if cluster in self.clusters_list:
            reclust = method(np.array(self.lim2clusters(cluster).coords))
            max_cluster = max(self.clusters_list) + 1
            nrr = 0  # recluster nr
            for nr, c in enumerate(self.clusters):
                if c == cluster:
                    self.clusters[nr] = reclust[nrr]+max_cluster
                    nrr += 1
        # number of cluster
        self.number_of_clustered_inlets = len(self.clusters)

    def recluster_outliers(self, method):
        if 0 in self.clusters_list:
            max_cluster = max(self.clusters_list)
            reclust = method(np.array(self.lim2clusters(0).coords))
            nrr = 0  # recluster nr
            for nr, c in enumerate(self.clusters):
                if c == 0:
                    rc = reclust[nrr]
                    nrr += 1
                    if rc > 0:
                        rc = rc + max_cluster
                    self.clusters[nr] = rc
            # number of cluster
            self.number_of_clustered_inlets = len(self.clusters)
            # renumber clusters
            self.renumber_clusters()

    def small_clusters_to_outliers(self, maxsize):
        for c in self.clusters_list:
            if c == 0:
                continue
            if self.clusters.count(c) <= maxsize:
                for nr, cc in enumerate(self.clusters):
                    if cc == c:
                        self.clusters[nr] = 0
        # renumber clusters
        self.renumber_clusters()

    def renumber_clusters(self):
        if 0 in self.clusters_list:
            new_numbers = range(len(self.clusters_list))
        else:
            new_numbers = range(1, len(self.clusters_list) + 1)
        old_numbers = self.clusters_list
        if old_numbers != new_numbers:
            for nr, c in enumerate(self.clusters):
                self.clusters[nr] = new_numbers[old_numbers.index(c)]
        self.sort_clusters()

    def sort_clusters(self):
        old_numbers = self.clusters_list
        sizes = self.clusters_size
        # now sort according to sizes but put 0, if present, at the begining
        if 0 in old_numbers:
            sizes.pop(old_numbers.index(0))
            old_numbers.pop(old_numbers.index(0))
            new_numbers = [0] + map(lambda i: old_numbers[i],np.argsort(sizes).tolist()[::-1])
        else:
            new_numbers = map(lambda i: old_numbers[i],np.argsort(sizes).tolist()[::-1])
        new_clusters = []
        for c in self.clusters:
            new_clusters.append(new_numbers[self.clusters_list.index(c)])
        self.clusters = new_clusters

    @property
    def clusters_list(self):
        return sorted(list(set(self.clusters)))

    @property
    @listify
    def clusters_centers(self):
        for c in self.clusters_list:
            yield np.mean(self.lim2clusters(c).coords, 0)

    @property
    def clusters_size(self):
        return map(self.clusters.count, self.clusters_list)

    @property
    @listify
    def clusters_std(self):
        for c, s in zip(self.clusters_list, self.clusters_size):
            if s == 1:
                yield 0.
            elif s == 2:
                yield pdist(self.lim2clusters(c).coords, 'euclidean')
            elif s > 2:
                yield np.std(pdist(self.lim2clusters(c).coords, 'euclidean'))

    @listify
    def spaths2ctypes(self, spaths):
        # surfin interin interout surfout
        for sp in spaths:
            ctypes = self.lim2spaths(sp).types
            clusters = self.lim2spaths(sp).clusters
            surfin, interin, interout, surfout = None, None, None, None
            for ct, cl in zip(ctypes, clusters):
                if ct in InletTypeCodes.all_surface and ct in InletTypeCodes.all_incoming:
                    surfin = cl
                if ct in InletTypeCodes.all_internal and ct in InletTypeCodes.all_incoming:
                    interin = cl
                if ct in InletTypeCodes.all_internal and ct in InletTypeCodes.all_outgoing:
                    interout = cl
                if ct in InletTypeCodes.all_surface and ct in InletTypeCodes.all_outgoing:
                    surfout = cl
            yield InletClusterExtendedType(surfin, interin, interout, surfout)

    def lim_to(self, what, towhat):
        if not is_iterable(towhat):
            towhat = [towhat]
        new_inlets = self.__class__([], onlytype=self.onlytype)
        new_inlets.number_of_clustered_inlets = self.number_of_clustered_inlets

        for inlet, ids, cluster, w in zip(self.inlets_list, self.inlets_ids, self.clusters, what):
            if w in towhat:
                new_inlets.inlets_list.append(inlet)
                new_inlets.inlets_ids.append(ids)
                new_inlets.clusters.append(cluster)

        return new_inlets

    def lim2spaths(self, spaths):
        if not is_iterable(spaths):
            spaths = [spaths]
        return self.lim_to(self.refs, [sp.id for sp in spaths])

    def lim2types(self, types):
        return self.lim_to(self.types, types)

    def lim2clusters(self, clusters):
        return self.lim_to(self.clusters, clusters)

    @listify
    def limspaths2(self, spaths):
        if not is_iterable(spaths):
            spaths = [spaths]
        refs = set(self.refs)
        for sp in spaths:
            if sp.id in refs:
                yield sp
