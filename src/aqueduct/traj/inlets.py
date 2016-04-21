from collections import namedtuple
import numpy as np
from aqueduct.utils.helpers import is_iterable, listify


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

    def __cmp__(self, other):
        if other is None:
            return 1
        if not isinstance(other,self.__class__):
            return 1

        result = 0
        base = max(max(self), max(other), len(self), len(other)) + 2

        def make_val(what):
            val = 0
            for nr, e in enumerate(tuple(what)[::-1]):
                if e is None:
                    val += (base ** nr) * base
                else:
                    val += (base ** nr) * e
            return val

        return make_val(self) - make_val(other)


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
        # 0 means outliers
        self.add_cluster_annotations(method(np.array(self.coords)))
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

    @property
    def clusters_list(self):
        return sorted(list(set(self.clusters)))

    @property
    @listify
    def clusters_centers(self):
        for c in self.clusters_list:
            yield np.mean(self.lim2clusters(c).coords,0)


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
    def limspaths2(self,spaths):
        if not is_iterable(spaths):
            spaths = [spaths]
        refs = set(self.refs)
        for sp in spaths:
            if sp.id in refs:
                yield sp