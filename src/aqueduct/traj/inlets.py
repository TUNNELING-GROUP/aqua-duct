from collections import namedtuple
import numpy as np
from aqueduct.utils.helpers import is_iterable

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


Inlet = namedtuple('Inlet', 'coords type reference')


class Inlets(object):
    # class for list of inlets
    def __init__(self, spaths, onlytype=InletTypeCodes.all_surface):

        self.onlytype = onlytype
        self.inlets_list = []
        self.clusters = []
        self.number_of_clustered_inlets = None

        for spath in spaths:
            self.extend_inlets(spath)

    def extend_inlets(self, spath, onlytype=None):

        if onlytype is None:
            onlytype = self.onlytype

        for inlet in spath.get_inlets():
            if onlytype is not None:
                if inlet.type not in onlytype:
                    continue
            self.inlets_list.append(inlet)

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

    def perform_clustering(self,method):
        self.add_cluster_annotations(method(np.array(self.coords)))
        self.number_of_clustered_inlets = len(self.clusters)

    @property
    def clusters_list(self):
        return sorted(list(set(self.clusters)))

    def lim_to(self, what, towhat):
        if not is_iterable(towhat):
            towhat = [towhat]
        new_inlets = self.__class__([], onlytype=self.onlytype)
        new_inlets.number_of_clustered_inlets = self.number_of_clustered_inlets

        for inlet, cluster,w in zip(self.inlets_list, self.clusters,what):
            if w in towhat:
                new_inlets.inlets_list.append(inlet)
                new_inlets.clusters.append(cluster)

        return new_inlets

    def lim2spaths(self, spaths):
        if not is_iterable(spaths):
            spaths = [spaths]
        return self.lim_to(self.refs,[sp.id for sp in spaths])

    def lim2types(self, types):
        return self.lim_to(self.types,types)

    def lim2clusters(self, clusters):
        return self.lim_to(self.clusters, clusters)




