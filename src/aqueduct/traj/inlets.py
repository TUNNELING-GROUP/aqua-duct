from collections import namedtuple


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
        self.clusters = None

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

    @property
    def coords(self):
        return [inlet.coords.tolist() for inlet in self.inlets_list]

    @property
    def types(self):
        return [inlet.type for inlet in self.inlets_list]

    @property
    def refs(self):
        return [inlet.reference for inlet in self.inlets_list]

