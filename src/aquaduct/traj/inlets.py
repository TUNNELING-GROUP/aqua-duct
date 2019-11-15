# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk,
#                          Magdalena Ługowska, Sandra Gołdowska <info@aquaduct.pl>
# Copyright (C) 2019  Tomasz Magdziarz <info@aquaduct.pl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging

logger = logging.getLogger(__name__)
import numpy as np
from collections import namedtuple, OrderedDict
from scipy.spatial.distance import pdist, squareform
import copy
from itertools import zip_longest
from functools import total_ordering

from aquaduct.utils.helpers import is_iterable, listify, lind
from aquaduct.utils import clui
from aquaduct.utils.maths import make_default_array
from aquaduct.geom.convexhull import SciPyConvexHull


class ProtoInletTypeCodes(object):
    surface = 'surface'
    internal = 'internal'

    incoming = 'inin'
    outgoing = 'inout'


class InletTypeCodes(ProtoInletTypeCodes):
    # TODO: Rework this and make it coherent with barber create_spheres.
    all_surface = tuple(((ProtoInletTypeCodes.surface, itype) for itype in
                         (ProtoInletTypeCodes.incoming, ProtoInletTypeCodes.outgoing)))
    all_internal = tuple(((ProtoInletTypeCodes.internal, itype) for itype in
                          (ProtoInletTypeCodes.incoming, ProtoInletTypeCodes.outgoing)))
    all_incoming = tuple(((itype, ProtoInletTypeCodes.incoming) for itype in
                          (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.internal)))
    all_outgoing = tuple(((itype, ProtoInletTypeCodes.outgoing) for itype in
                          (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.internal)))

    surface_incoming = (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.incoming)
    internal_incoming = (ProtoInletTypeCodes.internal, ProtoInletTypeCodes.incoming)
    internal_outgoing = (ProtoInletTypeCodes.internal, ProtoInletTypeCodes.outgoing)
    surface_outgoing = (ProtoInletTypeCodes.surface, ProtoInletTypeCodes.outgoing)


# clusers can be:
# None  - no such cluster
# nr - 0 means outliers
# because of passing paths one more type could be considered ???
# this might be also done in a different way


@total_ordering
class InletClusterGenericType(object):
    def __init__(self, inp, out):
        super(InletClusterGenericType, self).__init__()
        self.clusters = [inp, out]

    @property
    def input(self):
        return self.clusters[0]

    @property
    def output(self):
        return self.clusters[1]
        # return self.clusters[-1]

    @staticmethod
    def cluster2str(cl):
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

        def sort_for_cmp(a, b):

            # this function compares two elements and returns 0 for a == b, 1 for a > b , -1 for a < b

            if a == b:
                return 0
            elif a is None:
                if b > 0:
                    return 1
                elif b == 0:
                    return -1
            elif b is None:
                if a > 0:
                    return -1
                elif a == 0:
                    return 1
            elif a == 0:
                if b != 0:
                    return 1
            elif b == 0:
                if a != 0:
                    return -1
            elif a > b:
                return 1
            elif a < b:
                return -1

                # firstly only inputs are sorted, in case of self.input == other.input the outputs are sorted

        value = sort_for_cmp(self.input, other.input)
        if value == 0:
            return sort_for_cmp(self.output, other.output)
        else:
            return value

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __hash__(self):
        return hash(str(self))


class InletClusterExtendedType(InletClusterGenericType):
    def __init__(self, surfin, interin, interout, surfout):
        super(InletClusterExtendedType, self).__init__(surfin, surfout)
        self.clusters.extend([interin, interout])

    @property
    def generic(self):
        return InletClusterGenericType(*self.clusters[:2])


# InletOLD = namedtuple('Inlet', 'coords type reference frame') # for backward compatibility

class Inlet(object):
    def __init__(self, coords=None, inlet_type=None, reference=None, frame=None):
        # no none is allowed
        assert None not in [inlet_type, reference, frame], "Wrong Inlet init."
        self.coords = make_default_array(coords)
        self.type = inlet_type
        self.reference = reference
        self.frame = frame


class Inlets(object):
    # class for list of inlets
    def __init__(self, spaths, center_of_system=None, onlytype=InletTypeCodes.all_surface, passing=False, pbar=None):

        self.center_of_system = center_of_system
        self.onlytype = onlytype
        self.inlets_list = []
        self.inlets_ids = []
        self.clusters = []
        self.number_of_clustered_inlets = None
        self.spheres = []
        self.passing = passing

        self.tree = clui.SimpleTree()

        for spath in spaths:
            self.extend_inlets(spath)
            if pbar: pbar.next()

    def __iter__(self):
        return (inlet for inlet in self.inlets_list)

    def add_leaf_wrapper(self, name=None, message=None, toleaf=None):
        if name == 0:
            self.tree.add_leaf(name='(0)', message=message, toleaf=toleaf)
        else:
            self.tree.add_leaf(name=name, message=message, toleaf=toleaf)

    def resize_leaf_0(self):
        if 0 in self.clusters_list:
            if 0 not in self.tree.leafs_names:
                self.tree.add_leaf(name=0, message='size: %d' % self.clusters.count(0))
            else:
                self.tree.add_message(toleaf=0, message='size: %d' % self.clusters.count(0), replace=True)

    def add_message_wrapper(self, message=None, toleaf=None):
        self.tree.add_message(message=message, toleaf=toleaf)

    def extend_inlets(self, spath, onlytype=None):
        if not self.passing and spath.is_passing():
            return

        if onlytype is None:
            onlytype = self.onlytype

        added_list = []
        nr = len(self.inlets_list)
        for inlet in spath.get_inlets():
            if onlytype is not None:
                if inlet.type not in onlytype:
                    continue
            # is ther already clustering info?
            if len(self.clusters):
                self.clusters.append(0)
            self.inlets_list.append(inlet)
            self.inlets_ids.append(nr)
            added_list.append(nr)
            nr += 1
        return added_list

    def add_cluster_annotations(self, clusters):
        # this replaces clusters!
        assert len(clusters) == len(self.inlets_list)
        self.clusters = clusters
        self.tree = clui.SimpleTree()

    def add_outliers_annotations(self, new_clusters):
        assert len(new_clusters) == len(self.inlets_list)
        # this is meant for situation when some points from clusters are changed to outliers
        # lets loop over current clusters list
        for cluster, cluster_size in zip(self.clusters_list, self.clusters_size):
            if cluster == 0: continue
            # check if cluster was changed!
            if cluster not in new_clusters:
                # it was completly removed!
                self.add_message_wrapper(message='outliers detection', toleaf=cluster)
                self.add_leaf_wrapper(name=0, toleaf=cluster, message='size: %d' % cluster_size)
            elif self.clusters.count(cluster) != new_clusters.count(cluster):
                # some points were shifted to outliers
                self.add_message_wrapper(message='outliers detection', toleaf=cluster)
                self.add_leaf_wrapper(name=0, toleaf=cluster,
                                      message='size: %d' % (cluster_size - new_clusters.count(cluster)))
                # we need new cluster!
                new_cluster = max(new_clusters) + 1
                for nr, dummy in enumerate(new_clusters):
                    if new_clusters[nr] == cluster:
                        new_clusters[nr] = new_cluster
                self.add_leaf_wrapper(name=new_cluster, toleaf=cluster,
                                      message='size: %d' % (new_clusters.count(new_cluster)))
        self.clusters = new_clusters
        self.resize_leaf_0()

    def add_spheres(self, spheres):
        assert len(spheres) == len(self.inlets_list)
        self.spheres = spheres

    def get_inlets_references(self):
        return [inl.reference for inl in self.inlets_list]

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

    @property
    def refs_names(self):
        return [inlet.reference.name for inlet in self.inlets_list]

    def call_clustering_method(self, method, data, spheres=None):
        # this method runs clustering method against provided data
        # if center_of_system was set then use distance matrix...
        if self.center_of_system is not None:
            return method(np.array(data) - self.center_of_system, spheres=spheres)
        return method(np.array(data), spheres=spheres)

    def get_flat_tree(self, message=None):
        st = clui.SimpleTree(name='all', message='size: %d' % self.size)
        st.add_message(message=message)
        [st.add_leaf(name=leaf, message='size: %d' % csize) for leaf, csize in
         zip(self.clusters_list, self.clusters_size)]
        return st

    def perform_clustering(self, method):
        # this do clean clustering, all previous clusters are discarded
        # 0 means outliers
        self.add_cluster_annotations(self.call_clustering_method(method, self.coords, spheres=self.spheres))
        self.number_of_clustered_inlets = len(self.clusters)
        clui.message("New clusters created: %s" % (' '.join(map(str, sorted(set(self.clusters))))))
        # renumber clusters
        # self.renumber_clusters()
        # return clusters as simple tree
        self.tree = self.get_flat_tree(message=str(method))

    def perform_reclustering(self, method, skip_outliers=False, skip_size=None):
        # this do reclustering of all clusters, if no cluster exists perform_clustering is called
        if len(self.clusters) == 0:
            return self.perform_clustering(method)
        for cluster in self.clusters_list:
            if skip_outliers and cluster == 0:
                clui.message('Skipping outliers.')
                continue
            # check cluster size and skip if does not fit to skip_thershold function
            cluster_size = float(self.clusters.count(cluster)) / self.size
            if skip_size is not None:
                if skip_size(cluster_size):
                    clui.message('Cluster %d of size %0.3f skipped.' % (cluster, cluster_size))
                    continue
            clui.message('Cluster %d of size %0.3f submitted to reclustering.' % (cluster, cluster_size))
            self.recluster_cluster(method, cluster)
        # number of cluster
        self.number_of_clustered_inlets = len(self.clusters)
        # renumber clusters
        # self.renumber_clusters()
        # return clusters as simple tree

    # CLUSTER
    def recluster_cluster(self, method, cluster):
        if cluster in self.clusters_list:
            logger.debug('Reclustering %d cluster: initial number of clusters %d.' % (cluster, len(self.clusters_list)))
            reclust = self.call_clustering_method(method, self.lim2clusters(cluster).coords,
                                                  spheres=self.lim2clusters(cluster).spheres)
            if len(set(reclust)) <= 1:
                clui.message('No new clusters found.')
            else:
                # how many clusters? what aboout outliers?
                n_reclust = len(set(reclust))
                out_reclust = 0 in reclust
                clui.message('Cluster %d was split into %d clusters.' % (cluster, n_reclust))
                if out_reclust:
                    clui.message('Outliers were detected and will have annotation of the old cluster %d.' % cluster)
                else:
                    clui.message('No outliers were detected, the old cluster %d will be removed.' % cluster)
                # change numbers of reclust
                max_cluster = max(self.clusters_list)
                for nr, r in enumerate(reclust):
                    if r != 0:
                        reclust[nr] = r + max_cluster
                if cluster != 0:
                    self.add_message_wrapper(message=str(method), toleaf=cluster)
                [self.add_leaf_wrapper(name=leaf, toleaf=cluster, message=['size: %d' % reclust.count(leaf)]) for leaf
                 in sorted(list(set(reclust)))]
                if cluster == 0:
                    self.add_message_wrapper(message=['[RE]', str(method)], toleaf=cluster)
                if out_reclust:
                    clui.message('The old cluster %d will be split into new clusters: %s' % (
                        cluster, (' '.join(map(str, sorted(set(reclust))[1:])))))
                else:
                    clui.message('The old cluster %d will be split into new clusters: %s' % (
                        cluster, (' '.join(map(str, sorted(set(reclust)))))))
                # add new clusters
                nrr = 0
                for nr, c in enumerate(self.clusters):
                    if c == cluster:
                        self.clusters[nr] = reclust[nrr]
                        nrr += 1
            logger.debug('Reclustering %d cluster: final number of clusters %d.' % (cluster, len(self.clusters_list)))
        # number of cluster
        self.number_of_clustered_inlets = len(self.clusters)
        if cluster != 0:
            self.resize_leaf_0()

    def recluster_outliers(self, method):
        self.recluster_cluster(method, 0)
        # renumber clusters
        # self.renumber_clusters()

    def remove_inlets(self, cluster):
        # removes inlets completly - clusters and inlets
        n = self.clusters.count(cluster)
        to_remove = []
        for nr, c in enumerate(self.clusters):
            if c == cluster:
                to_remove.append(nr)
        for tr in to_remove[::-1]:
            self.clusters.pop(tr)
            self.inlets_list.pop(tr)
            self.inlets_ids.pop(tr)
        self.add_message_wrapper(message='cluster %d removed' % cluster, toleaf=cluster)

    def small_clusters_to_outliers(self, maxsize):
        new_out = 0
        for c in self.clusters_list:
            if c == 0:
                continue
            if self.clusters.count(c) <= maxsize:
                self.add_leaf_wrapper(name=0, toleaf=c, message='size: %d' % self.clusters.count(c))
                self.add_message_wrapper(message='|%d| to outliers' % maxsize, toleaf=c)
                for nr, cc in enumerate(self.clusters):
                    if cc == c:
                        self.clusters[nr] = 0
                        new_out += 1
                        # renumber clusters
                        # self.renumber_clusters()
        if new_out:
            if 0 not in self.tree.leafs_names:
                self.tree.add_leaf(name=0)
            self.add_message_wrapper(message=['|%d| to outliers' % maxsize, 'new size: %d' % self.clusters.count(0)],
                                     toleaf=0)
            # self.resize_leaf_0()

    def renumber_clusters(self):
        old_clusters = list(self.clusters)
        if 0 in self.clusters_list:
            new_numbers = list(range(len(self.clusters_list)))
        else:
            new_numbers = list(range(1, len(self.clusters_list) + 1))
        old_numbers = self.clusters_list
        if old_numbers != new_numbers:
            for nr, c in enumerate(self.clusters):
                self.clusters[nr] = new_numbers[old_numbers.index(c)]
        self.sort_clusters()
        # make old new dict
        ond = []
        for o, n in zip(old_clusters, self.clusters):
            if o not in ond:
                ond.append(o)
                self.add_message_wrapper(toleaf=o, message='new nr: %d' % n)

    def sort_clusters(self):
        old_numbers = self.clusters_list
        old_sizes = self.clusters_size
        # now sort according to sizes but put 0, if present, at the begining
        if 0 in old_numbers:
            zero_size = old_sizes.pop(old_numbers.index(0))
            zero_number = old_numbers.pop(old_numbers.index(0))  # which is zero!
            new_numbers = [zero_number] + [old_numbers[i] for i in np.argsort(old_sizes).tolist()[::-1]]
            new_sizes = [zero_size] + [old_sizes[i] for i in np.argsort(old_sizes).tolist()[::-1]]
            old_numbers = self.clusters_list
            # old_sizes = self.clusters_size
        else:
            new_numbers = [old_numbers[i] for i in np.argsort(old_sizes).tolist()[::-1]]
            new_sizes = [old_sizes[i] for i in np.argsort(old_sizes).tolist()[::-1]]
        # trans_dict = {n: o for o, n in zip(old_numbers, new_numbers)}
        trans_dict = dict((n, o) for o, n in
                          zip(old_numbers, new_numbers))  # more universal as dict comprehension may not work in <2.7
        new_clusters = []
        for c in self.clusters:
            new_clusters.append(trans_dict[c])
        self.clusters = new_clusters
        assert self.clusters_size == new_sizes

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
        return list(map(self.clusters.count, self.clusters_list))

    def join_clusters(self, clusters2join):
        new_cluster = max(self.clusters_list) + 1
        for nr, c in enumerate(self.clusters):
            if c in clusters2join:
                self.clusters[nr] = new_cluster
        for c in clusters2join:
            self.add_message_wrapper(toleaf=c, message='joined into: %d' % new_cluster)
        nc_size = self.clusters.count(new_cluster)
        self.add_leaf_wrapper(name=new_cluster, message='size: %d' % nc_size, toleaf=None)
        self.add_message_wrapper(message='made of %s' % ('+'.join(map(str, clusters2join))), toleaf=new_cluster)

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
            yield self.spath2ctype(sp)

    def spath2ctype(self, sp):
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
        return InletClusterExtendedType(surfin, interin, interout, surfout)

    def lim_to(self, what, towhat):
        if not is_iterable(towhat):
            towhat = [towhat]
        new_inlets = self.__class__([], onlytype=self.onlytype,
                                    center_of_system=self.center_of_system,
                                    passing=self.passing,
                                    pbar=None)  # do not propagate pabr
        new_inlets.number_of_clustered_inlets = self.number_of_clustered_inlets

        for inlet, ids, cluster, sphere, w in zip_longest(self.inlets_list, self.inlets_ids, self.clusters,
                                                           self.spheres, what):
            if w in towhat:
                new_inlets.inlets_list.append(inlet)
                new_inlets.inlets_ids.append(ids)
                new_inlets.clusters.append(cluster)
                if len(self.spheres) and (sphere or sphere == 0):
                    new_inlets.spheres.append(sphere)

        return new_inlets

    def lim2spaths(self, spaths):
        if not is_iterable(spaths):
            spaths = [spaths]
        return self.lim_to(self.refs, [sp.id for sp in spaths])

    def lim2rnames(self, rnames):
        if not is_iterable(rnames):
            rnames = [rnames]
        return self.lim_to(self.refs_names, rnames)

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

    def get_chull(self):
        return SciPyConvexHull(self.coords)
