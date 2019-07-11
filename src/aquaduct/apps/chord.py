#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Michał Banas
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

import colorsys
import itertools

import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np


def color_gen():
    colors = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
              "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324",
              "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"]
    for color in itertools.cycle(colors):
        yield color


def hex2rgb(color):
    """
    Convert HEX color to RGB format.
    :param color: String hex color.
    :return: Tuple with RGB values.
    """
    # https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python
    return tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))


def polar2point(angle, r):
    """
    Transform polar coordinates to cartesian coordinates.
    :param angle: Angle.
    :param r: Radius.
    """
    rad = np.deg2rad
    return np.array([r * np.sin(rad(angle)), r * np.cos(rad(angle))])


def generate_arc(r, sa, ea, max_angle=5, reversed_=False):
    """
    Generate arc vertrices with control points for quadratic Bezier curve.

    :param r: Radius
    :param sa: Start angle.
    :param ea: End angle
    :param max_angle: Max. angle for which control point will be calculated. If > 90 curves will be significantly distorted.
    :param reversed_: If True vertices will start from ea.

    If `ea - sa > max_angle`, then `ea - sa` will be divided by `max_angle` into parts and control point will be calculated for each part.
    """
    vertices = []

    angles = []
    n = ea - sa
    while True:
        if n - max_angle >= 0:
            angles.append(max_angle)

            if n - max_angle == 0:
                break

            n -= max_angle
        else:
            angles.append(n % max_angle)
            break

    if reversed_:
        sa = -ea
        angles = reversed(angles)

    vertices.append(polar2point(np.abs(sa), r))

    for angle in angles:
        cp_angle = np.abs(sa + angle / 2.)
        # https://stackoverflow.com/questions/1734745/how-to-create-circle-with-b%C3%A9zier-curves
        cp_length = r + r * np.pi * 0.1 * pow(angle / 90., 2)

        vertices.append(polar2point(cp_angle, cp_length))
        vertices.append(polar2point(np.abs(sa + angle), r))

        sa += angle

    return vertices


class Node(mpatches.PathPatch):
    def __init__(self, r, sa, ea, color):
        """
        Represent data on Chord circle.

        :param r: Radius.
        :param sa: Start angle.
        :param ea: End angle.
        :param color: Node color in HEX or matplotlib tuple format.
        """
        self.sa = sa
        self.ea = ea
        self.color = color

        # Remembers when last link was added
        self._link_angle = sa

        vertices = []
        codes = []

        # # Outer arc
        outer_arc = generate_arc(r, sa, ea)
        vertices.extend(outer_arc)

        codes.append(mpath.Path.MOVETO)
        codes.extend([mpath.Path.CURVE3] * (len(outer_arc) - 1))

        # Inner arc
        inner_arc = generate_arc(0.8 * r, sa, ea, reversed_=True)
        vertices.extend(inner_arc)

        codes.append(mpath.Path.LINETO)
        codes.extend([mpath.Path.CURVE3] * (len(inner_arc) - 1))

        # Closing poly
        vertices.append((0, 0))
        codes.append(mpath.Path.CLOSEPOLY)

        path = mpath.Path(vertices, codes)

        super(Node, self).__init__(path, facecolor=color, linewidth=0)

    def reserve_arc(self, angle):
        # Only if scaling is an optional
        # if self._link_angle + angle - self.sa > self.ea - self.sa:
        #     raise RuntimeError("Too much data for node. Use scaling option.")

        self._link_angle += angle

    def get_arc_offset(self):
        return self._link_angle


class Link(mpatches.PathPatch):
    def __init__(self, r, sa0, sa1, ea0, ea1, color):
        """
        Represent connection between two nodes.

        :param r: Radius.
        :param sa0: Source start angle.
        :param sa1: Source end angle.
        :param ea0: Destination start angle.
        :param ea1: Destination end angle.
        :param color: Link color in HEX or matplotlib tuple format.
        """
        vertices = []
        codes = []

        # Source arc
        source_arc = generate_arc(r, sa0, sa1)
        vertices.extend(source_arc)

        codes.append(mpath.Path.MOVETO)
        codes.extend([mpath.Path.CURVE3] * (len(source_arc) - 1))

        # Connection line
        vertices.append(polar2point(sa1, r))
        codes.append(mpath.Path.LINETO)

        vertices.append((0, 0))
        codes.append(mpath.Path.CURVE3)

        vertices.append(polar2point(ea0, r))
        codes.append(mpath.Path.CURVE3)

        # Dest arc
        dest_arc = generate_arc(r, ea0, ea1)
        vertices.extend(dest_arc)

        codes.extend([mpath.Path.CURVE3] * len(dest_arc))

        # Connection line
        vertices.append(polar2point(ea1, r))
        codes.append(mpath.Path.LINETO)

        vertices.append((0, 0))
        codes.append(mpath.Path.CURVE3)

        vertices.append(polar2point(sa0, r))
        codes.append(mpath.Path.CURVE3)

        path = mpath.Path(vertices, codes)
        super(Link, self).__init__(path, facecolor=color, linewidth=0, alpha=0.6)


class Arrow(mpatches.PathPatch):
    def __init__(self, r, sa, ea, color, max_angle=45):
        """
        Arrow patch for links.

        :param r: Radius.
        :param sa: Start angle.
        :param ea: End angle.
        :param color: Arrow color in HEX or matplotlib tuple format.
        :param max_angle:
        """
        vertices = []
        codes = []

        c = sa + (ea - sa) / 2

        if ea - sa < max_angle:
            max_angle = ea - sa

        arc = generate_arc(r, c - max_angle / 2, c + max_angle / 2)
        vertices.extend(arc)

        codes.append(mpath.Path.MOVETO)
        codes.extend([mpath.Path.CURVE3] * (len(arc) - 1))

        vertices.append(polar2point(sa + (ea - sa) / 2, 0.9 * r - 2))
        codes.append(mpath.Path.LINETO)

        vertices.append((0, 0))
        codes.append(mpath.Path.CLOSEPOLY)

        path = mpath.Path(vertices, codes)
        super(Arrow, self).__init__(path, facecolor=color, linewidth=0)


class Chord(object):
    def __init__(self, ax, r, nodes_sizes, links, labels, colors=[]):
        """
        
        :param ax:
        :param r:
        :param nodes_sizes:
        :param links:
        :param labels:
        :param colors:
        """
        self.nodes = []

        if colors:
            self.colors = colors
        else:
            cg = color_gen()
            self.colors = [next(cg) for _ in range(0, len(nodes_sizes))]

        sizes_sum = sum(nodes_sizes)

        # Nodes
        legend_node_index = []  # Keeps indexes of nodes which will be in legend
        sa = 0
        for i, size in enumerate(nodes_sizes):
            ea = sa + size * 360. / sizes_sum

            if ea - sa >= 10:
                angle = sa + (ea - sa) / 2
                pos = polar2point(angle, 1.05 * r)
                ax.text(pos[0], pos[1],
                        "{} ({:.2f}%)".format(labels[i], 100.0 * size / sizes_sum),
                        verticalalignment="center",
                        horizontalalignment="center",
                        fontsize=10,
                        rotation=360 - angle)
            else:
                legend_node_index.append(i)

            node = Node(r, sa, ea, self.colors[i])
            ax.add_patch(node)

            self.nodes.append(node)

            sa = ea

        # Create legend
        ax.legend([self.nodes[i] for i in legend_node_index],
                  ["{} ({:.2f}%)".format(labels[i], 100.0 * nodes_sizes[i] / sizes_sum) for i in legend_node_index],
                  bbox_to_anchor=(1.02, 1), loc="upper left")

        # Scales
        links_sum = [0.] * len(nodes_sizes)
        for link in links:
            links_sum[link["source"]] += link["value"]
            links_sum[link["dest"]] += link["value"]

        scales = []
        for i, (size, link_sum) in enumerate(zip(nodes_sizes, links_sum)):
            if size < link_sum:
                scales.append(0.9999 * size / link_sum)  # Scale is reduced to 99% due to float precision
            else:
                scales.append(1.)

        # Links
        for link in links:
            source_node = self.nodes[link["source"]]
            dest_node = self.nodes[link["dest"]]
            value = link["value"]
            sa0 = source_node.get_arc_offset()
            slink_arc = scales[link["source"]] * value * 360. / sizes_sum
            sa1 = sa0 + slink_arc

            source_node.reserve_arc(slink_arc)

            ea0 = dest_node.get_arc_offset()
            elink_arc = scales[link["dest"]] * value * 360. / sizes_sum
            ea1 = ea0 + elink_arc

            dest_node.reserve_arc(elink_arc)

            l = Link(0.8 * r, sa0, sa1, ea0, ea1, source_node.color)
            a = Arrow(0.81 * r, sa0, sa1, source_node.color)

            ax.add_patch(l)
            ax.add_patch(a)

            if sa1 - sa0 > 6:
                # Complimentary color
                rgb_color = hex2rgb(source_node.color.lstrip("#")) if isinstance(source_node.color,
                                                                                 str) else source_node.color
                hsl_color = list(colorsys.rgb_to_hls(*rgb_color))

                if not hsl_color[0]:
                    if 138 >= hsl_color[1] >= 118:  # For gray color automatically set white
                        hsl_color[1] = 255
                    else:
                        hsl_color[1] = (hsl_color[1] - 255) * -1.0
                else:
                    hsl_color[0] += 0.5

                complimentary_color = [c / 255 for c in colorsys.hls_to_rgb(*hsl_color)]

                # Arrow text
                angle = sa0 + (sa1 - sa0) / 2
                pos = polar2point(sa0 + (sa1 - sa0) / 2, 0.86 * r)
                ax.text(pos[0], pos[1],
                        link["value"],
                        verticalalignment="center",
                        horizontalalignment="center",
                        fontsize=9,
                        rotation=180 - angle,
                        color="white")


if __name__ == "__main__":
    fig, ax = plt.subplots(subplot_kw={"aspect": 1})
    ax.set_axis_off()
    ax.set_xlim(-110, 110)
    ax.set_ylim(-110, 110)

    labels = ["Data #1", "Data #2", "Data #3", "Data #4"]

    sizes = [5000, 6000, 7000]

    links = [dict(source=1, dest=0, value=2000),
             dict(source=1, dest=0, value=2000),
             dict(source=2, dest=0, value=4600),
             dict(source=0, dest=0, value=1000)]

    Chord(ax, 100, sizes, links, labels)

    fig.savefig("chord.png", format="png", dpi=2 ** 7)

    plt.show()
