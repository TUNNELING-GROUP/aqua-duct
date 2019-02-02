import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np


def polar2point(angle, r):
    rad = np.deg2rad
    return np.array([r * np.sin(rad(angle)), r * np.cos(rad(angle))])


def generate_arc(r, sa, ea, max_angle=45, reversed_=False):
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
        inner_arc = generate_arc(r - 1, sa, ea, reversed_=True)
        vertices.extend(inner_arc)

        codes.append(mpath.Path.LINETO)
        codes.extend([mpath.Path.CURVE3] * (len(inner_arc) - 1))

        # Closing poly
        vertices.append((0, 0))
        codes.append(mpath.Path.CLOSEPOLY)

        path = mpath.Path(vertices, codes)

        super(Node, self).__init__(path, facecolor=color, linewidth=0)

    def reserve_arc(self, angle):
        self._link_angle += angle

    def get_arc_offset(self):
        return self._link_angle


class Link(mpatches.PathPatch):
    def __init__(self, r, sa0, sa1, ea0, ea1, color):
        """
        Draw connection between two circle arcs.

        :param r: Radius.
        :param sa0: Source start angle.
        :param sa1: Source end angle.
        :param ea0: Destination start angle.
        :param ea1: Destination end angle.
        :param color: Color.
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
    def __init__(self, r, sa, ea, color):
        vertices = []
        codes = []

        arc = generate_arc(r, sa, ea)
        vertices.extend(arc)

        codes.append(mpath.Path.MOVETO)
        codes.extend([mpath.Path.CURVE3] * (len(arc) - 1))

        vertices.append(polar2point(sa + (ea - sa) / 2, r - 2))
        codes.append(mpath.Path.LINETO)

        vertices.append((0, 0))
        codes.append(mpath.Path.CLOSEPOLY)

        path = mpath.Path(vertices, codes)
        super(Arrow, self).__init__(path, facecolor=color, linewidth=0)


class Chord(object):
    def __init__(self, ax, r, nodes_sizes, links, labels, colors=None):
        self.nodes = []

        if colors:
            self.colors = colors
        else:
            self.colors = ["#FF33DD", "#554466", "#33FFDD", "#33DDFF"]

        sizes_sum = sum(nodes_sizes)

        # Nodes
        sa = 0
        for i, size in enumerate(nodes_sizes):
            ea = sa + size * 360. / sizes_sum

            angle = sa + (ea - sa) / 2
            pos = polar2point(angle, r + 0.5)
            ax.text(pos[0], pos[1],
                    labels[i],
                    verticalalignment="center",
                    horizontalalignment="center",
                    fontsize=8,
                    rotation=360 - angle)

            node = Node(r, sa, ea, self.colors[i])
            ax.add_patch(node)

            self.nodes.append(node)

            sa = ea

        # Links
        for link in links:
            source_node = self.nodes[link["source"]]
            dest_node = self.nodes[link["dest"]]
            value = link["value"]
            sa0 = source_node.get_arc_offset()
            sa1 = sa0 + value * 360. / sizes_sum

            source_node.reserve_arc(value * 360. / sizes_sum)

            ea0 = dest_node.get_arc_offset()
            ea1 = ea0 + value * 360. / sizes_sum

            dest_node.reserve_arc(value * 360. / sizes_sum)

            l = Link(r - 1, sa0, sa1, ea0, ea1, source_node.color)
            a = Arrow(r - 1, sa0, sa1, source_node.color)

            # # Arrow text
            # ax.text(pos[0], pos[1],
            #         labels[i],
            #         verticalalignment="center",
            #         horizontalalignment="center",
            #         fontsize=8,
            #         rotation=360-angle)

            ax.add_patch(l)
            ax.add_patch(a)


if __name__ == "__main__":
    fig, ax = plt.subplots(subplot_kw={"aspect": 1})
    ax.set_axis_off()
    ax.set_xlim(-11, 11)
    ax.set_ylim(-11, 11)

    labels = ["Data #1", "Data #2", "Data #3", "Data #4"]

    sizes = [120, 120, 120, 70]

    links = [dict(source=0, dest=1, value=10),
             dict(source=2, dest=2, value=30),
             dict(source=1, dest=2, value=10),
             dict(source=2, dest=3, value=20),
             dict(source=2, dest=0, value=10)]

    Chord(ax, 10, sizes, links, labels)

    plt.show()
