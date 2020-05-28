#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Michał Banas <-- an absolute tosser
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

from __future__ import absolute_import
from __future__ import print_function
import tkinter as tk
import base64
import csv
import re
import tkinter.messagebox
import tkinter.ttk
from io import BytesIO
from collections import defaultdict
from tkinter.filedialog import askopenfile

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes

import aquaduct
import aquaduct.apps.valveconfig.utils as utils
from aquaduct.apps.chord import Chord, color_gen
from aquaduct.apps.valveconfig import get_img


def is_float(value):
    try:
        float(value)
        return True
    except:
        return False


class DataException(Exception):
    pass


class FileDataProcessor(object):
    def __init__(self, filename):
        self.file = open(filename, "r")
        self.file_lines = [line.strip() for line in self.file.readlines()]

        self.table_end_pattern = "-+"

    def get_column_values(self, table_name, column_name):
        data = [row[column_name] for row in self._find_table(table_name)]

        return data

    def get_column_names(self, table_name):
        data = self._find_table(table_name)

        if len(data):
            return list(data[0].keys())
        else:
            raise DataException("Table \"{}\" is empty.".format(table_name))

    def _find_table(self, table_name):
        for i, line in enumerate(self.file_lines):
            if line.strip() == table_name:
                return self._parse_column(i)

        raise DataException("Table \"{}\" does not exist.".format(table_name))

    def _parse_column(self, start_line):
        data = []

        column_names = self.file_lines[start_line + 2].split()

        for line in self.file_lines[start_line + 4:]:
            if re.match(self.table_end_pattern, line):
                break

            values = []
            for value in line.split():
                if value.isdigit():
                    value = int(value)
                elif is_float(value):
                    value = float(value)

                values.append(value)

            dict_values = {}
            for i in range(0, len(values)):
                dict_values.update({column_names[i]: values[i]})

            data.append(dict_values)

        return data

    def __del__(self):
        self.file.close()


class CSVDataProcessor(object):
    def __init__(self, filename):
        self.file = open(filename, "r")
        self.csv_reader = csv.DictReader(self.file)

        self.column_names = list(next(self.csv_reader).keys())
        self.seek_file()

    def get_column_values(self, column_name):
        values = [float(row[column_name]) if is_float(row[column_name]) else int(row[column_name]) for row in
                  self.csv_reader]
        self.seek_file()

        return values

    def seek_file(self):
        """ Set file position at beginning and skip first row with column names """
        self.file.seek(0)
        next(self.csv_reader)

    def __del__(self):
        self.file.close()


# 1
def cluster_inlets(file_processor, suffix=""):
    fig, ax = plt.subplots()

    cluster_no = file_processor.get_column_values("Clusters summary - inlets" + suffix, "Cluster")

    y_in = file_processor.get_column_values("Clusters summary - inlets" + suffix, "INCOMING")
    y_out = file_processor.get_column_values("Clusters summary - inlets" + suffix, "OUTGOING")

    # Added 1 to make place for previous bar
    rects_bottom = ax.bar(list(range(1, len(cluster_no) + 1)), y_in, label="INCOMING")
    rects_top = ax.bar(list(range(1, len(cluster_no) + 1)), y_out, bottom=y_in, label="OUTGOING")

    y_lim = ax.get_ylim()[1]

    for rect_b, rect_t in zip(rects_bottom, rects_top):
        width_b = rect_b.get_width()
        height_b = rect_b.get_height()
        if height_b / y_lim >= 0.02:
            ax.text(rect_b.get_x() + width_b / 2., height_b / 2., height_b,
                    fontsize=7,
                    color=(1, 1, 1),
                    horizontalalignment="center",
                    verticalalignment="center")

        width_t = rect_t.get_width()
        height_t = rect_t.get_height()
        if height_t / y_lim >= 0.02:
            ax.text(rect_t.get_x() + width_t / 2., height_b + height_t / 2., height_t,
                    fontsize=7,
                    color=(1, 1, 1),
                    horizontalalignment="center",
                    verticalalignment="center")

    ax.set_title("Clusters size" + suffix)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Size")

    ax.set_xlim((0, len(cluster_no) + 1))
    ax.set_xticks(list(range(1, len(cluster_no) + 1)))
    ax.set_xticklabels(cluster_no)
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    return fig


# 2
def relative_clusters_flows(csv_processor, clusters_names, labels, colors):
    fig, ax = plt.subplots()

    x = csv_processor.get_column_values("# frame")

    all_paths_per_frame = csv_processor.get_column_values("amol_apaths_aclusts_walk")

    for cluster_name, label, color in zip(clusters_names, labels, colors):
        cluster_paths_per_frame = csv_processor.get_column_values(cluster_name)

        y = [float(cluster_paths) / all_paths if all_paths else 0
             for cluster_paths, all_paths in zip(cluster_paths_per_frame, all_paths_per_frame)]

        ax.plot(x, y, label=label, color=color)

    ax.set_title("Relative cluster flows")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Relative flow")

    ax.set_xlim((0, len(x)))
    ax.set_ylim((0, 1.05))  # added 0.05 because when there is 1, its not visible

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    return fig


# 3
def ligands_time(file_processor, molecule=None):
    plot_settings = dict(align="edge", height=1.0)

    title_suffix = "" if not molecule else " of {}".format(molecule)

    path_nr = file_processor.get_column_values("List of separate paths and properties", "Nr")

    path_nr = [y - 1 for y in path_nr]

    beg_values = file_processor.get_column_values("List of separate paths and properties", "BeginF")
    inp_values = file_processor.get_column_values("List of separate paths and properties", "InpF")
    obj_values = file_processor.get_column_values("List of separate paths and properties", "ObjF")
    out_values = file_processor.get_column_values("List of separate paths and properties", "OutF")

    # Get max frame to set plot width
    max_frame = max(file_processor.get_column_values("List of separate paths and properties", "EndF"))

    if molecule:
        residues_names = file_processor.get_column_values("List of separate paths and properties", "RES")

        for i, res_name in reversed(list(enumerate(residues_names))):
            if res_name != molecule:
                del path_nr[i]
                del beg_values[i]
                del inp_values[i]
                del obj_values[i]
                del out_values[i]

    # Sort all basing on BeginF
    beg_values, inp_values, obj_values, out_values = list(zip(*sorted(zip(beg_values,
                                                                     inp_values,
                                                                     obj_values,
                                                                     out_values))))

    inps_left = beg_values
    inps_width = inp_values

    objs_left = [sum(values) for values in zip(beg_values, inp_values)]
    objs_width = obj_values

    outs_left = [sum(values) for values in zip(beg_values, inp_values, obj_values)]
    outs_width = out_values

    fig = plt.figure(figsize=(6, 6))
    fig_dpi = fig.get_dpi()

    h = [Size.Fixed(1.0), Size.Fixed(max_frame/fig_dpi)]
    v = [Size.Fixed(1.0), Size.Fixed(len(path_nr)/fig_dpi)]

    fig.set_size_inches(max_frame/fig_dpi + 2.0, len(path_nr)/fig_dpi + 2.0)

    divider = Divider(fig, (0.0, 0.0, 1., 1.), h, v, aspect=False)
    # the width and height of the rectangle is ignored.

    ax = Axes(fig, divider.get_position())
    ax.set_axes_locator(divider.new_locator(nx=1, ny=1))

    fig.add_axes(ax)

    ax.barh(path_nr, inps_width, left=inps_left, color=(1, 0, 0), label="Incoming", **plot_settings)
    ax.barh(path_nr, objs_width, left=objs_left, color=(0, 1, 0), label="Object", **plot_settings)
    ax.barh(path_nr, outs_width, left=outs_left, color=(0, 0, 1), label="Outgoing", **plot_settings)

    ax.set_title("Molecule Entry Time Distribution" + title_suffix)
    ax.set_xlabel("Frame")
    ax.set_ylabel("Sorted paths")

    ax.set_ylim((0, max(path_nr)))
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    return fig


# 4
def chord_diagram_sizes(file_processor, labels={}, colors={}):
    fig, ax = plt.subplots(subplot_kw={"aspect": 1})
    ax.set_title("Cluster size")
    ax.set_axis_off()
    ax.set_xlim(-110, 110)
    ax.set_ylim(-110, 110)

    clusters_ids = file_processor.get_column_values("Clusters summary - inlets", "Cluster")
    clusters_ids = [str(id_) for id_ in clusters_ids]
    clusters_sizes = file_processor.get_column_values("Clusters summary - inlets", "Size")

    # TODO:  Move it to generate() method
    if not labels:
        labels = clusters_ids
    else:
        tmp_labels = []
        for id_ in clusters_ids:
            tmp_labels.append(labels[id_])

        labels = tmp_labels

    if colors:
        tmp_colors = []
        for id_ in clusters_ids:
            tmp_colors.append(colors[id_])

        colors = tmp_colors

    Chord(ax, 100, clusters_sizes, [], labels, colors)

    return fig


# 4
def chord_diagram_flows(file_processor, labels={}, colors={}, threshold=0.):
    fig, ax = plt.subplots(subplot_kw={"aspect": 1})
    ax.set_title("Intramolecular flows")
    ax.set_axis_off()
    ax.set_xlim(-110, 110)
    ax.set_ylim(-110, 110)

    flows = file_processor.get_column_values("Separate paths clusters types summary - mean lengths of paths", "CType")
    flows = [{flow.split(":")[0]: flow.split(":")[1]} for flow in flows]

    flows_sizes = file_processor.get_column_values("Separate paths clusters types summary - mean lengths of paths",
                                                   "Size")

    sizes_d = defaultdict(int)
    for flow, size in zip(flows, flows_sizes):
        source, dest = next(iter(flow.items()))

        if source != "N":
            source = int(source)

        if dest != "N":
            dest = int(dest)

        sizes_d[source] += size
        sizes_d[dest] += size

    clusters_ids = []
    sizes = []
    #cannot sort ints and strings together
    #for k, v in sorted(iter(sizes_d.items()), key=lambda x: x[0]):
    for k, v in (iter(sizes_d.items())):
        clusters_ids.append(str(k))
        sizes.append(v)

    flows = [{"source": clusters_ids.index(next(iter(flow.keys()))), "dest": clusters_ids.index(next(iter(flow.values()))),
              "value": size} for flow, size in zip(flows, flows_sizes) if 1. * size / sum(flows_sizes) >= threshold]

    if not labels:
        labels = clusters_ids
    else:
        tmp_labels = []
        threshold = 0.
        for id_ in clusters_ids:
            tmp_labels.append(labels[id_])

        labels = tmp_labels

    if colors:
        tmp_colors = []
        for id_ in clusters_ids:
            tmp_colors.append(colors[id_])

        colors = tmp_colors

    Chord(ax, 100, sizes, flows, labels, colors)

    return fig


# 7
def volume_scope_area(csv_processor):
    fig, ax = plt.subplots()
    ax.set_title("Scope area")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Area \u00C5^2$")

    x = csv_processor.get_column_values("# frame")
    y = csv_processor.get_column_values("scope_area")

    ax.set_xlim((min(x), max(x)))
    ax.set_ylim((min(y), max(y)))

    ax.plot(x, y)

    return fig


# 7
def volume_scope_volume(csv_processor):
    fig, ax = plt.subplots()
    ax.set_title("Scope volume")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Volume \u00C5^3$")

    x = csv_processor.get_column_values("# frame")
    y = csv_processor.get_column_values("scope_volume")

    ax.plot(x, y)

    ax.set_xlim((min(x), max(x)))
    ax.set_ylim((min(y), max(y)))

    return fig


# 7
def volume_object_area(csv_processor):
    fig, ax = plt.subplots()
    ax.set_title("Object area")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Area \u00C5^2$")

    x = csv_processor.get_column_values("# frame")
    y = csv_processor.get_column_values("object_area")

    ax.plot(x, y)

    ax.set_xlim((min(x), max(x)))
    ax.set_ylim((min(y), max(y)))

    return fig


# 7
def volume_object_volume(csv_processor):
    fig, ax = plt.subplots()
    ax.set_title("Object volume")
    ax.set_xlabel("Frame")
    ax.set_ylabel("Volume \u00C5^3$")

    x = csv_processor.get_column_values("# frame")
    y = csv_processor.get_column_values("object_volume")

    ax.plot(x, y)

    ax.set_xlim((min(x), max(x)))
    ax.set_ylim((min(y), max(y)))

    return fig


# 8
def cluster_area(file_processor, suffix=""):
    fig, ax = plt.subplots()

    ax.set_title("Clusters area" + suffix)
    ax.set_xlabel("Density")
    ax.set_ylabel("Area")

    column_names = file_processor.get_column_names("Clusters summary - areas" + suffix)
    clusters_id = file_processor.get_column_values("Clusters summary - areas" + suffix, "Cluster")

    x = [name for name in column_names if name.startswith("D")]
    x.sort(key=lambda x: int(x[1:]), reverse=True)

    # Values for each density
    column_values = [file_processor.get_column_values("Clusters summary - areas" + suffix, density) for density in x]

    # Transpose matrix to get values for each cluster
    y_set = list(zip(*column_values))

    for cluster_id, y in zip(clusters_id, y_set):
        ax.plot(x, y, label="Cluster {}".format(cluster_id))

    ax.set_xlim((0, len(x) - 1))
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    return fig


class StateFrame(tk.Frame):
    def __init__(self, parent, *args, **kw):
        tk.Frame.__init__(self, parent, *args, **kw)

        self.enabled = False

    def toggle(self):
        if self.enabled:
            self.enabled = False
            self.disable()
        else:
            self.enabled = True
            self.enable()

    def disable(self):
        self.enabled = False
        for child in self.winfo_children():
            child.configure(state=tk.DISABLED)

    def enable(self):
        self.enabled = True
        for child in self.winfo_children():
            child.configure(state=tk.NORMAL)


class Kraken(object):
    def __init__(self, parent):
        """
        Kraken App
        :param parent: Parent widget
        """
        self.parent = parent
        self.size = (600, 650)
        self.title = "Kraken"

        parent.title(self.title)
        parent.geometry("{}x{}".format(*self.size))

        logo = tk.PhotoImage(file=get_img("logo.gif"))

        logo_label = tkinter.ttk.Label(self.parent, image=logo, padding=-2)
        logo_label.image = logo
        logo_label.pack(padx=20, pady=20)

        # Used to auto positioning depending on length of version string
        version = "ver. " + aquaduct.version_nice()
        tkinter.ttk.Label(logo_label, text=version, background="white").place(relx=1, rely=1, x=-len(version) * 7, y=-20)

        self.main_frame = utils.VerticalScrolledFrame(self.parent)
        self.main_frame.pack(expand=1, fill="both")
        self.main_frame = self.main_frame.interior

        # General variables
        self.data_file = tk.StringVar()
        self.csv_file = tk.StringVar()
        self.dat_file = tk.StringVar()
        self.results_file = tk.StringVar()

        # Checkbuttons variables
        self.v1 = tk.BooleanVar()
        self.v2 = tk.BooleanVar()
        self.v3 = tk.BooleanVar()
        self.v4 = tk.BooleanVar()
        self.v7 = tk.BooleanVar()
        self.v8 = tk.BooleanVar()

        ###
        results_frame = tk.Frame(self.main_frame)
        results_frame.columnconfigure(0, weight=1)
        results_frame.columnconfigure(1, weight=1)
        results_frame.columnconfigure(2, weight=1)
        results_frame.pack(fill=tk.X, padx=100, pady=20)

        tkinter.ttk.Label(results_frame, text="Output file: ").grid(sticky="e", row=2, column=0)
        tkinter.ttk.Entry(results_frame, textvariable=self.results_file).grid(sticky="we", row=2, column=1)
        rload = tkinter.ttk.Button(results_frame, text="Load file", style="File.TButton")
        rload.grid(sticky="w", row=2, column=2)
        rload.bind("<Button-1>", lambda e: self.load_file(self.results_file))

        # Container for plots which use txt
        container_data = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        container_data.columnconfigure(0, weight=1)
        container_data.columnconfigure(1, weight=1)
        container_data.columnconfigure(2, weight=1)
        container_data.pack(fill=tk.X, padx=100, pady=10)

        tkinter.ttk.Label(container_data, text="Data file: ").grid(sticky="e", row=0, column=0)
        tkinter.ttk.Entry(container_data, textvariable=self.data_file).grid(sticky="we", row=0, column=1)
        data_load = tkinter.ttk.Button(container_data, text="Load file", style="File.TButton")
        data_load.grid(sticky="w", row=0, column=2)
        data_load.bind("<Button-1>", lambda e: self.load_file(self.data_file))

        # Container for plots which use csv
        container_csv = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        container_csv.columnconfigure(0, weight=1)
        container_csv.columnconfigure(1, weight=1)
        container_csv.columnconfigure(2, weight=1)
        container_csv.pack(fill=tk.X, padx=100, pady=10)

        tkinter.ttk.Label(container_csv, text="CSV file: ").grid(sticky="e", row=0, column=0)
        tkinter.ttk.Entry(container_csv, textvariable=self.csv_file).grid(sticky="we", row=0, column=1)
        csv_load = tkinter.ttk.Button(container_csv, text="Load file", style="File.TButton")
        csv_load.grid(sticky="w", row=0, column=2)
        csv_load.bind("<Button-1>", lambda e: self.load_file(self.csv_file))

        ### 1
        self.all1 = tk.BooleanVar(value=1)
        self.molecules1 = tk.StringVar()

        option1_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option1_frame.columnconfigure(0, weight=1)
        option1_frame.columnconfigure(1, weight=1)
        option1_frame.grid(sticky="ew", row=1, column=0, columnspan=3, padx=10, pady=10)

        cb1 = tkinter.ttk.Checkbutton(option1_frame, text="Clusters size", var=self.v1)
        cb1.pack(anchor="w")

        state1_frame = StateFrame(option1_frame)
        state1_frame.columnconfigure(0, weight=1)
        state1_frame.columnconfigure(1, weight=1)
        state1_frame.pack()

        cb1.configure(command=lambda: state1_frame.toggle())

        # Options
        tkinter.ttk.Label(state1_frame, text="For all").grid(sticky="w", row=0, column=1)
        tkinter.ttk.Checkbutton(state1_frame, var=self.all1).grid(sticky="e", row=0, column=0)

        tkinter.ttk.Label(state1_frame, text="For molecules: ").grid(row=1, column=0)
        tk.Entry(state1_frame, textvariable=self.molecules1).grid(row=1, column=1)

        state1_frame.disable()

        ### 3
        self.all3 = tk.BooleanVar(value=1)
        self.molecules3 = tk.StringVar()

        option3_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option3_frame.columnconfigure(0, weight=1)
        option3_frame.columnconfigure(1, weight=1)
        option3_frame.grid(sticky="ew", row=2, column=0, columnspan=3, padx=10, pady=10)

        cb3 = tkinter.ttk.Checkbutton(option3_frame, text="Molecule entry time distribution", var=self.v3)
        cb3.pack(anchor="w")

        state3_frame = StateFrame(option3_frame)
        state3_frame.columnconfigure(0, weight=1)
        state3_frame.columnconfigure(1, weight=1)
        state3_frame.pack()

        cb3.configure(command=lambda: state3_frame.toggle())

        # Options
        tkinter.ttk.Label(state3_frame, text="For all").grid(sticky="w", row=0, column=1)
        tkinter.ttk.Checkbutton(state3_frame, var=self.all3).grid(sticky="e", row=0, column=0)

        tkinter.ttk.Label(state3_frame, text="For molecules: ").grid(row=1, column=0)
        tk.Entry(state3_frame, textvariable=self.molecules3).grid(row=1, column=1)

        state3_frame.disable()

        ### 4
        self.chord_threshold = tk.StringVar()
        self.clusters_info4 = tk.StringVar()
        option4_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option4_frame.columnconfigure(0, weight=1)
        option4_frame.columnconfigure(1, weight=1)
        option4_frame.grid(sticky="ew", row=3, column=0, columnspan=3, padx=10, pady=10)

        cb4 = tkinter.ttk.Checkbutton(option4_frame, text="Intramolecular flows", var=self.v4)
        cb4.pack(anchor="w")

        state4_frame = StateFrame(option4_frame)
        state4_frame.columnconfigure(0, weight=1)
        state4_frame.columnconfigure(1, weight=1)
        state4_frame.pack()

        cb4.configure(command=lambda: state4_frame.toggle())

        tkinter.ttk.Label(state4_frame, text="Clusters info: ").grid(row=0, column=0)

        tk.Entry(state4_frame, textvariable=self.clusters_info4).grid(row=0, column=1)
        clust_info_load = tkinter.ttk.Button(state4_frame, command=lambda: self.load_file(self.clusters_info4),
                                     text="Load file", style="File.TButton")
        clust_info_load.grid(sticky="w", row=0, column=3)

        tkinter.ttk.Label(state4_frame, text="Threshold: ").grid(row=1, column=0)
        tk.Entry(state4_frame, textvariable=self.chord_threshold).grid(row=1, column=1)

        state4_frame.disable()

        ### 8
        self.all8 = tk.StringVar(value=1)
        self.molecules8 = tk.StringVar()

        option8_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option8_frame.columnconfigure(0, weight=1)
        option8_frame.columnconfigure(1, weight=1)
        option8_frame.grid(sticky="ew", row=4, column=0, columnspan=3, padx=10, pady=10)

        cb8 = tkinter.ttk.Checkbutton(option8_frame, text="Clusters area", var=self.v8)
        cb8.pack(anchor="w")

        state8_frame = StateFrame(option8_frame)
        state8_frame.columnconfigure(0, weight=1)
        state8_frame.columnconfigure(1, weight=1)
        state8_frame.pack()

        cb8.configure(command=lambda: state8_frame.toggle())

        # Options
        tkinter.ttk.Label(state8_frame, text="For all").grid(sticky="w", row=0, column=1)
        tkinter.ttk.Checkbutton(state8_frame, var=self.all8).grid(sticky="e", row=0, column=0)

        tkinter.ttk.Label(state8_frame, text="For molecules: ").grid(row=1, column=0)
        tk.Entry(state8_frame, textvariable=self.molecules8).grid(row=1, column=1)

        state8_frame.disable()

        ### 2
        self.clusters_info2 = tk.StringVar()

        option2_frame = tk.Frame(container_csv, bd=1, relief=tk.GROOVE)
        option2_frame.grid(sticky="ew", row=1, column=0, columnspan=3, padx=10, pady=10)
        cb2 = tkinter.ttk.Checkbutton(option2_frame, text="Relative cluster flows", var=self.v2)
        cb2.pack(anchor="w")

        state2_frame = StateFrame(option2_frame)
        state2_frame.columnconfigure(0, weight=1)
        state2_frame.columnconfigure(1, weight=1)
        state2_frame.pack()

        tkinter.ttk.Label(state2_frame, text="Clusters info: ").grid(row=0, column=0)
        clust_info_load2 = tkinter.ttk.Button(state2_frame, command=lambda: self.load_file(self.clusters_info2),
                                     text="Load file", style="File.TButton")
        clust_info_load2.grid(sticky="w", row=0, column=3)

        tk.Entry(state2_frame, textvariable=self.clusters_info2).grid(row=0, column=1)

        cb2.configure(command=lambda: state2_frame.toggle())

        state2_frame.disable()

        ### 7
        option7_frame = tk.Frame(container_csv, bd=1, relief=tk.GROOVE)
        option7_frame.grid(sticky="ew", row=2, column=0, columnspan=3, padx=10, pady=10)
        tkinter.ttk.Checkbutton(option7_frame, text="Volume per frame", var=self.v7).pack(anchor="w")

        ###
        self.generate_button = tkinter.ttk.Button(self.main_frame, text="Generate")
        self.generate_button.pack(pady=15)
        self.generate_button.bind("<Button-1>", lambda e: self.generate())

    def generate(self):
        if not True in [self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v7.get(), self.v8.get()]:
            tkinter.messagebox.showerror("Error", "No plot have been selected.")
            return

        if not self.results_file.get():
            tkinter.messagebox.showerror("Error", "You must specify results file.")
            return

        # Loading files
        traced_molecules = []
        if self.v1.get() or self.v3.get() or self.v4.get() or self.v8.get():
            if not self.data_file.get():
                tkinter.messagebox.showerror("Error", "Data file is not specified.")
                return

            try:
                f = FileDataProcessor(self.data_file.get())
            except Exception as e:
                tkinter.messagebox.showerror("Error", "Data file could not be opened.")
                print(e)

            # Fetching traced molecules names
            for line in f.file_lines:
                if line.startswith("Names of traced molecules:"):
                    traced_molecules = line.lstrip("Names of traced molecules:").split()

        if self.v2.get() or self.v7.get():
            if not self.csv_file.get():
                tkinter.messagebox.showerror("Error", "CSV file is not specified.")
                return

            try:
                c = CSVDataProcessor(self.csv_file.get())
            except Exception as e:
                tkinter.messagebox.showerror("Error", "CSV file could not be opened.")
                print(e)

        # Console log init
        log_window = tk.Toplevel(self.parent, width=100)
        log_console = tk.Text(log_window)
        log_console.pack(side=tk.LEFT, fill=tk.BOTH)

        log_scroll = tk.Scrollbar(log_window, command=log_console.yview)
        log_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        log_console["yscrollcommand"] = log_scroll.set

        self.generate_button.config(state=tk.DISABLED)
        self.generate_button.unbind("<Button-1>")

        self.parent.update_idletasks()  # Force console log to show before generating plots

        def close_log_window():
            log_window.destroy()
            self.generate_button.config(state=tk.NORMAL)
            self.generate_button.bind("<Button-1>", lambda e: self.generate())

        log_window.protocol("WM_DELETE_WINDOW", close_log_window)

        # Log console settings
        log_console.tag_config("error", foreground="red")
        log_console.tag_config("success", foreground="green")

        def log(*args, **kwargs):
            log_console.config(state=tk.NORMAL)
            log_console.insert(*args, **kwargs)
            log_console.config(state=tk.DISABLED)

        # Plotting
        plots = []

        if self.v1.get():
            log(tk.END, "{}\nGenerating Clusters size\n".format("-" * 30))

            if self.all1.get():
                log(tk.END, "* All ")

                try:
                    plot = BytesIO()
                    cluster_inlets(f).savefig(plot, format="png", bbox_inches="tight")
                    plots.append(plot)

                    log(tk.END, "\u2714\n", "success")
                except DataException as e:
                    log(tk.END, "\u2718\n", "error")
                    log(tk.END, "{}\n".format(e), "error")

            if self.molecules1.get():
                molecules = self.molecules1.get().replace(" ", "").upper().split(",")
                for molecule in molecules:
                    log(tk.END, "* {} ".format(molecule))

                    try:
                        plot = BytesIO()
                        cluster_inlets(f, suffix=" of {}".format(molecule)).savefig(plot, format="png",
                                                                                    bbox_inches="tight")
                        plots.append(plot)

                        log(tk.END, "\u2714\n", "success")
                    except DataException as e:
                        log(tk.END, "\u2718\n", "error")
                        log(tk.END, "{}\n".format(e), "error")

            log(tk.END, "Done.\n", "success")

        if self.v2.get():
            log(tk.END, "{}\nGenerating Relative cluster flows\n".format("-" * 30))

            labels_file = {}
            colors_file = {}
            if self.clusters_info2.get():
                with open(self.clusters_info2.get(), "r") as file:
                    for line in file.readlines():
                        id_, name, color = line.rstrip().split("\t")
                        labels_file.update({id_: name})
                        colors_file.update({id_: color})

            # Finding all clusters
            clusters = []
            ids = []
            for column_name in c.column_names:
                match = re.match("amol_apaths_(\d+)_walk$", column_name)
                if match:
                    clusters.append(match.group(0))
                    ids.append(int(match.group(1)))

            ids, clusters = list(zip(*sorted(zip(ids, clusters))))

            labels = [labels_file.pop(str(i), i) for i in ids]
            cg = color_gen()
            colors = [colors_file.pop(str(i), next(cg)) for i in ids]

            plot = BytesIO()
            relative_clusters_flows(c, clusters, labels, colors).savefig(plot, format="png", bbox_inches="tight")
            plots.append(plot)

            log(tk.END, "Done.\n", "success")

        if self.v3.get():
            log(tk.END, "{}\nGenerating Molecule entry time distribution\n".format("-" * 30))
            if self.all3.get():
                plot = BytesIO()
                ligands_time(f).savefig(plot, format="png")
                plots.append(plot)

            if self.molecules3.get():
                molecules = self.molecules3.get().replace(" ", "").upper().split(",")
                for molecule in molecules:
                    log(tk.END, "* {} ".format(molecule))

                    plot = BytesIO()
                    ligands_time(f, molecule).savefig(plot, format="png", bbox_inches="tight")
                    plots.append(plot)

                    log(tk.END, "\u2714\n", "success")

            log(tk.END, "Done.\n", "success")

        if self.v4.get():
            log(tk.END, "{}\nGenerating Intramolecular flows\n".format("-" * 30))

            labels = {}
            colors = {}
            if self.clusters_info4.get():
                with open(self.clusters_info4.get(), "r") as file:
                    for line in file.readlines():
                        id_, name, color = line.rstrip().split("\t")
                        labels.update({id_: name})
                        colors.update({id_: color})

            threshold = float(self.chord_threshold.get()) if self.chord_threshold.get() else 0.0

            log(tk.END, "* Clusters sizes ")
            plot1 = BytesIO()
            chord_diagram_sizes(f, labels, colors).savefig(plot1, format="png", dpi=2 ** 7, bbox_inches="tight")
            log(tk.END, "\u2714\n", "success")

            log(tk.END, "* Clusters flows ")
            plot2 = BytesIO()
            chord_diagram_flows(f, labels, colors, threshold).savefig(plot2, format="png", dpi=2 ** 7,
                                                                      bbox_inches="tight")
            log(tk.END, "\u2714\n", "success")

            plots.extend([plot1, plot2])

            log(tk.END, "Done.\n", "success")

        if self.v7.get():
            log(tk.END, "{}\nGenerating Volumes\n".format("-" * 30))

            plot1 = BytesIO()
            volume_scope_area(c).savefig(plot1, format="png", bbox_inches="tight")

            plot2 = BytesIO()
            volume_scope_volume(c).savefig(plot2, format="png", bbox_inches="tight")

            plot3 = BytesIO()
            volume_object_area(c).savefig(plot3, format="png", bbox_inches="tight")

            plot4 = BytesIO()
            volume_object_volume(c).savefig(plot4, format="png", bbox_inches="tight")

            plots.extend([plot1, plot2, plot3, plot4])

            log(tk.END, "Done.\n", "success")

        if self.v8.get():
            log(tk.END, "{}\nGenerating Clusters areas\n".format("-" * 30))

            if self.all8.get():
                plot = BytesIO()
                cluster_area(f).savefig(plot, format="png", bbox_inches="tight")
                plots.append(plot)

            if self.molecules8.get():
                molecules = self.molecules8.get().replace(" ", "").upper().split(
                    ",") if self.molecules8.get() else traced_molecules
                for molecule in molecules:
                    log(tk.END, "* {} ".format(molecule))

                    try:
                        plot = BytesIO()
                        cluster_area(f, suffix=" of {}".format(molecule)).savefig(plot, format="png")
                        plots.append(plot)
                        log(tk.END, "\u2714\n", "success")
                    except DataException as e:
                        log(tk.END, "\u2718\n", "error")
                        log(tk.END, "{}\n".format(e), "error")

            log(tk.END, "Done.\n", "success")

        # Save plots to file
        logo_base64 = base64.b64encode(open(get_img("logo.gif"), "rb").read()).decode("UTF-8")
        html = """<div style=\"text-align:center; margin-bottom: 40px\"><a href=\"http://www.tunnelinggroup.pl\"><img src=\"data:image/gif;base64,{}\"></a></div>""".format(logo_base64)
        
        for i, f in enumerate(plots):
            f.seek(0)
            html += """<div style=\"text-align: center\">
    <img style=\"max-width:100%\" src=\"data:image/png;base64,{}\">
</div>""".format(base64.b64encode(f.getvalue()).decode("UTF-8"))

        html += "<div style=\"text-align:center\"><h3>Document generated by <span style=\"color:#38b2fa\">Kraken</span>.</h3></div>"

        with open(self.results_file.get(), "w+") as f:
            f.write(html)
            log(tk.END, "Plots saved to {}.\n".format(self.results_file.get()), "success")

    def load_file(self, var):
        try:
            with askopenfile("r") as f:
                var.set(f.name)
        except AttributeError:  # In case of cancel selecting file
            pass


if __name__ == "__main__":
    root = tk.Tk()
    root.configure(background="white")
    root.resizable(1, 1)

    aq_icon = tk.PhotoImage(file=get_img("icon.gif"))
    root.tk.call('wm', 'iconphoto', root._w, aq_icon)

    s = tkinter.ttk.Style()
    s.theme_use("clam")

    s.configure("TLabel", padding=5)
    s.configure("File.TButton", padding=0, font=("TkDefaultFont", 8))  # Loading file button

    app = Kraken(root)
    root.mainloop()
