#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018  Michał Banas
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

import Tkinter as tk
import base64
import csv
import re
import tkMessageBox
import ttk
from cStringIO import StringIO
from tkFileDialog import askopenfile

import matplotlib.pyplot as plt

import aquaduct.apps.valveconfig.utils as utils
from aquaduct.apps.chord import Chord


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
            return list(data[0].iterkeys())
        else:
            raise DataException("Table \"{}\" is empty.".format(table_name))

    def _find_table(self, table_name):
        for i, line in enumerate(self.file_lines):
            if line.strip() == table_name:
                return self._parse_column(i)

        raise DataException("\"{}\" does not exist.".format(table_name))

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

        self.column_names = list(self.csv_reader.next().iterkeys())
        self.seek_file()

    def get_column_values(self, column_name):
        values = [int(row[column_name]) for row in self.csv_reader]
        self.seek_file()

        return values

    def seek_file(self):
        """ Set file position at beginning and skip first row with column names """
        self.file.seek(0)
        self.csv_reader.next()

    def __del__(self):
        self.file.close()


class DATDataProcessor(object):
    def __init__(self, filename):
        self.file = open(filename, "r")
        self.file_lines = [line.strip() for line in self.file.readlines()]

    def get_column_values(self, column_name):
        i = self.get_column_names().index(column_name)

        values = []
        for x in self.file_lines[1:]:
            value = x.split("\t")[i]
            if value.isdigit():
                value = int(value)
            elif is_float(value):
                value = float(value)

            values.append(value)

        return values

    def get_column_names(self):
        return self.file_lines[0].split("\t")

    def __del__(self):
        self.file.close()


# 1
def cluster_inlets(file_processor, suffix=""):
    fig, ax = plt.subplots()

    cluster_no = file_processor.get_column_values("Clusters summary - inlets" + suffix, "Cluster")

    y_in = file_processor.get_column_values("Clusters summary - inlets" + suffix, "INCOMING")
    y_out = file_processor.get_column_values("Clusters summary - inlets" + suffix, "OUTGOING")

    # Added 1 to make place for previous bar
    rects_bottom = ax.bar(range(1, len(cluster_no) + 1), y_in, label="INCOMING")
    rects_top = ax.bar(range(1, len(cluster_no) + 1), y_out, bottom=y_in, label="OUTGOING")

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

    ax.set_title("Clusters inlets" + suffix)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Size")

    ax.set_xlim((0, len(cluster_no) + 1))
    ax.set_xticks(range(1, len(cluster_no) + 1))
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
    ax.set_ylabel("???")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    ax.xaxis.set_major_locator(plt.AutoLocator())
    ax.set_xlim((0, len(x)))
    ax.set_ylim((0, 1.05))  # added 0.05 because when there is 1, its not visible

    return fig


# 3
def ligands_time(file_processor, molecule=None):
    # FIXME: Some data have "nan"
    fig, ax = plt.subplots()

    plot_settings = dict(align="edge", height=1.0)

    title_suffix = "" if not molecule else " of {}".format(molecule)
    ax.set_title("Ligands in time" + title_suffix)
    ax.set_xlabel("Frame")
    ax.set_ylabel("Separate path ID")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    path_nr = file_processor.get_column_values("List of separate paths and properties", "Nr")

    ax.set_ylim((0, max(path_nr)))

    path_nr = [y - 1 for y in path_nr]

    beg_values = file_processor.get_column_values("List of separate paths and properties", "BeginF")
    inp_values = file_processor.get_column_values("List of separate paths and properties", "InpF")
    obj_values = file_processor.get_column_values("List of separate paths and properties", "ObjF")
    out_values = file_processor.get_column_values("List of separate paths and properties", "OutF")

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
    beg_values, inp_values, obj_values, out_values = zip(*sorted(zip(beg_values,
                                                                     inp_values,
                                                                     obj_values,
                                                                     out_values)))

    inps_left = beg_values
    inps_width = [width if width - 1 >= 0 else 0 for width in inp_values]

    objs_left = [sum(values) for values in zip(beg_values, inp_values)]
    objs_width = [width - 1 for width in obj_values]

    outs_left = [sum(values) for values in zip(beg_values, inp_values, obj_values)]
    outs_width = [width - 1 for width in out_values]

    ax.barh(path_nr, inps_width, left=inps_left, color=(1, 0, 0), label="Incoming", **plot_settings)
    ax.barh(path_nr, objs_width, left=objs_left, color=(0, 1, 0), label="Object", **plot_settings)
    ax.barh(path_nr, outs_width, left=outs_left, color=(0, 0, 1), label="Outgoing", **plot_settings)

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")

    return fig


# 4
def chord_diagram(file_processor, labels={}, colors={}, threshold=0.):
    fig, ax = plt.subplots(subplot_kw={"aspect": 1})
    ax.set_title("Chord diagram")
    ax.set_axis_off()
    ax.set_xlim(-110, 110)
    ax.set_ylim(-110, 110)

    clusters_ids = file_processor.get_column_values("Clusters summary - inlets", "Cluster")
    clusters_ids = [str(id_) for id_ in clusters_ids]
    clusters_sizes = file_processor.get_column_values("Clusters summary - inlets", "Size")

    flows = file_processor.get_column_values("Separate paths clusters types summary - mean lengths of paths", "CType")
    flows = [{flow.split(":")[0]: flow.split(":")[1]} for flow in flows]

    flows_sizes = file_processor.get_column_values("Separate paths clusters types summary - mean lengths of paths",
                                                   "Size")

    # N size
    n_size = 0
    for flow, size in zip(flows, flows_sizes):
        if "N" in next(flow.iteritems()):
            n_size += size

    clusters_ids.append("N")
    clusters_sizes.append(n_size)

    flows = [{"source": clusters_ids.index(next(flow.iterkeys())), "dest": clusters_ids.index(next(flow.itervalues())),
              "value": size} for flow, size in zip(flows, flows_sizes) if 1. * size / sum(flows_sizes) >= threshold]

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

    Chord(ax, 100, clusters_sizes, flows, labels, colors)

    return fig


# 7
def volume_timeframe(dat_processor):
    raise NotImplemented()


# 8
def cluster_area(file_processor, suffix=""):
    fig, ax = plt.subplots()

    ax.set_title("Cluster areas" + suffix)
    ax.set_xlabel("Density")
    ax.set_ylabel("???")

    column_names = file_processor.get_column_names("Clusters summary - areas" + suffix)

    x = [name for name in column_names if name.startswith("D")]
    x.sort(key=lambda x: int(x[1:]), reverse=True)

    y = [file_processor.get_column_values("Clusters summary - areas" + suffix, density)[0] for density in x]

    ax.set_xlim((0, len(x) - 1))
    ax.plot(x, y)

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


class Octopus(object):
    def __init__(self, parent):
        """
        Octopus App
        :param parent: Parent widget
        """
        self.parent = parent
        self.size = (600, 650)
        self.title = "Octopus"

        parent.title(self.title)
        parent.geometry("{}x{}".format(*self.size))

        logo = tk.PhotoImage(file="../aquaduct/apps/valveconfig/logo.gif")

        logo_label = ttk.Label(self.parent, image=logo, padding=-2)
        logo_label.image = logo
        logo_label.pack(padx=20, pady=20)

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

        ttk.Label(results_frame, text="Results file: ").grid(sticky="e", row=2, column=0)
        ttk.Entry(results_frame, textvariable=self.results_file).grid(sticky="we", row=2, column=1)
        rload = ttk.Button(results_frame, text="Load file", style="File.TButton")
        rload.grid(sticky="w", row=2, column=2)
        rload.bind("<Button-1>", lambda e: self.load_file(self.results_file))

        # Container for plots which use txt
        container_data = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        container_data.columnconfigure(0, weight=1)
        container_data.columnconfigure(1, weight=1)
        container_data.columnconfigure(2, weight=1)
        container_data.pack(fill=tk.X, padx=100, pady=10)

        ttk.Label(container_data, text="Data file: ").grid(sticky="e", row=0, column=0)
        ttk.Entry(container_data, textvariable=self.data_file).grid(sticky="we", row=0, column=1)
        data_load = ttk.Button(container_data, text="Load file", style="File.TButton")
        data_load.grid(sticky="w", row=0, column=2)
        data_load.bind("<Button-1>", lambda e: self.load_file(self.data_file))

        # Container for plots which use csv
        container_csv = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        container_csv.columnconfigure(0, weight=1)
        container_csv.columnconfigure(1, weight=1)
        container_csv.columnconfigure(2, weight=1)
        container_csv.pack(fill=tk.X, padx=100, pady=10)

        ttk.Label(container_csv, text="CVS file: ").grid(sticky="e", row=0, column=0)
        ttk.Entry(container_csv, textvariable=self.csv_file).grid(sticky="we", row=0, column=1)
        csv_load = ttk.Button(container_csv, text="Load file", style="File.TButton")
        csv_load.grid(sticky="w", row=0, column=2)
        csv_load.bind("<Button-1>", lambda e: self.load_file(self.csv_file))

        # Container for plots which use dat
        container_dat = tk.Frame(self.main_frame, bd=1, relief=tk.SUNKEN)
        container_dat.columnconfigure(0, weight=1)
        container_dat.columnconfigure(1, weight=1)
        container_dat.columnconfigure(2, weight=1)
        container_dat.pack(fill=tk.X, padx=100, pady=10)

        ttk.Label(container_dat, text="DAT file: ").grid(sticky="e", row=0, column=0)
        ttk.Entry(container_dat, textvariable=self.dat_file).grid(sticky="we", row=0, column=1)
        dat_load = ttk.Button(container_dat, text="Load file", style="File.TButton")
        dat_load.grid(sticky="w", row=0, column=2)
        dat_load.bind("<Button-1>", lambda e: self.load_file(self.dat_file))

        ### 1
        option1_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option1_frame.columnconfigure(0, weight=1)
        option1_frame.columnconfigure(1, weight=1)
        option1_frame.grid(sticky="ew", row=1, column=0, columnspan=3, padx=10, pady=10)

        cb1 = ttk.Checkbutton(option1_frame, text="Inlets per cluster", var=self.v1)
        cb1.pack(anchor="w")

        state1_frame = StateFrame(option1_frame)
        state1_frame.columnconfigure(0, weight=1)
        state1_frame.columnconfigure(1, weight=1)
        state1_frame.pack()

        cb1.configure(command=lambda: state1_frame.toggle())

        ttk.Label(state1_frame, text="Molecules: ").grid(row=0, column=0)
        self.molecules1 = tk.Entry(state1_frame)
        self.molecules1.grid(row=0, column=1)

        state1_frame.disable()

        ### 3
        option3_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option3_frame.columnconfigure(0, weight=1)
        option3_frame.columnconfigure(1, weight=1)
        option3_frame.grid(sticky="ew", row=2, column=0, columnspan=3, padx=10, pady=10)

        cb3 = ttk.Checkbutton(option3_frame, text="Ligands per time", var=self.v3)
        cb3.pack(anchor="w")

        state3_frame = StateFrame(option3_frame)
        state3_frame.columnconfigure(0, weight=1)
        state3_frame.columnconfigure(1, weight=1)
        state3_frame.pack()

        cb3.configure(command=lambda: state3_frame.toggle())

        ttk.Label(state3_frame, text="Molecules: ").grid(row=0, column=0)
        self.molecules3 = tk.Entry(state3_frame)
        self.molecules3.grid(row=0, column=1)

        state3_frame.disable()

        ### 4
        option4_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option4_frame.columnconfigure(0, weight=1)
        option4_frame.columnconfigure(1, weight=1)
        option4_frame.grid(sticky="ew", row=3, column=0, columnspan=3, padx=10, pady=10)

        cb4 = ttk.Checkbutton(option4_frame, text="Chord diagram", var=self.v4)
        cb4.pack(anchor="w")

        state4_frame = StateFrame(option4_frame)
        state4_frame.columnconfigure(0, weight=1)
        state4_frame.columnconfigure(1, weight=1)
        state4_frame.pack()

        cb4.configure(command=lambda: state4_frame.toggle())

        ttk.Label(state4_frame, text="Clusters info: ").grid(row=0, column=0)
        self.clusters_info = tk.Entry(state4_frame)
        self.clusters_info.grid(row=0, column=1)

        ttk.Label(state4_frame, text="Threshold: ").grid(row=1, column=0)
        self.chord_threshold = tk.Entry(state4_frame)
        self.chord_threshold.grid(row=1, column=1)

        state4_frame.disable()

        ### 8
        option8_frame = tk.Frame(container_data, bd=1, relief=tk.GROOVE)
        option8_frame.columnconfigure(0, weight=1)
        option8_frame.columnconfigure(1, weight=1)
        option8_frame.grid(sticky="ew", row=4, column=0, columnspan=3, padx=10, pady=10)

        cb8 = ttk.Checkbutton(option8_frame, text="Clusters area", var=self.v8)
        cb8.pack(anchor="w")

        state8_frame = StateFrame(option8_frame)
        state8_frame.columnconfigure(0, weight=1)
        state8_frame.columnconfigure(1, weight=1)
        state8_frame.pack()

        cb8.configure(command=lambda: state8_frame.toggle())

        ttk.Label(state8_frame, text="Molecules: ").grid(row=0, column=0)
        self.molecules8 = tk.Entry(state8_frame)
        self.molecules8.grid(row=0, column=1)

        state8_frame.disable()

        ### 2
        option2_frame = tk.Frame(container_csv, bd=1, relief=tk.GROOVE)
        option2_frame.grid(sticky="ew", row=1, column=0, columnspan=3, padx=10, pady=10)
        cb2 = ttk.Checkbutton(option2_frame, text="Relative cluster flows", var=self.v2)
        cb2.pack(anchor="w")

        state2_frame = StateFrame(option2_frame)
        state2_frame.columnconfigure(0, weight=1)
        state2_frame.columnconfigure(1, weight=1)
        state2_frame.pack()

        ttk.Label(state2_frame, text="Colors file: ").grid(row=0, column=0)
        self.colors_file2 = tk.Entry(state2_frame)
        self.colors_file2.grid(row=0, column=1)

        cb2.configure(command=lambda: state2_frame.toggle())

        state2_frame.disable()

        ### 7
        option7_frame = tk.Frame(container_dat, bd=1, relief=tk.GROOVE)
        option7_frame.grid(sticky="ew", row=1, column=0, columnspan=3, padx=10, pady=10)
        ttk.Checkbutton(option7_frame, text="Volume per frame", var=self.v7).pack(anchor="w")

        ###
        self.generate_button = ttk.Button(self.main_frame, text="Generate")
        self.generate_button.pack(pady=15)
        self.generate_button.bind("<Button-1>", lambda e: self.generate())

    def generate(self):
        if not True in [self.v1.get(), self.v2.get(), self.v3.get(), self.v4.get(), self.v7.get(), self.v8.get()]:
            tkMessageBox.showerror("Error", "No plot have been selected.")
            return

        if not self.results_file.get():
            tkMessageBox.showerror("Error", "You must specify results file.")
            return

        # Loading files
        traced_molecules = []
        if self.v1.get() or self.v3.get() or self.v4.get() or self.v8.get():
            if not self.data_file.get():
                tkMessageBox.showerror("Error", "Data file is not specified.")
                return

            try:
                f = FileDataProcessor(self.data_file.get())
            except Exception as e:
                tkMessageBox.showerror("Error", "Data file could not be opened.")
                print e

            # Fetching traced molecules names
            for line in f.file_lines:
                if line.startswith("Names of traced molecules:"):
                    traced_molecules = line.lstrip("Names of traced molecules:").split()

        if self.v2.get():
            if not self.csv_file.get():
                tkMessageBox.showerror("Error", "CSV file is not specified.")
                return

            try:
                c = CSVDataProcessor(self.csv_file.get())
            except Exception as e:
                tkMessageBox.showerror("Error", "CSV file could not be opened.")
                print e

        if self.v7.get():
            if not self.dat_file.get():
                tkMessageBox.showerror("Error", "DAT file is not specified.")
                return

            try:
                d = DATDataProcessor(self.dat_file)
            except Exception as e:
                tkMessageBox.showerror("Error", "DAT file could not be opened.")
                print e

        # Console log init
        log_window = tk.Toplevel(self.parent, width=100)
        log_console = tk.Text(log_window)
        log_console.pack(fill=tk.BOTH)

        self.generate_button.config(state=tk.DISABLED)
        self.generate_button.unbind("<Button-1>")

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
            log(tk.END, "{}\nGenerating inlets per cluster\n".format("-" * 30))

            if not self.molecules1.get() and len(traced_molecules) == 1:
                plot = StringIO()
                cluster_inlets(f).savefig(plot, format="png", bbox_inches="tight")
                plots.append(plot)
            else:
                molecules = self.molecules1.get().replace(" ", "").upper().split(
                    ",") if self.molecules1.get() else traced_molecules
                for molecule in molecules:
                    log(tk.END, "* {} ".format(molecule))

                    plot = StringIO()
                    cluster_inlets(f, suffix=" of {}".format(molecule)).savefig(plot, format="png", bbox_inches="tight")
                    plots.append(plot)

                    log(tk.END, u"\u2714\n", "success")

            log(tk.END, "Done.\n", "success")

        if self.v2.get():
            log(tk.END, "{}\nGenerating relative cluster flows\n".format("-" * 30))

            plot = StringIO()
            # Finding all clusters
            clusters = []
            ids = []
            for column_name in c.column_names:
                match = re.match("amol_apaths_(\d+)_walk$", column_name)
                if match:
                    clusters.append(match.group(0))
                    ids.append(int(match.group(1)))

            ids, clusters = zip(*sorted(zip(ids, clusters)))
            labels = ["Cluster " + str(i) for i in ids]

            # TODO: Fetching colors from file
            colors = ["#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0",
                      "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", "#9a6324",
                      "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080"]

            relative_clusters_flows(c, clusters, labels, colors).savefig(plot, format="png", bbox_inches="tight")
            plots.append(plot)

            log(tk.END, "Done.\n", "success")

        if self.v3.get():
            log(tk.END, "{}\nGenerating ligands per time\n".format("-" * 30))
            if not self.molecules3.get():
                plot = StringIO()
                ligands_time(f).savefig(plot, format="png", bbox_inches="tight")
                plots.append(plot)
            else:
                molecules = self.molecules3.get().replace(" ", "").upper().split(",")
                for molecule in molecules:
                    log(tk.END, "* {} ".format(molecule))

                    plot = StringIO()
                    ligands_time(f, molecule).savefig(plot, format="png", bbox_inches="tight")
                    plots.append(plot)

                    log(tk.END, u"\u2714\n", "success")

            log(tk.END, "Done.\n", "success")

        if self.v4.get():
            log(tk.END, "{}\nGenerating flows between tunnels\n".format("-" * 30))

            labels = {}
            colors = {}
            if self.clusters_info.get():
                with open(self.clusters_info.get(), "r") as file:
                    for line in file.readlines():
                        id_, name, color = line.rstrip().split("\t")
                        labels.update({id_: name})
                        colors.update({id_: color})

            threshold = float(self.chord_threshold.get()) if self.chord_threshold.get() else 0.0

            plot = StringIO()
            chord_diagram(f, labels, colors, threshold).savefig(plot, format="png", dpi=2 ** 7)
            plots.append(plot)

            log(tk.END, "Done.\n", "success")

        if self.v7.get():
            log(tk.END, "{}\nGenerating volume per time\n".format("-" * 30))
            log(tk.END, "Done.\n", "success")

        if self.v8.get():
            log(tk.END, "{}\nGenerating clusters areas\n".format("-" * 30))

            molecules = self.molecules8.get().replace(" ", "").upper().split(
                ",") if self.molecules8.get() else traced_molecules
            for molecule in molecules:
                log(tk.END, "* {} ".format(molecule))

                try:
                    plot = StringIO()
                    cluster_area(f, suffix=" of {}".format(molecule)).savefig(plot, format="png")
                    plots.append(plot)
                    log(tk.END, u"\u2714\n", "success")
                except DataException as e:
                    log(tk.END, u"\u2718\n", "error")
                    log(tk.END, "{}\n".format(e), "error")

            log(tk.END, "Done.\n", "success")

        # Save plots to file
        logo_base64 = base64.b64encode(open("../aquaduct/apps/valveconfig/logo.gif", "rb").read())
        html = """<div style="text-align:center; margin-bottom: 40px"><a href="http://www.tunnelinggroup.pl"><img src="data:image/gif;base64,{}"></a></div>""".format(
            logo_base64)
        for i, f in enumerate(plots):
            html += """<div style="text-align: center">
    <img src=\"data:image/png;base64,{}\">
</div>""".format(base64.encodestring(f.getvalue()))

        html += "<div style=\"text-align:center\"><h3>Document generated by <span style=\"color:#38b2fa\">Octopus</span>.</h3></div>"

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

    s = ttk.Style()
    s.theme_use("clam")

    s.configure("TLabel", padding=5)
    s.configure("File.TButton", padding=0, font=("TkDefaultFont", 8))  # Loading file button

    app = Octopus(root)
    root.mainloop()
