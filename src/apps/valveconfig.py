#!/usr/bin/env python3.7
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

import tkinter as tk
import hashlib
import traceback
import os
import shutil
import tempfile
import tkinter.messagebox
import tkinter.ttk
import webbrowser
from configparser import ConfigParser, NoOptionError, NoSectionError
from collections import OrderedDict, defaultdict
from tkinter.filedialog import askopenfile, askdirectory, asksaveasfile

import aquaduct
import aquaduct.apps.valveconfig.defaults as defaults
import aquaduct.apps.valveconfig.utils as utils
from aquaduct.apps.valveconfig import get_img


class ValveConfigApp(object):
    def __init__(self, parent):
        """ Valve Configurator App """
        self.size = (600, 650)

        self.parent = parent

        self.values = OrderedDict()

        # Used to hide that frame, after all entries are shown
        self.init_frame = tkinter.ttk.Frame(self.parent)
        self.notebook = tkinter.ttk.Notebook(self.parent)

        self.level = None

        self.hiding_frames = {}

        self.cluster_frame_index = None

        self.recluster_frame_index = None

        # Used for identifying scroll frame by interior frame for removing it from notebook
        self.scrolled_frames = {}

        self.frames = []

        # Dict with input frames of required entries. Used to change color if is unfilled
        self.required_entries = defaultdict(dict)

        self.config_filename = tk.StringVar()
        # Used to determine if max_level difference message should be displayed
        self.config_loaded = False

        self.values_hash = self.get_values_hash()

        self.title = "Valve Configurator"

        parent.title(self.title)
        parent.geometry("600x650")

        self.parent.protocol("WM_DELETE_WINDOW", self.on_window_close)

        self.init_gui()

    def init_gui(self):
        """ Prepare initial frame """
        # Logo
        logo = tk.PhotoImage(file=get_img("logo.gif"))

        logo_label = tkinter.ttk.Label(self.parent, image=logo, padding=-2)
        logo_label.image = logo
        logo_label.pack(padx=20, pady=20)

        # Used to auto positioning depending on length of version string
        version = "ver. " + aquaduct.version_nice()
        tkinter.ttk.Label(logo_label, text=version, background="white").place(relx=1, rely=1, x=-len(version) * 7, y=-20)

        # Frame with loading configuration file
        load_frame = tkinter.ttk.LabelFrame(self.init_frame, text="Configuration file")
        load_frame.columnconfigure(0, weight=1)
        load_frame.columnconfigure(1, weight=1)

        tkinter.ttk.Button(load_frame, text="Load config", command=self.open_config_file).grid(sticky="E", row=0, column=0)
        tkinter.ttk.Button(load_frame, text="New config", command=self.create_new_config_file_dialog).grid(sticky="W",
                                                                                                   row=0,
                                                                                                   column=1,
                                                                                                   padx=5,
                                                                                                   pady=2)
        tkinter.ttk.Entry(load_frame, textvariable=self.config_filename, state="readonly").grid(sticky="EW",
                                                                                        row=1,
                                                                                        column=0,
                                                                                        columnspan=2,
                                                                                        padx=5,
                                                                                        pady=2)

        # Level selection
        level_frame = tkinter.ttk.LabelFrame(self.init_frame, text="Configuration level")

        self.level = tk.StringVar()

        # List of sorted levels from easiest to hardest
        self.levels = []
        for level_name, _ in reversed(sorted(iter(defaults.LEVELS.items()), key=lambda k_v: (k_v[1], k_v[0]))):
            self.levels.append(level_name)

        tkinter.ttk.OptionMenu(level_frame, self.level, self.levels[0], *self.levels).pack(pady=20)

        load_frame.pack(fill=tk.X, padx=100, pady=20)
        level_frame.pack(fill=tk.X, padx=100, pady=20)

        def show_new_gui():
            self.prepare_section_frames()
            self.create_interface()

        tkinter.ttk.Button(self.init_frame, text="Forward", command=show_new_gui).pack(pady=50)

        self.init_frame.pack(expand=1, fill="both")

    def prepare_section_frames(self):
        """ Parse DEFAULTS and depends on it create Notebook tabs and entries """
        for section in defaults.DEFAULTS:
            if section._nested: # Skip nested frames which were added in runtime.
                continue

            section_name = section.config_name
            section_name_long = section.name
            entries = section.entries

            if section.abs_level is not None:
                if section.abs_level != defaults.LEVELS[self.level.get()]:
                    continue

            if defaults.LEVELS[self.level.get()] > section.level:
                continue

            if section_name == "clustering":
                self.cluster_frame_index = len(self.frames)

            if section_name == "reclustering":
                self.recluster_frame_index = len(self.frames)

            scroll_frame = utils.VerticalScrolledFrame(self.notebook)

            frame = scroll_frame.interior
            frame.columnconfigure(0, weight=1)
            frame.columnconfigure(1, weight=1)

            self.frames.append(frame)
            self.scrolled_frames[frame] = scroll_frame

            self.notebook.add(scroll_frame, text=section_name_long)

            tkinter.ttk.Label(frame,
                      text=section_name_long,
                      style="Title.TLabel",
                      background=utils.get_widget_bg(self.parent)).grid(row=0, column=0, columnspan=2)
            tkinter.ttk.Separator(frame, orient=tk.HORIZONTAL).grid(sticky="EW", row=1, column=0, columnspan=3)

            self.entry_filler(frame, section_name, entries)

        if self.cluster_frame_index:
            for section_name in self.get_recursive_clustering_sections("clustering"):
                self.append_entries(section_name)

                # Increment max_level
                self.add_max_level(1)

            cluster_add_button = tkinter.ttk.Button(self.frames[self.cluster_frame_index], text="Add clustering section")
            # Setting row=1000 let skip calculating position of button each time new section is added
            cluster_add_button.grid(row=1000, column=0, columnspan=2, pady=20)

            cluster_add_section_callback = utils.CallbackWrapper(self.callback_add_section, "clustering")
            cluster_add_button.bind("<Button-1>", cluster_add_section_callback)

        # if "reclustering" in self.hiding_frames:
        if self.recluster_frame_index:
            for section_name in self.get_recursive_clustering_sections("reclustering"):
                self.append_entries(section_name)

            recluster_add_button = tkinter.ttk.Button(self.frames[self.recluster_frame_index], text="Add reclustering section")
            # Setting row=1000 let skip calculating position of button each time new section is added
            recluster_add_button.grid(row=1000, column=0, columnspan=2, pady=20)

            recluster_add_section_callback = utils.CallbackWrapper(self.callback_add_section, "reclustering")
            recluster_add_button.bind("<Button-1>", recluster_add_section_callback)

        # Refresh all hiding frames
        self.refresh_menus()

        self.init_frame.pack_forget()

        self.notebook.pack(expand=1, fill="both")

        # After preparing all section load values from config
        if self.config_filename.get() != "":
            self.load_config_values(self.config_filename.get())

        self.values_hash = self.get_values_hash()

    def create_interface(self):
        """ Creates window menu """
        menu_bar = tk.Menu(self.parent)

        # Level menu
        level_menu = tk.Menu(menu_bar)

        def change_level_callback(level_name):
            self.level.set(level_name)
            self.recreate_gui()

        for level_name in self.levels:
            level_cb = utils.CallbackWrapper(change_level_callback, level_name)
            level_menu.add_command(label=level_name, command=level_cb)

        # File menu
        file_menu = tk.Menu(menu_bar, tearoff=0)
        file_menu.add_command(label="New", command=self.create_new_config_file_dialog)

        def open_config():
            if self.open_config_file():
                self.load_config_values(self.config_filename.get())

                self.values_hash = self.get_values_hash()

        file_menu.add_command(label="Open", command=open_config)

        def save_callback(*args):
            if self.config_filename.get() == "":
                f = asksaveasfile("w")
                f.close()
                if f and self.save_config(f.name):
                    self.parent.title("{} - {}".format(self.title, f.name))
                    self.config_filename.set(f.name)
                else:
                    os.remove(f.name)
            else:
                self.save_config(self.config_filename.get())

        file_menu.add_command(label="Save", command=save_callback)

        def save_as_callback(*args):
            f = asksaveasfile("w")
            if f and self.save_config(f.name) and self.config_filename.get() == "":
                self.parent.title("{} - {}".format(self.title, f.name))
                self.config_filename.set(f.name)

        file_menu.add_command(label="Save as", command=save_as_callback)

        file_menu.add_separator()

        reset_cb = utils.CallbackWrapper(self.load_defaults, reset=True)
        file_menu.add_command(label="Reset to default", command=reset_cb)
        file_menu.add_cascade(label="Level", menu=level_menu)
        file_menu.add_command(label="Exit", command=self.on_window_close)
        menu_bar.add_cascade(label="File", menu=file_menu)

        # Run menu
        run_menu = tk.Menu(menu_bar, tearoff=0)
        run_menu.add_command(label="Valve", command=self.valve_run_dialog)
        # run_menu.add_command(label="Pond", command=self.pond_run_dialog)
        menu_bar.add_cascade(label="Run", menu=run_menu)

        # Help menu
        help_menu = tk.Menu(menu_bar, tearoff=0)
        help_menu.add_command(label="Documentation", command=self.open_docs)
        help_menu.add_command(label="About", command=self.about)
        menu_bar.add_cascade(label="Help", menu=help_menu)

        self.parent.config(menu=menu_bar)

    def callback_add_section(self, section_name):
        """
        Callback for button to create new recursive clustering in clustering/reclustering frame

        :param section_name: Config section name where new frames will be appended, allowed are "clustering" or "reclustering"
        """
        # Find free name for new section
        index = 0
        while section_name + str(index) in self.values:
            index += 1

        self.append_entries(section_name + str(index))

        if section_name.startswith("clustering"):
            self.add_max_level(1)

        self.refresh_menus()

    def callback_remove_section(self, section_name, frame):
        """
        Callback for button to remove existing recursive clustering in clustering/reclustering frame

        :param section_name: Config section name to remove, allowed are "clustering" or "reclustering"
        :param frame: frame Which contains options related to that section
        """
        del self.values[section_name]
        del self.hiding_frames[section_name]

        # Remove from MENUS list
        for i, item in enumerate(defaults.MENUS):
            if item.startswith(section_name):
                del defaults.MENUS[i]

        if section_name.startswith("clustering"):
            self.add_max_level(-1)

        frame.grid_forget()

    def get_recursive_clustering_sections(self, section_name):
        """
        Finds recursively all section used in recursive_clustering options

        :param section_name: Config section name from which fetching will start.
        :return Uniqe section names without section_name param.
        :rtype: list
        """
        if section_name != "clustering" and section_name != "reclustering":
            raise RuntimeError(
                "There is no possibility to get recursive clustering from {} section".format(section_name))

        try:
            with open(self.config_filename.get(), "r") as config_file:
                config = ConfigParser()
                config.readfp(config_file)

        except IOError:  # File is empty
            return []

        clustering_sections = []
        try:
            clustering_section = config.get(section_name, "recursive_clustering")
            while clustering_section not in clustering_sections and clustering_section != "None":
                clustering_sections.append(clustering_section)
                clustering_section = config.get(clustering_section, "recursive_clustering")
        except (NoOptionError, NoSectionError):
            pass

        return list(set(clustering_sections).difference([section_name]))

    def append_entries(self, section_name):
        """
        Append new frame with new clustering or reclustering options.

        :param section_name: Config section name where new frames will be appended, allowed are "clustering" or "reclustering".
        """
        if section_name.startswith("clustering") and not self.cluster_frame_index:
            return

        if section_name.startswith("reclustering") and not self.recluster_frame_index:
            return

        if section_name.startswith("clustering"):
            default_section_name = "clustering"
            frame = self.frames[self.cluster_frame_index]
        elif section_name.startswith("reclustering"):
            default_section_name = "reclustering"
            frame = self.frames[self.recluster_frame_index]
        else:
            raise RuntimeError(
                "Appending entries to sections other than clustering or reclusteriation is not allowed")

        # Calculate row
        row = len(frame.grid_slaves(None, 0)) + 2

        # Set option menu from recursive clustering to control hiding frames
        defaults.MENUS.append("{}:method".format(section_name))

        default_entries = defaults.get_default_section(default_section_name)

        inner_frame = tkinter.ttk.Frame(frame, style="I.TFrame")
        inner_frame.columnconfigure(0, weight=1)
        inner_frame.columnconfigure(1, weight=1)
        inner_frame.grid(sticky="EW", row=row, column=0, columnspan=2, padx=60, pady=10)

        tkinter.ttk.Label(inner_frame, text="Recursive clustering", style="ITitle.TLabel").grid(
            row=0,
            column=0,
            pady=5,
            columnspan=2)
        tkinter.ttk.Separator(inner_frame, orient=tk.HORIZONTAL).grid(sticky="EW", row=1, column=0, columnspan=2)

        self.entry_filler(inner_frame, section_name, default_entries.entries)

        remove_button = tkinter.ttk.Button(inner_frame, text="Remove")
        remove_button.grid(row=999, column=0, columnspan=2, pady=5)

        remove_section_callback = utils.CallbackWrapper(self.callback_remove_section, section_name, inner_frame)
        remove_button.bind("<Button-1>", remove_section_callback)

    def entry_filler(self, parent, section_name, entries):
        """
        Creates entries depending on DEFAULTS values.

        When no entry added it hide notebook tab.

        Disables name option for clustering, reclustering and derivatives section.

        :param parent: Parent widget.
        :param section_name: Config section name.
        :param entries: List of entries.

        If section is nested, it's level is not checked. It allow to nest section that have higher level,
        but it cannot be visible in notebook.
        """
        # Keeps label frames with actual row
        group_frames = {}

        _, row = parent.grid_size()

        entries_appended = 0
        for entry in entries:
            # If its primary reclustering or clustering section disable changing name of that section
            # Additionaly it disable inlets_clustering:max_level if config level is higher than max_level level
            state = tk.NORMAL
            if section_name == "clustering" and entry.config_name == "name":
                state = tk.DISABLED
            elif section_name == "reclustering" and entry.config_name == "name":
                state = tk.DISABLED
            elif section_name == "inlets_clustering" and entry.config_name == "max_level":
                if defaults.LEVELS[self.level.get()] != entry.level:
                    state = tk.DISABLED

            # Skip entries that dont match level except inlets_clustering:max_level
            # and are not nested section
            if defaults.LEVELS[self.level.get()] > entry.level and not isinstance(entry, defaults.DefaultSection):
                if not (section_name == "inlets_clustering" and entry.config_name == "max_level"):
                    continue

            entries_appended += 1

            if isinstance(entry, defaults.DefaultSection):
                nested_frame = tkinter.ttk.LabelFrame(parent, text=entry.name)
                nested_frame.grid_columnconfigure(0, weight=1)
                nested_frame.grid_columnconfigure(1, weight=1)

                """
                Adding in runtime section to DEFAULTS sections to let them be visible by 
                defaults.get_default_entry or defaults.get_default_section functions during 
                reloading values when level was changed.
                """
                entry._nested = True
                defaults.DEFAULTS.append(entry)
                self.entry_filler(nested_frame, entry.config_name, entry.entries)

                if nested_frame.grid_size()[1]:
                    nested_frame.grid(row=parent.grid_size()[1], column=0, columnspan=2, pady=15, ipadx=30)
            elif entry.optionmenu_value:
                if section_name not in self.hiding_frames:
                    self.hiding_frames[section_name] = {}

                if entry.optionmenu_value not in self.hiding_frames[section_name]:
                    self.hiding_frames[section_name][entry.optionmenu_value] = utils.HidingFrame(parent,
                                                                                                 parent.grid_size()[1],
                                                                                                 text=entry.optionmenu_value.capitalize() + " options",
                                                                                                 style="HF.TFrame")

                    # It must be showed first, to allocate space in frame, otherwise grid_size wont return proper value
                    self.hiding_frames[section_name][entry.optionmenu_value].show()

                hiding_frame = self.hiding_frames[section_name][entry.optionmenu_value]

                if section_name not in self.values:
                    self.values[section_name] = {}

                if entry.optionmenu_value not in self.values[section_name]:
                    self.values[section_name][entry.optionmenu_value] = {}

                self.values[section_name][entry.optionmenu_value][entry.config_name] = utils.entry_factory(hiding_frame,
                                                                                                           hiding_frame.grid_size()[1],
                                                                                                           entry.name,
                                                                                                           entry.default_values,
                                                                                                           entry.help_text,
                                                                                                           state,
                                                                                                           info_text=entry.info_text,
                                                                                                           warning_text=entry.warning_text)

                if entry.required:
                    self.required_entries[section_name][entry.config_name] = self.values[section_name][
                        entry.config_name]
            else:
                if entry.group_label:
                    if entry.group_label not in group_frames:
                        label_frame = tkinter.ttk.LabelFrame(parent, text=entry.group_label)
                        label_frame.grid_columnconfigure(0, weight=1)
                        label_frame.grid_columnconfigure(1, weight=1)
                        label_frame.grid(row=parent.grid_size()[1], column=0, columnspan=2, pady=15, ipadx=30)

                        group_frames[entry.group_label] = [label_frame, 0]

                    entry_parent = group_frames[entry.group_label][0]
                    entry_row = group_frames[entry.group_label][1]

                    group_frames[entry.group_label][1] += 1
                else:
                    entry_parent = parent

                entry_widget = utils.entry_factory(entry_parent,
                                                   entry_parent.grid_size()[1],
                                                   entry.name,
                                                   entry.default_values,
                                                   entry.help_text,
                                                   state,
                                                   info_text=entry.info_text,
                                                   warning_text=entry.warning_text)

                if (section_name.startswith("clustering") or
                    section_name.startswith("reclustering")) and entry.config_name == "name":
                    entry_widget.set(section_name)

                if defaults.is_menu(section_name, entry.config_name):
                    opt_menu_callback = utils.CallbackWrapper(self.option_menu_changed, section_name, entry_widget)
                    entry_widget.input_var.trace("w", opt_menu_callback)

                if section_name not in self.values:
                    self.values[section_name] = {}

                self.values[section_name][entry.config_name] = entry_widget

                if entry.required:
                    self.required_entries[section_name][entry.config_name] = entry_widget

        # Hide tab if unused
        if entries_appended == 0:
            if parent in self.scrolled_frames:  # Only for sections that are part of notebook
                self.notebook.forget(self.scrolled_frames[parent])

    def refresh_menus(self):
        """
        Sets visible only hiding frames that are chosen by option menu.
        """
        for menu in defaults.MENUS:
            section_name, entry_name = menu.split(":")

            if self.values.get(section_name, None):
                self.option_menu_changed(section_name, self.values[section_name][entry_name])

    def re_clustering_exists(self, section_name):
        """
        Returns true if clustering or reclustering sections exists

        It omit ability to create recursive sections

        :param section_name: Section name
        :type str
        :return: True if section exists, otherwise False
        :rtype bool
        """
        return section_name in self.hiding_frames

    def option_menu_changed(self, section, entry):
        """
        Callback for option menu that control hiding frames

        :param section: Section name where hiding frame is located.
        :param entry: Option menu widget, which control hiding frames in section.
        """
        for method, hidden_frame in self.hiding_frames[section].items():
            if method == entry.get():
                hidden_frame.show()
            else:
                hidden_frame.grid_forget()

    def add_max_level(self, value):
        """
        Add value to max_level

        :param value: Number which will be added to max_level.
        """
        current_value = self.values["inlets_clustering"]["max_level"].get()
        self.values["inlets_clustering"]["max_level"].set(current_value + value)

    def open_config_file(self):
        """ Show dialog to choose config file and save its name """
        filetypes = (("text files", "*.txt"), ("config files", "*.cfg"), ("all files", "*.*"))
        try:
            with askopenfile("r", filetypes=filetypes) as config_file:
                self.parent.title("{} - {}".format(self.title, config_file.name))
                self.config_filename.set(config_file.name)
        except AttributeError:
            return False

        return True

    def recreate_gui(self, *args):
        """
        Recreate GUI.

        Saves current values to temporary file.

        :param args: Event informations.
        """
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfile.close()
        tmp_filename = tmpfile.name

        # Save values to temporary file.
        self.save_config(tmpfile.name, required_checking=False)

        self.cluster_frame_index = None
        self.recluster_frame_index = None

        for tab_id in self.notebook.tabs():
            self.notebook.forget(tab_id)

        self.hiding_frames = {}

        self.prepare_section_frames()
        self.load_defaults()  # Reset values
        self.load_config_values(tmp_filename)

        os.remove(tmp_filename)

    def get_values_hash(self):
        """ Compute hash from values. """
        hash = hashlib.md5()
        hash.update((str(self.values)).encode('utf-8'))
        # for section_name in self.values:
        #     for value in self.values[section_name].itervalues():
        #         hash.update(str(value.get()))

        return hash.digest()

    def load_config_values(self, config_filename):
        """
        Load config values to values dictionary from file.

        Additionally it ask user if inlets_clustering:max_level in config
        is different from number of found recursive_clustering sections.

        :param config_filename: Name of file with configuration.
        """
        # Config can be loaded when UI was created.
        # Checking if there is no additional clustering or reclustering sections and add them.
        if self.re_clustering_exists("clustering"):
            for section_name in self.get_recursive_clustering_sections("clustering"):
                if not self.re_clustering_exists(section_name):
                    self.append_entries(section_name)

        with open(config_filename, "r") as config_file:
            config = ConfigParser()
            config.readfp(config_file)

            for section_name in self.values:
                for entry_name in self.values[section_name]:

                    # Additional for is required to handle clustering/reclustering methods
                    entries = [entry_name] if not isinstance(self.values[section_name][entry_name], dict) else self.values[section_name][entry_name]
                    for entry_name in entries:
                        try:
                            config_value = config.get(section_name, entry_name)

                            #  Inform user that config max_level differ from number of found recursive clustering
                            result = True
                            if section_name == "inlets_clustering" and entry_name == "max_level" and self.config_loaded:
                                if self.values["inlets_clustering"]["max_level"].get() != int(config_value):
                                    result = tkinter.messagebox.askyesno("Max level",
                                                                   "Config max_level value differs from number of found clustering sections.\n"
                                                                   "Do you want to keep value from configuration file?")

                            if not result:
                                continue

                            optionmenu_value = defaults.get_default_entry(section_name, entry_name).optionmenu_value
                            if optionmenu_value:
                                try:
                                    method = config.get(section_name, "method")
                                except NoOptionError:
                                    method = defaults.get_default_entry(section_name, "method")

                                self.values[section_name][method][entry_name].set(config_value)
                            else:
                                self.values[section_name][entry_name].set(config_value)
                        except NoOptionError:
                            continue
                        except NoSectionError:
                            break

        self.config_loaded = True

    def load_defaults(self, reset=False, *args):
        """
        Load defaults values to values dictionary.

        :param reset: If set to True confirmation will be needed to reset.
        """
        if reset:
            result = tkinter.messagebox.askyesno("Reset values", "Are you sure to reset all values?")

            if not result:
                return

        for section in defaults.DEFAULTS:
            if defaults.LEVELS[self.level.get()] > section.level:
                continue

            for entry_section, entry in section.iter_entries():
                # Skip entries that dont match level except inlets_clustering:max_level
                if defaults.LEVELS[self.level.get()] > entry.level:
                    if not (entry_section.config_name == "inlets_clustering" and entry.config_name == "max_level"):
                        continue

                if entry.optionmenu_value:
                    self.values[entry_section.config_name][entry.optionmenu_value][entry.config_name].set(entry.default_value)
                else:
                    self.values[entry_section.config_name][entry.config_name].set(entry.default_value)

        # Delete appended clustering and reclustering sections
        if self.cluster_frame_index:
            cluster_frame = self.frames[self.cluster_frame_index]
            for widget in cluster_frame.grid_slaves():
                if isinstance(widget, tkinter.ttk.Frame):
                    widget.grid_forget()

        if self.recluster_frame_index:
            recluster_frame = self.frames[self.recluster_frame_index]
            for widget in recluster_frame.grid_slaves():
                if isinstance(widget, tkinter.ttk.Frame):
                    widget.grid_forget()

        self.refresh_menus()

        # Delete rest informations about appended sections
        for section_name in self.values.keys():
            if section_name.startswith("clustering") or section_name.startswith("clustering"):
                if section_name != "clustering" and section_name != "reclustering":
                    del self.values[section_name]
                    del self.hiding_frames[section_name]

                    # Remove from MENUS list
                    for i, item in enumerate(defaults.MENUS):
                        if item.startswith(section_name):
                            del defaults.MENUS[i]

    def get_active_frame_name(self, section_name):
        """
        Return currently gridded HidingFrame in section_name.

        :param section_name: Name of section from which active HidingFrame will be returned.
        :return: Name of chosen option menu value.
        """
        for method, frame in self.hiding_frames[section_name].items():
            if any(frame.grid_info()):
                return method

    def save_config(self, config_filename, required_checking=True):
        """
        Save all values to config file.

        :param config_filename: Name of file where options will be saved.
        :param required_checking: If False checking required fields and saving dialog is skipped. Needed to save options before changing level.
        """
        # Unhighlight all frames
        for section_name in self.required_entries:
            for entry_name in self.required_entries[section_name]:
                self.required_entries[section_name][entry_name].unhighlight()

        if required_checking:
            # Check if all field all filled before saving
            tabs_id = self.notebook.tabs()
            first_found = [None, None, None]  # Keeps info of first found unfilled entry
            for i, section_name in enumerate(self.values):
                for entry_name in self.values[section_name]:
                    if isinstance(self.values[section_name][entry_name], dict):
                        for entry_name_optionmenu, entry_value_optionmenu in self.values[section_name][entry_name].items():
                            if not entry_value_optionmenu.get() and entry_name_optionmenu in self.required_entries[section_name]:
                                section_full_name = defaults.get_default_section(section_name).name
                                entry_full_name = defaults.get_default_entry(section_name, entry_value_optionmenu).name[:-2]

                                if not first_found[0]:
                                    first_found[0] = section_full_name
                                    first_found[1] = entry_full_name
                                    first_found[2] = tabs_id[i]

                                self.required_entries[section_name][entry_name].highlight()
                    else:
                        if not self.values[section_name][entry_name].get() and entry_name in self.required_entries[section_name]:
                            section_full_name = defaults.get_default_section(section_name).name
                            entry_full_name = defaults.get_default_entry(section_name, entry_name).name[:-2]

                            if not first_found[0]:
                                first_found[0] = section_full_name
                                first_found[1] = entry_full_name
                                first_found[2] = tabs_id[i]

                            self.required_entries[section_name][entry_name].highlight()

            if first_found[2]:  # 0 and 1 indexes is not checked because it can be just an empty string
                tkinter.messagebox.showerror("Unfilled field",
                                       "Field \"{}\" in \"{}\" must be specified.".format(entry_full_name,
                                                                                          section_full_name))
                self.notebook.select(first_found[2])
                return False

        # Create config file backup
        shutil.copy(config_filename, config_filename + ".bak")
        try:
            with open(config_filename, "w+") as config_file:
                config = ConfigParser()

                for section in defaults.DEFAULTS:
                    if section.additional:
                        continue

                    if section.config_name == "clustering" and not self.re_clustering_exists(section.config_name):
                        continue

                    if section.config_name == "reclustering" and not self.re_clustering_exists(section.config_name):
                        continue

                    if not config.has_section(section.config_name):
                        config.add_section(section.config_name)

                    for entry_section, entry in section.iter_entries():
                        value = None

                        # Try to get value from entry
                        try:
                            if entry.optionmenu_value:
                                value = self.values[entry_section.config_name][entry.optionmenu_value][entry.config_name].get()
                            else:
                                value = self.values[entry_section.config_name][entry.config_name].get()
                        except KeyError:
                            value = entry.default_value

                        if value == "" or (not isinstance(value, bool) and value == 0):
                            continue

                        if entry.optionmenu_value:
                            # Save options only from chosen options in optionmenu
                            if entry.optionmenu_value == self.get_active_frame_name(section.config_name):
                                config.set(section.config_name, entry.config_name, str(value))
                        else:
                            config.set(section.config_name, entry.config_name, str(value))

                # Add created clustering and reclustering sections
                clustering_section_names = [key for key in list(self.values.keys()) if
                                            key.startswith("clustering") and not key == "clustering"]
                reclustering_section_names = [key for key in list(self.values.keys()) if
                                              key.startswith("reclustering") and not key == "reclustering"]

                for section_name in clustering_section_names + reclustering_section_names:
                    section_config_name = self.values[section_name]["name"].get()
                    config.add_section(section_config_name)

                    for entry_section, entry in defaults.get_default_section(section_name).iter_entries():
                        value = None

                        # Try to get value from entry
                        try:
                            if entry.optionmenu_value:
                                value = self.values[section_name][entry.optionmenu_value][entry.config_name].get()
                            else:
                                value = self.values[section_name][entry.config_name].get()
                        except KeyError:
                            value = entry.default_value

                        if value == "" or value == 0:
                            continue

                        if entry.optionmenu_value:
                            # Save options only from chosen options in optionmenu
                            if entry.optionmenu_value == self.get_active_frame_name(section_name):
                                config.set(section_config_name, entry.config_name, str(value))
                        else:
                            config.set(section_config_name, entry.config_name, str(value))

                config.write(config_file)
                self.values_hash = self.get_values_hash()
        except Exception as e:
            # In case of error restore config from backup file
            tkinter.messagebox.showinfo("Exeption",
                                  "There was error during saving configuration. See console for more information.\n{}".format(
                                      str(e)))
            print((traceback.print_exc()))
            shutil.copy(config_filename + ".bak", config_filename)

        os.remove(config_filename + ".bak")

        if required_checking:
            tkinter.messagebox.showinfo("Saved", "Saving complete")

        return True

    def open_docs(self):
        webbrowser.open("http://www.aquaduct.pl/apidocs/valve/valve_config.html")

    def create_new_config_file_dialog(self):
        """
        Show dialog with create or choose existing file options.
        """
        window = tk.Toplevel(self.parent, padx=20)
        window.title("New file")

        directory_name = tk.StringVar()
        file_name = tk.StringVar()

        def select_dir():
            directory_name.set(askdirectory())
            window.lift()  # Move to foreground

        tkinter.ttk.Label(window, text="Directory:").grid(row=0, column=0)
        tkinter.ttk.Entry(window, textvariable=directory_name, state="readonly").grid(row=0, column=1)
        tkinter.ttk.Button(window, text="Choose directory", command=select_dir, style="File.TButton").grid(row=0, column=2)

        tkinter.ttk.Label(window, text="File name:").grid(row=1, column=0)
        tkinter.ttk.Entry(window, textvariable=file_name).grid(row=1, column=1)

        def create(dir_var, filename_var):
            path = dir_var.get() + os.path.sep + filename_var.get()
            with open(path, "w+"):
                self.config_filename.set(path)
                self.parent.title("{} - {}".format(self.title, path))

            window.destroy()

        create_callback = utils.CallbackWrapper(create, directory_name, file_name)
        tkinter.ttk.Button(window, text="Create", command=create_callback).grid(row=3, column=0, columnspan=3)

    def about(self):
        window = tk.Toplevel(self.parent)
        window.title("About")

        logo = tk.PhotoImage(file=get_img("logo.gif"))

        logo_label = tkinter.ttk.Label(window, image=logo)
        logo_label.image = logo
        logo_label.pack()

        content = """Aqua-Duct {}

ValveConfigurator

Copyright Tunneling Group \xa9 2018""".format(aquaduct.version_nice())

        tkinter.ttk.Label(window, text=content, justify=tk.CENTER).pack()
        tkinter.ttk.Label(window, text=aquaduct.__author__, justify=tk.CENTER).pack()

        hyperlink = tk.Label(window, text="http://www.tunnelinggroup.pl/", foreground="blue", justify=tk.CENTER)
        hyperlink.pack()
        hyperlink_cb = utils.CallbackWrapper(webbrowser.open_new, "http://www.tunnelinggroup.pl/")
        hyperlink.bind("<Button-1>", hyperlink_cb)

    def valve_run_dialog(self):
        """ Open dialog with options to run Valve """
        if self.config_filename.get() == "":
            tkinter.messagebox.showinfo("Config file", "You have to save or open existing config.")
            return
        elif self.values_hash != self.get_values_hash():
            result = tkinter.messagebox.askyesno("Unsaved config", "There are unsaved changed. Are you sure to continue?")

            if not result:
                return

        window = tk.Toplevel(self.parent, padx=20, pady=15)
        window.title("Run Valve")

        scrolled_frame = utils.VerticalScrolledFrame(window)
        scrolled_frame.pack(expand=1, fill=tk.BOTH)
        frame = scrolled_frame.interior
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        command = tk.StringVar()
        args = {}

        def update_command(*a):
            cmd = "valve_run "
            for arg, var in args.items():
                if isinstance(var.get(), bool):
                    if var.get():
                        cmd += arg + " "
                elif var.get():
                    cmd += "{} {} ".format(arg, var.get())

            command.set(cmd)

        for row, entry in enumerate(defaults.VALVE_DEFAULTS.entries):
            var = utils.entry_factory(frame, row, entry.name, entry.default_values, entry.help_text)
            args[entry.config_name] = var
            var.input_var.trace("w", update_command)

        utils.Text(frame, textvariable=command, width=60, height=2).grid(row=999,
                                                                         column=0,
                                                                         columnspan=2,
                                                                         sticky="ew",
                                                                         pady=10)

        # Update field with opened config file
        args["-c"].set(self.config_filename.get())

        run_cb = utils.CallbackWrapper(os.system, command.get() + " & disown")
        run_button = tkinter.ttk.Button(frame, text="Run")
        run_button.bind("<Button-1>", run_cb)
        run_button.grid(row=1000, column=0, columnspan=2)

    def pond_run_dialog(self):
        """ Open dialog with options to run Valve """
        if self.config_filename.get() == "":
            tkinter.messagebox.showinfo("Config file", "You have to save or open existing config.")
            return
        elif self.values_hash != self.get_values_hash():
            result = tkinter.messagebox.askyesno("Unsaved config", "There are unsaved changed. Are you sure to continue?")

            if not result:
                return

        window = tk.Toplevel(self.parent, padx=20, pady=15)
        window.title("Run Pond")

        scrolled_frame = utils.VerticalScrolledFrame(window)
        scrolled_frame.pack(expand=1, fill=tk.BOTH)
        frame = scrolled_frame.interior
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        command = tk.StringVar()
        args = {}

        def update_command(*a):
            cmd = "pond_run "
            for arg, var in args.items():
                if isinstance(var.get(), bool):
                    if var.get():
                        cmd += arg + " "
                elif var.get():
                    cmd += "{} {} ".format(arg, var.get())

            command.set(cmd)

        for row, entry in enumerate(defaults.POND_DEFAULTS.entries):
            var = utils.entry_factory(frame, row, entry.name, entry.default_values, entry.help_text)
            args[entry.config_name] = var
            var.input_var.trace("w", update_command)

        utils.Text(frame, textvariable=command, width=60, height=2).grid(row=999,
                                                                         column=0,
                                                                         columnspan=2,
                                                                         sticky="ew",
                                                                         pady=10)

        # Update field with opened config file
        args["-c"].set(self.config_filename.get())

        run_cb = utils.CallbackWrapper(os.system, command.get())
        run_button = tkinter.ttk.Button(frame, text="Run")
        run_button.bind("<Button-1>", run_cb)
        run_button.grid(row=1000, column=0, columnspan=2)

    def on_window_close(self):
        if self.get_values_hash() != self.values_hash:
            result = tkinter.messagebox.askyesno("Unsaved changes", "You have unsaved changes. Are you want to quit?")
            if result:
                self.parent.destroy()
        else:
            self.parent.destroy()


if __name__ == "__main__":
    root = tk.Tk()
    root.configure(background="white")
    root.resizable(1, 1)

    aq_icon = tk.PhotoImage(file=get_img("icon.gif"))
    root.tk.call('wm', 'iconphoto', root._w, aq_icon)

    s = tkinter.ttk.Style()

    s.theme_use("clam")

    ###
    # Global styles
    ###
    font_name = "TkDefaultFont"
    font_size = 8
    s.configure("TLabel", padding=5, background="#dcdad5", font=(font_name, font_size))
    s.configure("TButton", padding=5, background="#dcdad5", font=(font_name, font_size))
    s.configure("TMenubutton", padding=5, background="#dcdad5", font=(font_name, font_size))

    # Top label in every self.notebook tab
    s.configure("Title.TLabel", font=("TkDefaultFont", 16))

    ###
    # self.notebook style
    ###
    s.configure("TNotebook.Tab", width=20, anchor=tk.CENTER, font=("TkDefaultFont", 8))

    # s.configure("Tself.notebook", tabposition="wn")
    # s.map("Tself.notebook.Tab",
    # width=[("selected", 20)],
    # padding=[("selected", [0, 7])]
    # )

    ###
    # For recursive clustering frames
    ###
    s.configure("I.TFrame", background="#E6E4DF", relief=tk.SUNKEN)
    s.configure("I.TFrame.Label", background="#E6E4DF")
    s.configure("ITitle.TLabel", background="#E6E4DF", font=("TkDefaultFont", 14))

    ###
    # For hiding frames
    ###
    s.configure("HF.TFrame", background="#F0EEE9", relief=tk.RAISED)
    s.configure("HF.TFrame.Label", background="#F0EEE9", font=("TkDefaultFont", 9, "bold"))

    ###
    # Other
    ###
    s.configure("File.TButton", padding=0, font=("TkDefaultFont", 8))  # Loading file button
    s.configure("Configured.TLabel", padding=0, font=("TkDefaultFont", 12))  # Loading file button
    s.configure("R.TButton", background="red3", bordercolor="red3")  # Reset button

    ###
    # Highlighted
    ###
    s.configure("Highlighted.TFrame", background="red")

    app = ValveConfigApp(root)
    root.mainloop()
