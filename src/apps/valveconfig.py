#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2018 Michał Banas
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
import hashlib
import os
import shutil
import tempfile
import tkMessageBox
import ttk
from ConfigParser import ConfigParser, NoOptionError, NoSectionError
from collections import OrderedDict, defaultdict
from tkFileDialog import askopenfile, askdirectory

import aquaduct.apps.valveconfig.defaults as defaults
import aquaduct.apps.valveconfig.utils as utils
from aquaduct import version_nice


class ValveConfigApp(object):
    def __init__(self, parent):
        """ Valve Configurator App """
        self.size = (600, 650)

        self.parent = parent

        self.values = OrderedDict()

        # Used to hide that frame, after all entries are shown
        self.init_frame = ttk.Frame(self.parent)
        self.notebook = ttk.Notebook(self.parent)

        self.level = None

        self.hiding_frames = {}

        self.cluster_frame_index = None
        # Row of the last element in clustering tab
        self.cluster_row = 0

        self.recluster_frame_index = None
        # Row of the last element in reclustering tab
        self.recluster_row = 0

        # Used for identifying scroll frame by interior frame for removing it from notebook
        self.scrolled_frames = {}

        self.frames = []

        # Dict with input frames of required entries. Used to change color if is unfilled
        self.required_entries = defaultdict(dict)

        self.config_filename = tk.StringVar()

        self.values_hash = self.get_values_hash()

        parent.title("Valve Configurator")
        parent.geometry("600x650")

        self.parent.protocol("WM_DELETE_WINDOW", self.on_window_close)

        self.init_gui()

    def init_gui(self):
        """ Prepare initial frame """
        # Logo
        logo = tk.PhotoImage(data=utils.LOGO_ENCODED)

        logo_label = ttk.Label(self.parent, image=logo, padding=-2)
        logo_label.image = logo
        logo_label.pack(padx=20, pady=20)

        # Used to auto positioning depending on length of version string
        version = "ver. " + version_nice()
        ttk.Label(logo_label, text=version, background="white").place(relx=1, rely=1, x=-len(version) * 7, y=-20)

        # Frame with loading configuration file
        load_frame = ttk.LabelFrame(self.init_frame, text="Configuration file")
        load_frame.columnconfigure(0, weight=1)
        load_frame.columnconfigure(1, weight=1)

        ttk.Button(load_frame, text="Load config", command=self.open_config_file).grid(sticky="E", row=0, column=0)
        ttk.Button(load_frame, text="New file", command=self.create_new_config_file).grid(sticky="W",
                                                                                          row=0,
                                                                                          column=1,
                                                                                          padx=5,
                                                                                          pady=2)
        ttk.Entry(load_frame, textvariable=self.config_filename, state="readonly").grid(sticky="EW",
                                                                                        row=1,
                                                                                        column=0,
                                                                                        columnspan=2,
                                                                                        padx=5,
                                                                                        pady=2)

        # Level selection
        level_frame = ttk.LabelFrame(self.init_frame, text="Configuration level")

        self.level = tk.StringVar()

        # List of sorted levels from easiest to hardest
        self.levels = []
        for level_name, level in reversed(sorted(defaults.LEVELS.iteritems(), key=lambda (k, v): (v, k))):
            self.levels.append(level_name)

        ttk.OptionMenu(level_frame, self.level, self.levels[0], *self.levels).pack(pady=20)

        load_frame.pack(fill=tk.X, padx=100, pady=20)
        level_frame.pack(fill=tk.X, padx=100, pady=20)

        def show_new_gui():
            self.prepare_section_frames()
            self.create_interface()

        ttk.Button(self.init_frame, text="Forward", command=show_new_gui).pack(pady=50)

        self.init_frame.pack(expand=1, fill="both")

    def prepare_section_frames(self):
        """ Parse DEFAULTS and depends on it create Notebook tabs and entries """
        for section in defaults.DEFAULTS:
            section_name = section.config_name
            section_name_long = section.name
            entries = section.entries

            if defaults.LEVELS[self.level.get()] > section.level:
                continue

            if section_name == "clusterization":
                self.cluster_frame_index = len(self.frames)

            if section_name == "reclusterization":
                self.recluster_frame_index = len(self.frames)

            scroll_frame = utils.VerticalScrolledFrame(self.notebook)

            frame = scroll_frame.interior
            frame.columnconfigure(0, weight=1)
            frame.columnconfigure(1, weight=1)

            self.frames.append(frame)
            self.scrolled_frames[frame] = scroll_frame

            self.notebook.add(scroll_frame, text=section_name_long)

            ttk.Label(frame,
                      text=section_name_long,
                      style="Title.TLabel",
                      background=utils.get_widget_bg(self.parent)).grid(row=0, column=0, columnspan=2)
            ttk.Separator(frame, orient=tk.HORIZONTAL).grid(sticky="EW", row=1, column=0, columnspan=3)

            self.entry_filler(frame, section_name, entries, 2)

        # Set row number to row of hidden frame of default clustering section. Used to append new frames
        # TODO: It may not work when other widgets are rendered first
        self.cluster_row = next(self.hiding_frames["clusterization"].itervalues()).row + 1

        if "reclusterization" in self.hiding_frames:
            self.recluster_row = next(self.hiding_frames["reclusterization"].itervalues()).row + 1

        for section_name in self.get_recursive_clustering_sections("clusterization"):
            self.append_entries(section_name)

            # Increment max_level
            self.add_max_level(1)

        for section_name in self.get_recursive_clustering_sections("reclusterization"):
            self.append_entries(section_name)

        cluster_add_button = ttk.Button(self.frames[self.cluster_frame_index], text="Add clustering section")
        # Setting row=1000 let skip calculating position of button each time new section is added
        cluster_add_button.grid(row=1000, column=0, columnspan=2, pady=20)

        cluster_add_section_callback = utils.CallbackWrapper(self.callback_add_section, "clusterization")
        cluster_add_button.bind("<Button-1>", cluster_add_section_callback)

        if "reclusterization" in self.hiding_frames:
            recluster_add_button = ttk.Button(self.frames[self.recluster_frame_index], text="Add reclustering section")
            # Setting row=1000 let skip calculating position of button each time new section is added
            recluster_add_button.grid(row=1000, column=0, columnspan=2, pady=20)

            recluster_add_section_callback = utils.CallbackWrapper(self.callback_add_section, "reclusterization")
            recluster_add_button.bind("<Button-1>", recluster_add_section_callback)

        # Refresh all hiding frames
        self.refresh_menus()

        self.init_frame.pack_forget()

        self.notebook.pack(expand=1, fill="both")

        # After preparing all section load values from config
        if self.config_filename.get() != "":
            self.load_config_values(self.config_filename.get())
        else:
            self.values_hash = self.get_values_hash()

    def create_interface(self):
        """ Creates bottom interface """
        bottom_frame = ttk.Frame(self.parent)
        bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

        save_button = ttk.Button(bottom_frame, text="Save")

        def save_callback(*args):
            if self.config_filename.get() == "":
                self.open_config_file()

                # If creating new file was canceled
                if self.config_filename.get() == "":
                    return

            self.save_config(self.config_filename.get())

        save_button.bind("<Button-1>", save_callback)
        save_button.pack(side=tk.LEFT, pady=5, padx=40)

        ttk.OptionMenu(bottom_frame, self.level, self.levels[0], *self.levels).pack(side=tk.LEFT)
        self.level.trace("w", self.recreate_gui)

        reset_button = ttk.Button(bottom_frame, text="Reset to default", style="R.TButton")

        reset_callback = utils.CallbackWrapper(self.load_defaults, reset=True)
        reset_button.bind("<Button-1>", reset_callback)
        reset_button.pack(side=tk.RIGHT, pady=5, padx=40)

    def callback_add_section(self, section_name):
        """
        Callback for button to create new recursive clusterization in clusterization/reclusterization frame

        :param section_name: Config section name where new frames will be appended, allowed are "clusterization" or "reclusterization"
        """
        # Find free name for new section
        index = 0
        while section_name + str(index) in self.values:
            index += 1

        self.append_entries(section_name + str(index))

        if section_name.startswith("clusterization"):
            self.add_max_level(1)

        self.refresh_menus()

    def callback_remove_section(self, section_name, frame):
        """
        Callback for button to remove existing recursive clusterization in clusterization/reclusterization frame

        :param section_name: Config section name to remove, allowed are "clusterization" or "reclusterization"
        :param frame: frame Which contains options related to that section
        """
        del self.values[section_name]
        del self.hiding_frames[section_name]

        # Remove from MENUS list
        for i, item in enumerate(defaults.MENUS):
            if item.startswith(section_name):
                del defaults.MENUS[i]

        if section_name.startswith("clusterization"):
            self.add_max_level(-1)

        frame.grid_forget()

    def get_recursive_clustering_sections(self, section_name):
        """
        Finds recursively all section used in recursive_clusterization options

        :param section_name: Config section name from which fetching will start.
        :return Section names.
        :rtype: list
        """
        if section_name != "clusterization" and section_name != "reclusterization":
            raise RuntimeError(
                "There is no possibility to get recursive clusterization from {} section".format(section_name))

        try:
            with open(self.config_filename.get(), "r") as config_file:
                config = ConfigParser()
                config.readfp(config_file)

                clustering_sections = []
        except IOError:  # File is empty
            return []

        try:
            clustering_section = config.get(section_name, "recursive_clusterization")
            while clustering_section not in clustering_sections:
                clustering_sections.append(clustering_section)
                clustering_section = config.get(clustering_section, "recursive_clusterization")
        except (NoOptionError, NoSectionError):
            pass

        return clustering_sections

    def append_entries(self, section_name):
        """
        Append new frame with new clusterization or reclusterization options.

        :param section_name: Config section name where new frames will be appended, allowed are "clusterization" or "reclusterization".
        """
        if section_name.startswith("clusterization"):
            default_section_name = "clusterization"
            frame = self.frames[self.cluster_frame_index]
            row = self.cluster_row
        elif section_name.startswith("reclusterization"):
            default_section_name = "reclusterization"
            frame = self.frames[self.recluster_frame_index]
            row = self.recluster_row
        else:
            raise RuntimeError(
                "Appending entries to sections other than clusterization or reclusteriation is not allowed")

        # Set option menu from recursive clustering to control hiding frames
        defaults.MENUS.append("{}:method".format(section_name))

        default_entries = defaults.get_default_section(default_section_name)

        inner_frame = ttk.Frame(frame, style="I.TFrame")
        inner_frame.columnconfigure(0, weight=1)
        inner_frame.columnconfigure(1, weight=1)
        inner_frame.grid(sticky="EW", row=row, column=0, columnspan=2, padx=60, pady=10)

        ttk.Label(inner_frame, text="Recursive clusterization", style="ITitle.TLabel").grid(
            row=0,
            column=0,
            pady=5,
            columnspan=2)
        ttk.Separator(inner_frame, orient=tk.HORIZONTAL).grid(sticky="EW", row=1, column=0, columnspan=2)

        self.entry_filler(inner_frame, section_name, default_entries.entries, 2)

        remove_button = ttk.Button(inner_frame, text="Remove")
        remove_button.grid(row=999, column=0, columnspan=2, pady=5)

        remove_section_callback = utils.CallbackWrapper(self.callback_remove_section, section_name, inner_frame)
        remove_button.bind("<Button-1>", remove_section_callback)

        if default_section_name == "clusterization":
            self.cluster_row += 1
        else:
            self.recluster_row += 1

    def entry_filler(self, parent, section_name, entries, row):
        """
        Creates entries depending on DEFAULTS values.

        When no entry added it hide notebook tab.

        Disables name option for clusterization, reclusterization and derivatives section.

        :param parent: Parent widget.
        :param section_name: Config section name.
        :param entries: List of entries.
        :param row: Row number at which first entry will be shown.
        """
        # Keeps label frames with actual row
        group_frames = {}

        entries_appended = 0
        for entry in entries:
            # If its primary reclustering or clustering section disable changing name of that section
            # Additionaly it disable inlets_clusterization:max_level if its on different level than chose
            state = tk.NORMAL
            if section_name == "clusterization" and entry.config_name == "name":
                state = tk.DISABLED
            elif section_name == "reclusterization" and entry.config_name == "name":
                state = tk.DISABLED
            elif section_name == "inlets_clusterization" and entry.config_name == "max_level":
                if defaults.LEVELS[self.level.get()] != entry.level:
                    state = tk.DISABLED

            # Skip entries that dont match level except inlets_clusterization:max_level
            if defaults.LEVELS[self.level.get()] > entry.level:
                if not (section_name == "inlets_clusterization" and entry.config_name == "max_level"):
                    continue

            entries_appended += 1
            if entry.optionmenu_value:
                if section_name not in self.hiding_frames:
                    self.hiding_frames[section_name] = {}

                if entry.optionmenu_value not in self.hiding_frames[section_name]:
                    if len(self.hiding_frames[section_name]) != 0:
                        hiding_frame_row = next(self.hiding_frames[section_name].itervalues()).row
                    else:
                        hiding_frame_row = row

                    self.hiding_frames[section_name][entry.optionmenu_value] = utils.HidingFrame(parent,
                                                                                                 hiding_frame_row,
                                                                                                 text=entry.optionmenu_value.capitalize() + " options",
                                                                                                 style="HF.TFrame")

                hiding_frame = self.hiding_frames[section_name][entry.optionmenu_value]

                if section_name not in self.values:
                    self.values[section_name] = {}

                self.values[section_name][entry.config_name] = utils.entry_factory(hiding_frame,
                                                                                   hiding_frame.inner_row,
                                                                                   entry.name,
                                                                                   entry.default_values,
                                                                                   entry.help_text,
                                                                                   state,
                                                                                   info_text=entry.info_text,
                                                                                   warning_text=entry.warning_text)

                hiding_frame.inner_row += 1

                if entry.required:
                    self.required_entries[section_name][entry.config_name] = self.values[section_name][
                        entry.config_name]
            else:
                if entry.group_label:
                    if entry.group_label not in group_frames:
                        label_frame = ttk.LabelFrame(parent, text=entry.group_label)
                        label_frame.grid_columnconfigure(0, weight=1)
                        label_frame.grid_columnconfigure(1, weight=1)
                        label_frame.grid(row=row, column=0, columnspan=2, pady=15, ipadx=30)

                        group_frames[entry.group_label] = [label_frame, 0]

                    entry_parent = group_frames[entry.group_label][0]
                    entry_row = group_frames[entry.group_label][1]

                    group_frames[entry.group_label][1] += 1
                else:
                    entry_parent = parent
                    entry_row = row

                entry_widget = utils.entry_factory(entry_parent,
                                                   entry_row,
                                                   entry.name,
                                                   entry.default_values,
                                                   entry.help_text,
                                                   state,
                                                   info_text=entry.info_text,
                                                   warning_text=entry.warning_text)

                if (section_name.startswith("clusterization") or
                    section_name.startswith("reclusterization")) and entry.config_name == "name":
                    entry_widget.set(section_name)

                if defaults.is_menu(section_name, entry.config_name):
                    opt_menu_callback = utils.CallbackWrapper(self.option_menu_changed, section_name, entry_widget)
                    entry_widget.input_var.trace("w", opt_menu_callback)

                if section_name not in self.values:
                    self.values[section_name] = {}

                self.values[section_name][entry.config_name] = entry_widget

                row += 1

                if entry.required:
                    self.required_entries[section_name][entry.config_name] = entry_widget

        # Hide tab if unused
        if entries_appended == 0:
            self.notebook.forget(self.scrolled_frames[parent])

    def refresh_menus(self):
        """
        Sets visible only hiding frames that are chosen by option menu.
        """
        for menu in defaults.MENUS:
            section_name, entry_name = menu.split(":")

            if defaults.LEVELS[self.level.get()] <= defaults.get_default_section(section_name).level:
                self.option_menu_changed(section_name, self.values[section_name][entry_name])

    def option_menu_changed(self, section, entry):
        """
        Callback for option menu that control hiding frames

        :param section: Section name where hiding frame is located.
        :param entry: Option menu widget, which control hiding frames in section.
        """
        for method, hidden_frame in self.hiding_frames[section].iteritems():
            if method == entry.get():
                hidden_frame.show()
            else:
                hidden_frame.grid_forget()

    def add_max_level(self, value):
        """
        Increments or decrements inlets_clusterization:max_level by value.

        :param value: Number which will be added to max_level.
        """
        current_value = self.values["inlets_clusterization"]["max_level"].get()
        self.values["inlets_clusterization"]["max_level"].set(current_value + value)

    def open_config_file(self):
        """ Show dialog to choose config file and save its name """
        filetypes = (("text files", "*.txt"), ("config files", "*.cfg"), ("all files", "*.*"))
        try:
            with askopenfile("r", filetypes=filetypes) as config_file:
                self.config_filename.set(config_file.name)
        except AttributeError:
            pass

    def recreate_gui(self, *args):
        """
        Recreate GUI.

        Saves current values to temporary file.

        :param args: Event informations.
        """
        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmp_filename = tmpfile.name

        # Save values to temporary file.
        self.save_config(tmpfile.name, required_checking=False)

        for tab_id in self.notebook.tabs():
            self.notebook.forget(tab_id)

        self.prepare_section_frames()
        self.load_defaults()  # Reset values
        self.load_config_values(tmp_filename)

        os.remove(tmp_filename)

    def get_values_hash(self):
        """ Compute hash from values. """
        hash = hashlib.md5()
        for section_name in self.values:
            for value in self.values[section_name].itervalues():
                hash.update(str(value.get()))

        return hash.digest()

    def load_config_values(self, config_filename):
        """
        Load config values to values dictionary from file.

        Additionally it ask user if inlets_clusterization:max_level in config
        is different from number of found recursive_clusterization sections.

        :param config_filename: Name of file with configuration.
        """
        with open(config_filename, "r") as config_file:
            config = ConfigParser()
            config.readfp(config_file)

            for section_name in self.values:
                for entry_name in self.values[section_name]:

                    try:
                        config_value = config.get(section_name, entry_name)

                        #  Inform user that config max_level differ from number of found recursive clusterization
                        result = True
                        if section_name == "inlets_clusterization" and entry_name == "max_level":
                            if self.values["inlets_clusterization"]["max_level"] != config_value:
                                result = tkMessageBox.askyesno("Max level",
                                                               "Config max_level value differs from number of clusterization section found.\n"
                                                               "Do you want to keep value from configure file?")

                        if not result:
                            continue

                        self.values[section_name][entry_name].set(config_value)
                    except NoOptionError:
                        continue
                    except NoSectionError:
                        break

        self.values_hash = self.get_values_hash()

    def load_defaults(self, reset=False, *args):
        """
        Load defaults values to values dictionary.

        :param reset: If set to True confirmation will be neeeded to reset.
        """
        if reset:
            result = tkMessageBox.askyesno("Reset values", "Are you sure to reset all values?")

            if not result:
                return

        for section in defaults.DEFAULTS:
            if defaults.LEVELS[self.level.get()] > section.level:
                continue

            for entry in section.entries:
                # Skip entries that dont match level except inlets_clusterization:max_level
                if defaults.LEVELS[self.level.get()] > entry.level:
                    if not (section.config_name == "inlets_clusterization" and entry.config_name == "max_level"):
                        continue

                self.values[section.config_name][entry.config_name].set(entry.default_value)

        # Delete appended clustering and reclustering sections
        if self.cluster_frame_index:
            cluster_frame = self.frames[self.cluster_frame_index]
            for widget in cluster_frame.grid_slaves():
                if isinstance(widget, ttk.Frame):
                    widget.grid_forget()

        if self.recluster_frame_index:
            recluster_frame = self.frames[self.recluster_frame_index]
            for widget in recluster_frame.grid_slaves():
                if isinstance(widget, ttk.Frame):
                    widget.grid_forget()

        # Due to that earlier loop hide hiding frames too its necessary to refresh them
        self.refresh_menus()

        # Delete rest informations about appended sections
        for section_name in self.values.iterkeys():
            if section_name.startswith("clusterization") or section_name.startswith("clusterization"):
                if section_name != "clusterization" and section_name != "reclusterization":
                    del self.values[section_name]
                    del self.hiding_frames[section_name]

                    # Remove from MENUS list
                    for i, item in enumerate(defaults.MENUS):
                        if item.startswith(section_name):
                            del defaults.MENUS[i]

    def create_new_config_file(self):
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

        ttk.Label(window, text="Directory:").grid(row=0, column=0)
        ttk.Entry(window, textvariable=directory_name, state="readonly").grid(row=0, column=1)
        ttk.Button(window, text="Choose directory", command=select_dir, style="File.TButton").grid(row=0, column=2)

        ttk.Label(window, text="File name:").grid(row=1, column=0)
        ttk.Entry(window, textvariable=file_name).grid(row=1, column=1)

        def create(dir_var, filename_var):
            path = dir_var.get() + os.path.sep + filename_var.get()
            with open(path, "w+"):
                self.config_filename.set(path)

        create_callback = utils.CallbackWrapper(create, directory_name, file_name)
        ttk.Button(window, text="Create", command=create_callback).grid(row=3, column=0, columnspan=3)

    def get_active_frame_name(self, section_name):
        """
        Return currently visible frame.

        :param section_name: Section name in which hiding frames will be considered.
        :return: Name of chosen option menu value.
        """
        for method, frame in self.hiding_frames[section_name].iteritems():
            if any(frame.grid_info()):
                return method

    def save_config(self, config_filename, required_checking=True):
        """
        Save all values to config file.

        :param config_filename: Name of file where options will be saved.
        :param required_checking: If False checking required fields and saving dialog is skipped.
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
                    if not self.values[section_name][entry_name].get() and entry_name in self.required_entries[
                        section_name]:
                        section_full_name = defaults.get_default_section(section_name).name
                        entry_full_name = defaults.get_default_entry(section_name, entry_name).name[:-2]

                        if not first_found[0]:
                            first_found[0] = section_full_name
                            first_found[1] = entry_full_name
                            first_found[2] = tabs_id[i]

                        self.required_entries[section_name][entry_name].highlight()

            if first_found[2]:  # 0 and 1 indexes is not checked because it can be just an empty string
                tkMessageBox.showerror("Unfilled field",
                                       "Field \"{}\" in \"{}\" must be specified.".format(entry_full_name,
                                                                                          section_full_name))
                self.notebook.select(first_found[2])
                return

        # Create config file backup
        shutil.copy(config_filename, config_filename + ".bak")
        try:
            with open(config_filename, "w+") as config_file:
                config = ConfigParser()
                config.readfp(config_file)

                for section_name in self.values:
                    if not config.has_section(section_name):
                        config.add_section(section_name)

                    for option_name, value in self.values[section_name].iteritems():
                        default_value = defaults.get_default_entry(section_name, option_name).default_value

                        # Skip empty strings and entries which are not visible
                        if value.get() != "" and value.get() != default_value:
                            if not defaults.get_default_entry(section_name, option_name).optionmenu_value:
                                config.set(section_name, option_name, value.get())
                            # This take care of saving options only from visible frame
                            else:
                                if defaults.get_default_entry(section_name,
                                                              option_name).optionmenu_value == self.get_active_frame_name(
                                    section_name):
                                    config.set(section_name, option_name, value.get())

                config.write(config_file)
                self.values_hash = self.get_values_hash()
        except Exception as e:
            # In case of error restore config from backup file
            tkMessageBox.showinfo("Exeption", "There was error during saving configuration.\n{}".format(str(e)))
            shutil.copy(config_filename + ".bak", config_filename)

        if required_checking:
            tkMessageBox.showinfo("Saved", "Saving complete")

    def on_window_close(self):
        if self.get_values_hash() != self.values_hash:
            result = tkMessageBox.askyesno("Unsaved changes", "You have unsaved changes. Are you want to quit?")
            if result:
                self.parent.destroy()
        else:
            self.parent.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    root.configure(background="white")
    root.resizable(1, 1)

    s = ttk.Style()

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
    # For recursive clusterization frames
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
