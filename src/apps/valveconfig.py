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
import tkMessageBox
import ttk
from ConfigParser import ConfigParser, NoOptionError, NoSectionError
from tkFileDialog import askopenfile

import aquaduct.apps.valveconfig.utils as utils


class ValveConfigApp(object):
    def __init__(self, parent):
        """ Valve Configurator App """
        self.size = (600, 650)

        self.parent = parent

        self.values = {}

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

        self.config_filename = tk.StringVar()

        parent.title("Valve Configurator")
        parent.geometry("600x650")

        self.init_gui()

    def init_gui(self):
        """ Prepare initial frame """
        # Logo
        logo = tk.PhotoImage(data=utils.LOGO_ENCODED)

        logo_label = ttk.Label(self.parent, image=logo, padding=-2)
        logo_label.image = logo
        logo_label.pack(padx=20, pady=20)

        # Frame with loading configuration file
        load_frame = ttk.LabelFrame(self.init_frame, text="Configuration file")

        ttk.Button(load_frame, text="Load config", command=self.open_config_file).pack(anchor=tk.CENTER,
                                                                                       side=tk.LEFT,
                                                                                       pady=20, padx=20)
        ttk.Entry(load_frame, textvariable=self.config_filename, state="readonly").pack(anchor=tk.CENTER,
                                                                                        expand=1,
                                                                                        fill=tk.X,
                                                                                        side=tk.RIGHT,
                                                                                        pady=20,
                                                                                        padx=5)

        # Level selection
        level_frame = ttk.LabelFrame(self.init_frame, text="Configuration level")

        self.level = tk.StringVar()
        ttk.OptionMenu(level_frame, self.level, "Easy", "Easy", "Normal").pack(pady=20)

        load_frame.pack(fill=tk.X, padx=100, pady=20)
        level_frame.pack(fill=tk.X, padx=100, pady=20)

        ttk.Button(self.init_frame, text="Forward", command=self.prepare_section_frames).pack(pady=50)

        self.init_frame.pack(expand=1, fill="both")

    def prepare_section_frames(self):
        """ Parse DEFAULTS and depends on it create Notebook tabs and entries """
        if self.config_filename.get() is "":
            tkMessageBox.showerror("Error", "You must choose file")
            return

        for section in utils.DEFAULTS:
            section_name = section["config_name"]
            section_name_long = section["name_long"]
            entries = section["entries"]

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

        recluster_add_button = ttk.Button(self.frames[self.recluster_frame_index], text="Add reclustering section")
        # Setting row=1000 let skip calculating position of button each time new section is added
        recluster_add_button.grid(row=1000, column=0, columnspan=2, pady=20)

        recluster_add_section_callback = utils.CallbackWrapper(self.callback_add_section, "reclusterization")
        recluster_add_button.bind("<Button-1>", recluster_add_section_callback)

        # Refresh all hiding frames
        self.refresh_menus()

        self.init_frame.pack_forget()

        self.notebook.pack(expand=1, fill="both")

        bottom_frame = ttk.Frame(self.parent)
        bottom_frame.pack(fill=tk.X)

        save_button = ttk.Button(bottom_frame, text="Save")
        save_button.pack(pady=5)

        save_button.bind("<Button-1>", self.save_config)

        # After preparing all section load values from config
        self.load_config_values()

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
        for i, item in enumerate(utils.MENUS):
            if item.startswith(section_name):
                del utils.MENUS[i]

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
        Append new frame with copy of the parent's entires.

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
        utils.MENUS.append("{}:method".format(section_name))

        defaults = utils.get_default_section(default_section_name)

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

        self.entry_filler(inner_frame, section_name, defaults["entries"], 2)

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

        Disables inlets_clusterization:max_level when it don't match the configuration level.

        :param parent: Parent widget.
        :param section_name: Config section name.
        :param entries: List of entries.
        :param row: Row number at which first entry will be shown.
        """
        entries_appended = 0
        for entry in entries:

            entry_name, entry_name_long, default, help, method, level = entry

            # If its primary reclustering or clustering section disable changing name of that section
            # Additionaly it disable inlets_clusterization:max_level if its on different level than chose
            state = tk.NORMAL
            if section_name == "clusterization" and entry_name == "name":
                state = tk.DISABLED
            elif section_name == "reclusterization" and entry_name == "name":
                state = tk.DISABLED
            elif section_name == "inlets_clusterization" and entry_name == "max_level":
                if utils.LEVELS[self.level.get()] != level:
                    state = tk.DISABLED

            # Skip entries that dont match level except inlets_clusterization:max_level
            if utils.LEVELS[self.level.get()] > level:
                if not (section_name == "inlets_clusterization" and entry_name == "max_level"):
                    continue

            entries_appended += 1

            if method:
                if section_name not in self.hiding_frames:
                    self.hiding_frames[section_name] = {}

                if method not in self.hiding_frames[section_name]:
                    if len(self.hiding_frames[section_name]) != 0:
                        hiding_frame_row = next(self.hiding_frames[section_name].itervalues()).row
                    else:
                        hiding_frame_row = row

                    self.hiding_frames[section_name][method] = utils.HidingFrame(parent,
                                                                                 hiding_frame_row,
                                                                                 text=method.capitalize() + " options",
                                                                                 style="HF.TFrame")

                hiding_frame = self.hiding_frames[section_name][method]

                if section_name not in self.values:
                    self.values[section_name] = {}

                self.values[section_name][entry_name] = utils.entry_factory(hiding_frame,
                                                                            hiding_frame.inner_row,
                                                                            entry_name_long,
                                                                            default, help, state)

                hiding_frame.inner_row += 1
            else:
                entry = utils.entry_factory(parent, row, entry_name_long, default, help, state)

                if (section_name.startswith("clusterization") or section_name.startswith(
                        "reclusterization")) and entry_name == "name":
                    entry.set(section_name)

                if utils.is_menu(section_name, entry_name):
                    opt_menu_callback = utils.CallbackWrapper(self.option_menu_changed, section_name, entry)
                    entry.input_var.trace("w", opt_menu_callback)

                if section_name not in self.values:
                    self.values[section_name] = {}

                self.values[section_name][entry_name] = entry

                row += 1

        # Hide tab if unused
        if entries_appended == 0:
            self.notebook.forget(self.scrolled_frames[parent])

    def refresh_menus(self):
        """
        Sets visible only hiding frames that are chosen by option menu.
        """
        for menu in utils.MENUS:
            section_name, entry_name = menu.split(":")

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

    def load_config_values(self):
        """
        Load values to values dictionary for each option from config file.


        Additionally it ask user if inlets_clusterization:max_level in config is different from found recursive_clusterization sections.
        """
        with open(self.config_filename.get(), "r") as config_file:
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

    def get_active_frame_name(self, section_name):
        """
        Return currently visible frame.

        :param section_name: Section name in which hiding frames will be considered.
        :return: Name of chosen option menu value.
        """
        for method, frame in self.hiding_frames[section_name].iteritems():
            if any(frame.grid_info()):
                return method

    def save_config(self, e):
        """
        Save all values to previously specified config file.

        Used as callback to Save button.

        :param e: Event informations.
        """
        if not self.values["global"]["top"].get():
            tkMessageBox.showerror("Topology file name", "You must specify topology file")
            return

        if not self.values["global"]["trj"].get():
            tkMessageBox.showerror("Trajectory file name", "You must specify at least one trajectory file")
            return

        with open(self.config_filename.get(), "w+") as config_file:
            config = ConfigParser()
            config.readfp(config_file)

            for section in self.values:
                if not config.has_section(section):
                    config.add_section(section)

                for option, value in self.values[section].iteritems():
                    default_value = utils.get_default_entry(section, option)[2]
                    if isinstance(default_value, tuple):
                        default_value = default_value[0]
                    elif isinstance(default_value, list):
                        if isinstance(default_value[0], tuple):
                            default_value = default_value[0][0]
                        # If second value is bool its checkbox.
                        # If checkbox is False there is no matter what is in input widget
                        elif isinstance(default_value[1], bool):
                            default_value = default_value[1]
                            if default_value:
                                default_value = default_value[0]
                        else:
                            default_value = default_value[0]

                    # Skip empty strings and entries which are not visible
                    if value.get() != "" and value.get() != default_value:
                        if not utils.get_default_entry(section, option)[4]:
                            config.set(section, option, value.get())
                        # This take care of saving options only from visible frame
                        else:
                            if utils.get_default_entry(section, option)[4] == self.get_active_frame_name(section):
                                config.set(section, option, value.get())

            config.write(config_file)

        tkMessageBox.showinfo("Saved", "Saving complete")


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

    app = ValveConfigApp(root)
    root.mainloop()
