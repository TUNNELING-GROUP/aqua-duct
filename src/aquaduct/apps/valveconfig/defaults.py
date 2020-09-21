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

import copy

class longstr(str):
    """
    Class used to specify type of default value.

    Represents Text field.
    """
    pass


class filetype(str):
    """
    Class used to specify type of default value.

    Represents Entry with file loading button.
    """
    pass


class manyfiletype(str):
    """
    Class used to specify type of default value.

    Represents Entry with file loading button which duplicate itself when previous is loaded.
    """
    pass


class dirtype(str):
    """
    Class used to specify type of default value.

    Represents Entry with dir loading button.
    """
    pass


class DefaultSection(object):
    def __init__(self, config_name, name, level, abs_level=None, additional=False):
        """
        Contrains info about section necessary to create it.

        :param config_name: Name of section in config.
        :param name: Brief label text, which will be displayed in Notebook tab.
        :param level: Entry level.
        :param abs_level: If specified section will be shown only when corresponding level was chosen.
        :param additional: If True values from that section wont be considered as default. Useful when section is mix of other sections.
        """
        self.config_name = config_name
        self.name = name
        self.level = level
        self.abs_level = abs_level
        self.additional = additional

        self.entries = []

        self._nested = False

    def add_entry(self, entry):
        """
        Allow to add new entry to the section or nest other section.
        Nested section will be showed in LabelFrame.
        :param entry: Entry or section
        :type entry DefaultEntry, DefaultSection
        """
        if isinstance(entry, DefaultEntry) or isinstance(entry, DefaultSection):
            self.entries.append(entry)
        else:
            raise TypeError("entry must be DefaultEntry or DefaultSection type.")

    def iter_entries(self):
        """
        Iterates through section entries.
        """
        def iter_inner(section):
            for entry in section.entries:
                if isinstance(entry, DefaultSection):
                    iter_inner(entry)
                else:
                    yield section, entry

        return iter_inner(self)


class DefaultEntry(object):
    def __init__(self, config_name, name, default_values, help_text, level, group_label=None, info_text=None,
                 warning_text=None, optionmenu_value=None, required=None):
        """
        Contains info about entry necessary to create it.

        :param config_name: Name of option in config.
        :param name: Brief label text, which will be displayed near widget.
        :param default_values: List of default values.
        :param help_text: Tooltip text.
        :param level: Entry level. Check LEVELS dict for adjust it.
        :param group_label: Used to group labels into frames. Content is a title of frame.
        :param info_text: If present information icon with content of that variable will be displayed.
        :param warning_text: If present warning icon with content of that variable will be displayed.
        :param optionmenu_value: OptionMenu value, which will cause display appropriate hiding frame.
        """
        self.config_name = config_name
        self.name = name

        if not isinstance(default_values, list):
            raise TypeError("Defaults values should be in the list")

        self.default_values = default_values
        self.help_text = help_text
        self.level = level
        self.group_label = group_label

        if info_text and warning_text:
            raise ValueError("There is no posibility to specify information and warning text.")

        self.info_text = info_text
        self.warning_text = warning_text
        self.optionmenu_value = optionmenu_value
        self.required = required

    @property
    def default_value(self):
        """
        Determines default value based on specified default values list.

        :return: Entry default value.
        """
        if len(self.default_values) == 1:
            if isinstance(self.default_values[0], tuple):
                default_value = self.default_values[0][0]
            elif isinstance(self.default_values[0], list):
                default_value = self.default_values[0][0]
            else:
                default_value = self.default_values[0]
        elif len(self.default_values) == 2:
            if isinstance(self.default_values[0], tuple):
                default_value = self.default_values[0][0]
            elif isinstance(self.default_values[0], list):
                default_value = self.default_values[0][0]
            # If second value is bool its checkbox.
            # If checkbox is False there is no matter what is in input widget
            elif isinstance(self.default_values[1], bool):
                default_value = self.default_values[1]
                if default_value:
                    default_value = self.default_values[0]
            else:
                default_value = self.default_values[0]

        return default_value


def is_menu(section, option):
    """
    Determines if option in section control hiding frames.

    :return: True if option control hiding frames, otherwise False.
    :rtype: bool
    """
    s = ":".join((section, option))
    for item in MENUS:
        if s == item:
            return True

    return False


def get_default_entry(section_name, option_name):
    """
    Return default entries.

    :param section_name: Name of section where option is located.
    :param option_name: Option name which default values are demanded.
    :return: DefaultEntry
    """
    for entry in get_default_section(section_name).entries:
        if option_name == entry.config_name:
            return entry

    raise RuntimeError("Entry {} in {} section does not exists.".format(option_name, section_name))


def get_default_section(section_name):
    """
    Return default section informations.

    :param section_name: Name of section which informations are demaned.
    :return: Default section informations.
    :rtype: DefaultSection
    """
    if section_name.startswith("clustering"):
        section_name = "clustering"

    if section_name.startswith("reclustering"):
        section_name = "reclustering"

    for section in DEFAULTS:
        if section.config_name == section_name and not section.additional:
            return section

    raise RuntimeError("Section {} does not exists.".format(section_name))


"""
Possible combinations:
str() -> Regular Entry
float() -> Smaller Entry
int() -> Smaller entry
bool() -> Checkbutton
filetype() -> Entry with file loading button
manyfiletype() -> Entry with file loading button with ability to duplicate itself
tuple() -> Option menu widget
any widget, bool() -> widget with checkbox
widget, float() -> create "widget_value(float_value)"

To create hiding frame set optionmenu_value to specific value and add info into MENUS list

Grouped entries will be displayed at place where first entry of that group will occur
Grouping does not work in entries with specified optionmenu_value
"""
# @formatter:off
DEFAULTS = []


###
# Initial, Easy
###
initial_section = DefaultSection(config_name="global", name="Initial", level=2, abs_level=2, additional=True)
initial_section.add_entry(DefaultEntry(config_name="top",
                                       name="Topology file: ",
                                       default_values=[manyfiletype()],
                                       help_text="Path to topology file. Aqua-Duct supports PDB, PRMTOP, PFS topology files.",
                                       level=2,
                                       required=1))
initial_section.add_entry(DefaultEntry(config_name="trj",
                                       name="Trajectory file: ",
                                       default_values=[manyfiletype()],
                                       help_text="Path to trajectory file. Aqua-Duct supports NC and DCD trajectory files.",
                                       level=2,
                                       required=1))
initial_section.add_entry(DefaultEntry(config_name="cache_dir",
                                       name="Cache directory: ",
                                       default_values=[dirtype("cache")],
                                       help_text="Allows to set path to the directory for cache data.",
                                       level=2))
initial_section.add_entry(DefaultEntry(config_name="cache_type",
                                       name="Cache type: ",
                                       default_values=[dirtype("full")],
                                       help_text="Allows to set type of coordinates cache.",
                                       level=2))

traceable_residues_nested = DefaultSection(config_name="traceable_residues", name="Traceable residues", level=2,
                                           additional=True)
traceable_residues_nested.add_entry(DefaultEntry(config_name="scope",
                                                 name="Scope: ",
                                                 default_values=["protein"],
                                                 help_text="Definition of Scope of interest.",
                                                 level=2,
                                                 required=1))
traceable_residues_nested.add_entry(DefaultEntry(config_name="object",
                                                 name="Object: ",
                                                 default_values=[str()],
                                                 help_text="Definition of Object of interest.",
                                                 level=2,
                                                 required=1))

initial_section.add_entry(traceable_residues_nested)

separate_paths_nested = DefaultSection(config_name="separate_paths", name="Separate paths", level=2, additional=True)
separate_paths_nested.add_entry(DefaultEntry(config_name="auto_barber",
                                             name="Auto Barber: ",
                                             default_values=[str()],
                                             help_text="This option allows to select molecular entity used in Auto Barber procedure. ",
                                             level=2))
separate_paths_nested.add_entry(DefaultEntry(config_name="execute",
                                              name="Execute: ",
                                              default_values=[("runonce", "run", "skip")],
                                              help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run calculations are always performed and if dump is set dump file is saved. If it is set to runonce calculations are performed if there is no dump file specified by dump option. If it is present calculations are skipped and data is loaded from the file. If it is set to skip calculations are skip and if dump is set data is loaded from the file.",
                                              level=1,
                                              info_text=" "))
initial_section.add_entry(separate_paths_nested)

DEFAULTS.append(initial_section)

###
# Clustering, Easy
###
clust_easy_section = DefaultSection(config_name="inlets_clustering", name="Clustering", level=2, abs_level=2, additional=True)
clust_easy_section.add_entry(DefaultEntry(config_name="detect_outliers",
                                          name="Detect outliers: ",
                                          default_values=[["False", "Auto"]],
                                          help_text="If set, detection of outliers is executed. It could be set as a floating point distance threshold or set to Auto.",
                                          level=1))
clust_easy_section.add_entry(DefaultEntry(config_name="singletons_outliers",
                                          name="Singletons outliers: ",
                                          default_values=["False"],
                                          help_text="Maximal size of cluster to be considered as outliers. If set to number > 0 clusters of that size are removed and their objects are moved to outliers.",
                                          level=2))
clust_easy_section.add_entry(DefaultEntry(config_name="max_level",
                                          name="Max level: ",
                                          default_values=[0],
                                          help_text="Maximal number of recursive clustering levels.",
                                          level=0,
                                          info_text=" "))

DEFAULTS.append(clust_easy_section)

###
# Global, Common
###
global_section = DefaultSection(config_name="global", name="General options", level=1)
global_section.add_entry(DefaultEntry(config_name="top",
                                      name="Topology file: ",
                                      default_values=[manyfiletype()],
                                      help_text="Path to topology file. Aqua-Duct supports PDB, PRMTOP, PFS topology files.",
                                      level=1,
                                      required=1))
global_section.add_entry(DefaultEntry(config_name="trj",
                                      name="Trajectory file: ",
                                      default_values=[manyfiletype()],
                                      help_text="Path to trajectory file. Aqua-Duct supports NC and DCD trajectory files.",
                                      level=1,
                                      required=1))
global_section.add_entry(DefaultEntry(config_name="twoway",
                                      name="Two-way scanning: ",
                                      default_values=[True],
                                      help_text="Try to use two-way scanning in the stage II.",
                                      level=0,
                                      info_text=" "))
global_section.add_entry(DefaultEntry(config_name="sandwich",
                                      name="Sandwich: ",
                                      default_values=[False],
                                      help_text="If set True trajectories are read as layers.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="waterfall",
                                      name="Waterfall: ",
                                      default_values=[False],
                                      help_text="If set True ends of trajectories are ends of paths.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="min_frame",
                                      name="Minimal frame: ",
                                      default_values=[0],
                                      help_text="Minimal number of frame to be read from trajectory data.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="max_frame",
                                      name="Maximal frame: ",
                                      default_values=[str()],
                                      help_text="Maximal number of frame to be read from trajectory data. If set None trajectory is read to the last possible frame.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="step_frame",
                                      name="Frame step: ",
                                      default_values=[1],
                                      help_text="Step used in reading trajectory. Default value of 1 stands for reading every frame. If it is greater than 1, only every step-value frame is read.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="sps",
                                      name="Single precision storage: ",
                                      default_values=[True],
                                      help_text="Try to store data in single precission storage.",
                                      level=0))
global_section.add_entry(DefaultEntry(config_name="cache_dir",
                                      name="Cache directory: ",
                                      default_values=[dirtype("cache")],
                                      help_text="Allows to set path to the directory for cache data.",
                                      level=1))
global_section.add_entry(DefaultEntry(config_name="cache_mem",
                                      name="Memory cache: ",
                                      default_values=[False],
                                      help_text="If set True, all data will be cached in RAM.",
                                      level=1))
DEFAULTS.append(global_section)

###
# Tracking section for normal level, mixin of traceable residues and separate paths
###
tracking_section = DefaultSection(config_name="traceable_residues",
                                  name="Tracking",
                                  level=1,
                                  abs_level=1,
                                  additional=True)

traceable_residues_section = DefaultSection(config_name="traceable_residues", name="Traceable residues", level=0)
traceable_residues_section.add_entry(DefaultEntry(config_name="execute",
                                                  name="Execute: ",
                                                  default_values=[("runonce", "run", "skip")],
                                                  help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run calculations are always performed and if dump is set dump file is saved. If it is set to runonce calculations are performed if there is no dump file specified by dump option. If it is present calculations are skipped and data is loaded from the file, If it is set to skip calculations are skip and if dump is set data is loaded from the file.",
                                                  level=1,
                                                  info_text=" "))
traceable_residues_section.add_entry(DefaultEntry(config_name="dump",
                                                  name="Dump file: ",
                                                  default_values=["1_traceable_residues_data.dump"],
                                                  help_text="File name of dump data. It is used to save results of calculations or to load previously calculated data - this depends on execute option.",
                                                  level=0))
traceable_residues_section.add_entry(DefaultEntry(config_name="scope",
                                                  name="Scope: ",
                                                  default_values=[str()],
                                                  help_text="Definition of Scope of interest.",
                                                  level=1,
                                                  required=1))
traceable_residues_section.add_entry(DefaultEntry(config_name="scope_convexhull",
                                                  name="Convex hull scope: ",
                                                  default_values=[True],
                                                  help_text="Flag to set if Scope is direct or convex hull definition.",
                                                  level=1))
traceable_residues_section.add_entry(DefaultEntry(config_name="scope_everyframe",
                                                  name="Everyframe scope: ",
                                                  default_values=[False],
                                                  help_text="Flag to set Scope evaluation mode. If set True Scope is evaluated in every frame. This make sense if the definition is complex and depends on distances between molecular entities.",
                                                  level=0,
                                                  warning_text="Could be time-consuming."))
traceable_residues_section.add_entry(DefaultEntry(config_name="scope_convexhull_inflate",
                                                  name="Scope convex hull inflate: ",
                                                  default_values=[str()],
                                                  help_text="Increase (or if negative - decrease) size of the scope convex hull.",
                                                  level=0,
                                                  warning_text="Could be time-consuming."))
traceable_residues_section.add_entry(DefaultEntry(config_name="object",
                                                  name="Object: ",
                                                  default_values=[str()],
                                                  help_text="Definition of Object of interest.",
                                                  level=1,
                                                  required=1))
traceable_residues_section.add_entry(DefaultEntry(config_name="add_passing",
                                                  name="Add passing: ",
                                                  default_values=[str()],
                                                  help_text="Definition of molecules that should be added to traced molecules even if they were not present in Object.",
                                                  level=0,
                                                  warning_text="Could be time-consuming."))
DEFAULTS.append(traceable_residues_section)

tracking_section.add_entry(copy.copy(traceable_residues_section))
tracking_section.add_entry(separate_paths_nested)

tracking_analysis_nested_section = DefaultSection(config_name="analysis",
                                                  name="Calculate scope and object size",
                                                  level=1,
                                                  abs_level=1,
                                                  additional=True)
tracking_analysis_nested_section.add_entry(DefaultEntry(config_name="calculate_scope_object_size",
                                                        name="Calculate scope and object size: ",
                                                        default_values=[False],
                                                        help_text="If set to True volumes and areas of object and scope approximated by convex hulls will be calculated for each analyzed frames and saved in output CSV file.",
                                                        level=1))
tracking_analysis_nested_section.add_entry(DefaultEntry(config_name="scope_chull",
                                                        name="Scope hull definition: ",
                                                        default_values=["protein"],
                                                        help_text="Scope convex hull definition used in calculating volume and area.",
                                                        level=1))
tracking_analysis_nested_section.add_entry(DefaultEntry(config_name="object_chull",
                                                        name="Object convex hull: ",
                                                        default_values=[str()],
                                                        help_text="Object convex hull definition used in calculating volume and area.",
                                                        level=1))
tracking_analysis_nested_section.add_entry(DefaultEntry(config_name="cluster_area",
                                                        name="Cluster area: ",
                                                        default_values=[True],
                                                        help_text="If set True, clusters’ areas are estimated with kernel density estimation method (KDE).",
                                                        level=1))

tracking_section.add_entry(tracking_analysis_nested_section)
DEFAULTS.append(tracking_section)

raw_paths_section = DefaultSection(config_name="raw_paths", name="Raw paths", level=1)
raw_paths_section.add_entry(DefaultEntry(config_name="execute",
                                         name="Execute: ",
                                         default_values=[("runonce", "run", "skip")],
                                         help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run calculations are always performed and if dump is set dump file is saved. If it is set to runonce calculations are performed if there is no dump file specified by dump option. If it is present calculations are skipped and data is loaded from the file, If it is set to skip calculations are skip and if dump is set data is loaded from the file.",
                                         level=0,
                                         info_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="dump",
                                         name="Dump file: ",
                                         default_values=["2_raw_paths_data.dump"],
                                         help_text="File name of dump data. It is used to save results of calculations or to load previously calculated data - this depends on execute option.",
                                         level=0))
raw_paths_section.add_entry(DefaultEntry(config_name="scope",
                                         name="Scope: ",
                                         default_values=[str()],
                                         help_text="Definition of Scope of interest. If None value form previous stage is used.",
                                         level=0,
                                         info_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="scope_convexhull",
                                         name="Convex hull scope: ",
                                         default_values=[str()],
                                         help_text="Flag to set if the Scope is direct or convex hull definition.",
                                         level=0,
                                         warning_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="scope_everyframe",
                                         name="Everyframe scope: ",
                                         default_values=[False],
                                         help_text="Flag to set Scope evaluation mode. If set True Scope is evaluated in every frame. This make sense if the definition is complex and depends on distances between molecular entities. If None value from previous stage is used.",
                                         level=0,
                                         warning_text="Could be time-consuming."))
raw_paths_section.add_entry(DefaultEntry(config_name="scope_convexhull_inflate",
                                         name="Scope convex hull inflate: ",
                                         default_values=[str()],
                                         help_text="Increase (or if negative - decrease) size of the scope convex hull. If None, value from previous stage is used.",
                                         level=0,
                                         warning_text="Could be time-consuming."))
raw_paths_section.add_entry(DefaultEntry(config_name="object",
                                         name="Object: ",
                                         default_values=[str()],
                                         help_text="Definition of Object of interest. If None value from the previous stage is used.",
                                         level=0,
                                         info_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="clear_in_object_info",
                                         name="Clear in object info: ",
                                         default_values=[False],
                                         help_text="If it is set to True information on occupation of Object site by traceable residues calculated in the previous stage is cleared and have to be recalculated. This is useful if definition of Object was changed.",
                                         level=0,
                                         warning_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="discard_singletons",
                                         name="Discard singletons: ",
                                         default_values=[1],
                                         help_text="If > 0, discards paths of given frame length.",
                                         level=0,
                                         info_text=" "))
raw_paths_section.add_entry(DefaultEntry(config_name="discard_empty_paths",
                                         name="Discard empty paths: ",
                                         default_values=[True],
                                         help_text="If set to True, empty paths are discarded.",
                                         level=0,
                                         info_text=" "))
DEFAULTS.append(raw_paths_section)

separate_paths_section = DefaultSection(config_name="separate_paths", name="Separate paths", level=0)
separate_paths_section.add_entry(DefaultEntry(config_name="execute",
                                              name="Execute: ",
                                              default_values=[("runonce", "run", "skip")],
                                              help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run calculations are always performed and if dump is set dump file is saved. If it is set to runonce calculations are performed if there is no dump file specified by dump option. If it is present calculations are skipped and data is loaded from the file. If it is set to skip calculations are skip and if dump is set data is loaded from the file.",
                                              level=1,
                                              info_text=" "))
separate_paths_section.add_entry(DefaultEntry(config_name="dump",
                                              name="Dump file: ",
                                              default_values=["3_separate_paths_data.dump"],
                                              help_text="File name of dump data. It is used to save results of calculations or to load previously calculated data - this depends on execute option.",
                                              level=0))
separate_paths_section.add_entry(DefaultEntry(config_name="discard_empty_paths",
                                              name="Discard empty paths: ",
                                              default_values=[True],
                                              help_text="If set to True empty paths are discarded.",
                                              level=0,
                                              info_text=" "))
separate_paths_section.add_entry(DefaultEntry(config_name="sort_by_id",
                                              name="Sort by ID: ",
                                              default_values=[True],
                                              help_text="If set to True separate paths are sorted by ID. Otherwise they are sorted in order of appearance.",
                                              level=0))
separate_paths_section.add_entry(DefaultEntry(config_name="discard_short_paths",
                                              name="Discard short paths: ",
                                              default_values=[20],
                                              help_text="This option allows to discard paths which are shorter than the threshold which is defined as total number of frames.",
                                              level=0,
                                              info_text=" "))
separate_paths_section.add_entry(DefaultEntry(config_name="discard_short_object",
                                              name="Discard short object: ",
                                              default_values=[2.0],
                                              help_text="This option allows to discard paths which objects are shorter than the threshold which is defined as total length in metric units.",
                                              level=0,
                                              warning_text=" "))
separate_paths_section.add_entry(DefaultEntry(config_name="discard_short_logic",
                                              name="Discard short logic: ",
                                              default_values=[("or", "and")],
                                              help_text="If both \"Discard short paths\" and \"Discard short object\" options are used, this option allows to set combination logic. If it is set to \"or\", a path is discarded if any of discard criterion is met. If it is set \"and\", both criteria have to be met to discard path.",
                                              level=0))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber",
                                              name="Auto Barber: ",
                                              default_values=[str()],
                                              help_text="This option allows to select molecular entity used in Auto Barber procedure. ",
                                              level=1,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber_mincut",
                                              name="Auto Barber mincut: ",
                                              default_values=[str()],
                                              help_text="Minimal radius of spheres used in Auto Barber. If a sphere has radius smaller then this value it is not used in AutoBarber procedure. This option can be switched off by setting it to None.",
                                              level=0,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut",
                                              name="Auto Barber maxcut: ",
                                              default_values=[1.2],
                                              help_text="Maximal radius of spheres used in Auto Barber. If a sphere has radius greater then this value it is not used in AutoBarber procedure. This option can be switched off by setting it to None.",
                                              level=0,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber_mincut_level",
                                              name="Auto Barber mincut level:",
                                              default_values=[True],
                                              help_text="If set True spheres of radius smaller than mincut are resized to mincut value.",
                                              level=0,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut_level",
                                              name="Auto Barber maxcut level:",
                                              default_values=[True],
                                              help_text="If set True spheres of radius greater than maxcut are resized to maxcut value.",
                                              level=0,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="auto_barber_tovdw",
                                              name="Auto Barber VdW: ",
                                              default_values=[True],
                                              help_text="Correct cutting sphere by decreasing its radius by VdW radius of the closest atom.",
                                              level=0,
                                              group_label="Auto Barber"))
separate_paths_section.add_entry(DefaultEntry(config_name="allow_passing_paths",
                                              name="Allow passing paths: ",
                                              default_values=[False],
                                              help_text="If set True paths that do not enter the object are detected and added to the rest of paths as ‘passing’ paths.",
                                              level=0,
                                              warning_text="Required for hotspots and energy profiles calculations."))
separate_paths_section.add_entry(DefaultEntry(config_name="calculate_coo",
                                              name="Calculate center of object: ",
                                              default_values=[False],
                                              help_text="If set True center of object is calculated by averaging coordinates of all paths in the object area.",
                                              level=0,
                                              warning_text="Calculation of CoO may be very long for some systems."))
separate_paths_section.add_entry(DefaultEntry(config_name="separate_barber",
                                              name="Separate barber: ",
                                              default_values=[True],
                                              help_text="Apply AutoBarber for each type of traced molecules separately.",
                                              level=0))

DEFAULTS.append(separate_paths_section)

inlets_clustering_section = DefaultSection(config_name="inlets_clustering", name="Inlets clustering",
                                           level=1)
inlets_clustering_section.add_entry(DefaultEntry(config_name="execute",
                                                 name="Execute: ",
                                                 default_values=[("runonce", "run", "skip")],
                                                 help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run calculations are always performed and if dump is set dump file is saved. If it is set to runonce calculations are performed if there is no dump file specified by dump option. If it is present calculations are skipped and data is loaded from the file. If it is set to skip calculations are skip and if dump is set data is loaded from the file.",
                                                 level=1,
                                                 info_text=" "))
inlets_clustering_section.add_entry(DefaultEntry(config_name="dump",
                                                 name="Dump file: ",
                                                 default_values=["4_inlets_clustering_data.dump"],
                                                 help_text="File name of dump data. It is used to save results of calculations or to load previously calculated data - this depends on execute option.",
                                                 level=0))
inlets_clustering_section.add_entry(DefaultEntry(config_name="recluster_outliers",
                                                 name="Recluster outliers: ",
                                                 default_values=[False],
                                                 help_text="If set to True reclustering of outliers is executed according to the method defined in reclustering section.",
                                                 level=0,
                                                 info_text=" "))
inlets_clustering_section.add_entry(DefaultEntry(config_name="detect_outliers",
                                                 name="Detect outliers: ",
                                                 default_values=[["False", "Auto"]],
                                                 help_text="If set, detection of outliers is executed. It could be set as a floating point distance threshold or set to Auto.",
                                                 level=1))
inlets_clustering_section.add_entry(DefaultEntry(config_name="singletons_outliers",
                                                 name="Singletons outliers: ",
                                                 default_values=["False"],
                                                 help_text="Maximal size of cluster to be considered as outliers. If set to number > 0 clusters of that size are removed and their objects are moved to outliers.",
                                                 level=1))
inlets_clustering_section.add_entry(DefaultEntry(config_name="max_level",
                                                 name="Max level: ",
                                                 default_values=[0],
                                                 help_text="Maximal number of recursive clustering levels.",
                                                 level=0,
                                                 info_text=" "))
inlets_clustering_section.add_entry(DefaultEntry(config_name="create_master_paths",
                                                 name="Create master paths: ",
                                                 default_values=[True],
                                                 help_text="If set to True master paths are created (fast CPU and big RAM recommended; 50k frames long simulation may need ca 20GB of memory)",
                                                 level=1,
                                                 warning_text="Required for energy profiles calculations and could be time-consuming"))
inlets_clustering_section.add_entry(DefaultEntry(config_name="master_paths_amount",
                                                 name="Master paths amount",
                                                 default_values=[str()],
                                                 help_text="Allows to limit number of single paths used for master paths calculations. If it is a number in range ``(0,1)``, then it is interpreted as a percent number of paths to be used. It is is a integer number ``>= 1`` it is an absoulte number of used paths.",
                                                 level=0))
inlets_clustering_section.add_entry(DefaultEntry(config_name="separate_master",
                                                 name="Separate master paths",
                                                 default_values=[False],
                                                 help_text="If set to True, master paths are generated for each of traced molecules' types separately.",
                                                 level=0))
inlets_clustering_section.add_entry(DefaultEntry(config_name="separate_master_all",
                                                 name="Separate master paths all",
                                                 default_values=[True],
                                                 help_text="If separate_master is used and this option is set True, master paths will be also generated for all traced molecules",
                                                 level=0))
inlets_clustering_section.add_entry(DefaultEntry(config_name="exclude_passing_in_clustering",
                                                 name="Exclude passing in clustering: ",
                                                 default_values=[True],
                                                 help_text="If set to True passing paths are not clustered with normal paths.",
                                                 level=0,
                                                 group_label="Passing paths",
                                                 info_text=" "))
inlets_clustering_section.add_entry(DefaultEntry(config_name="add_passing_to_clusters",
                                                 name="Add passing to clusters: ",
                                                 default_values=[str()],
                                                 help_text="Allows to run procedure for adding passing paths inlets to clusters with Auto Barber method. To enable this the option should be set to molecular entity that will be used by Auto Barber.",
                                                 level=0,
                                                 group_label="Passing paths",
                                                 info_text=" "))
inlets_clustering_section.add_entry(DefaultEntry(config_name="join_clusters",
                                                 name="Join clusters: ",
                                                 default_values=[str()],
                                                 help_text="This option allows to join selected clusters. Clusters’ IDs joined with + character lists clusters to be joined together. Several such blocks separated by space can be used. For example, if set to 1+3+4 5+6 clusters 1, 3, and 4 will be joined in one cluster and cluster 5, and 6 will be also joined in another cluster.",
                                                 level=1,
                                                 group_label="Post Clustering Optimization"))
inlets_clustering_section.add_entry(DefaultEntry(config_name="renumber_clusters",
                                                 name="Renumber clusters: ",
                                                 default_values=[False],
                                                 help_text="If set True, clusters have consecutive numbers starting from 1 (or 0 if outliers are present) starting from the bigest cluster.",
                                                 level=1,
                                                 group_label="Post Clustering Optimization"))
inlets_clustering_section.add_entry(DefaultEntry(config_name="inlets_center",
                                                 name="Central point for inlets: ",
                                                 default_values=[['cos','coo']],
                                                 help_text="Central point for inlets influencs cluster area calculations and cluster countours.",
                                                 warning_text="If `coo` is selected `calculate_coo` should be set `True`. ",
                                                 level=1))
inlets_clustering_section.add_entry(DefaultEntry(config_name="remove_inlets",
                                                 name="Remove inlets from specified clusters: ",
                                                 default_values=[""],
                                                 help_text="Inlets in specified clusters are removed - corresponding paths are changed.",
                                                 warning_text="",
                                                 level=1))
inlets_clustering_section.add_entry(DefaultEntry(config_name="clustering_order",
                                                 name="Order of clustering operations: ",
                                                 default_values=[['old-school','aquarius']],
                                                 help_text="Order in which clustering related operations are executed.",
                                                 warning_text="",
                                                 level=1))



DEFAULTS.append(inlets_clustering_section)

clustering_section = DefaultSection(config_name="clustering", name="Clustering", level=1)
clustering_section.add_entry(DefaultEntry(config_name="name",
                                          name="Name: ",
                                          default_values=["clustering"],
                                          help_text="Used to refer other clustering method in \"Recursive clustering\" option",
                                          level=2))
clustering_section.add_entry(DefaultEntry(config_name="method",
                                          name="Method: ",
                                          default_values=[
                                              ("barber", "dbscan", "affprop", "meanshift", "birch", "kmeans")],
                                          help_text="Name of clustering method. ",
                                          level=2))

# Barber options
clustering_section.add_entry(DefaultEntry(config_name="auto_barber",
                                          name="Auto barber: ",
                                          default_values=[str()],
                                          help_text="This option allows to select molecular entity used in Auto Barber procedure.",
                                          level=2,
                                          optionmenu_value="barber"))
clustering_section.add_entry(DefaultEntry(config_name="auto_barber_mincut",
                                          name="Auto Barber mincut: ",
                                          default_values=[str()],
                                          help_text="Minimal radius of spheres used in Auto Barber. If a sphere has radius smaller than this value, it is not used to cut. This option can be switched off by setting it to None.",
                                          level=0,
                                          optionmenu_value="barber"))
clustering_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut",
                                          name="Auto Barber maxcut: ",
                                          default_values=[2.8],
                                          help_text="Maximal radius of spheres used in Auto Barber. If a sphere has radius greater than this value, it is not used to cut. This option can be switched off by setting it to None.",
                                          level=0,
                                          optionmenu_value="barber"))
clustering_section.add_entry(DefaultEntry(config_name="auto_barber_mincut_level",
                                          name="Auto Barber mincut level: ",
                                          default_values=[True],
                                          help_text="If set True, spheres of radius less than mincut are resized to mincut value.",
                                          level=0,
                                          optionmenu_value="barber"))
clustering_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut_level",
                                          name="Auto Barber maxcut level:",
                                          default_values=[True],
                                          help_text="If set True, spheres of radius greater than maxcut are resized to maxcut value.",
                                          level=0,
                                          optionmenu_value="barber"))
clustering_section.add_entry(DefaultEntry(config_name="auto_barber_tovdw",
                                          name="Auto Barber VdW: ",
                                          default_values=[True],
                                          help_text="If set True, cutting of spheres is corrected by decreasing its radius by VdW radius of the closest atom.",
                                          level=0,
                                          optionmenu_value="barber"))

# Dbscan options
clustering_section.add_entry(DefaultEntry(config_name="eps",
                                          name="Maximum distance: ",
                                          default_values=[float()],
                                          help_text="The maximum distance between two samples for them to be considered as in the same neighborhood.",
                                          level=0,
                                          optionmenu_value="dbscan"))
clustering_section.add_entry(DefaultEntry(config_name="min_samples",
                                          name="Maxium samples: ",
                                          default_values=[int()],
                                          help_text="The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.",
                                          level=0, optionmenu_value="dbscan"))
clustering_section.add_entry(DefaultEntry(config_name="metric",
                                          name="Metric: ",
                                          default_values=[("euclidean", "cityblock", "cosine", "manhattan")],
                                          help_text="The metric to use when calculating distance between instances in a feature array.",
                                          level=0,
                                          optionmenu_value="dbscan"))
clustering_section.add_entry(DefaultEntry(config_name="algorithm",
                                          name="Algorithm: ",
                                          default_values=[("auto", "ball_tree", "kd_tree", "brute")],
                                          help_text="The algorithm to be used by the NearestNeighbors module to compute pointwise distances and find nearest neighbors.",
                                          level=0,
                                          optionmenu_value="dbscan"))
clustering_section.add_entry(DefaultEntry(config_name="leaf_size",
                                          name="Leaf size: ",
                                          default_values=[int()],
                                          help_text="Leaf size passed to BallTree or cKDTree.",
                                          level=0,
                                          optionmenu_value="dbscan"))

# Affprop options
clustering_section.add_entry(DefaultEntry(config_name="damping",
                                          name="Damping factor: ",
                                          default_values=[float()],
                                          help_text="Damping factor between 0.5 and 1.",
                                          level=0,
                                          optionmenu_value="affprop"))
clustering_section.add_entry(DefaultEntry(config_name="convergence_iter",
                                          name="Maximum no effect iterations: ",
                                          default_values=[int()],
                                          help_text="Number of iterations with no change in the number of estimated clusters that stops the convergence.",
                                          level=0,
                                          optionmenu_value="affprop"))
clustering_section.add_entry(DefaultEntry(config_name="max_iter",
                                          name="Maximum number of iterations: ",
                                          default_values=[int()],
                                          help_text="Maximum number of iterations.",
                                          level=0,
                                          optionmenu_value="affprop"))
clustering_section.add_entry(DefaultEntry(config_name="preference",
                                          name="Preference: ",
                                          default_values=[float()],
                                          help_text="Points with larger values of preferences are more likely to be chosen as exemplars.",
                                          level=0,
                                          optionmenu_value="affprop"))

# Meanshift options
clustering_section.add_entry(DefaultEntry(config_name="bandwidth",
                                          name="Bandwidth: ",
                                          default_values=["Auto"],
                                          help_text="Bandwidth used in the RBF kernel. If Auto or None automatic method for bandwidth estimation is used.",
                                          level=2,
                                          optionmenu_value="meanshift"))
clustering_section.add_entry(DefaultEntry(config_name="cluster_all",
                                          name="Cluster all points: ",
                                          default_values=[bool()],
                                          help_text="If true, then all points are clustered, even those orphans that are not within any kernel.",
                                          level=0,
                                          optionmenu_value="meanshift"))
clustering_section.add_entry(DefaultEntry(config_name="bin_seeding",
                                          name="bin_seeding",
                                          default_values=[bool()],
                                          help_text="If true, initial kernel locations are not locations of all points, but rather the location of the discretized version of points, where points are binned onto a grid whose coarseness corresponds to the bandwidth.",
                                          level=0,
                                          optionmenu_value="meanshift"))
clustering_section.add_entry(DefaultEntry(config_name="min_bin_freq",
                                          name="min_bin_freq",
                                          default_values=[int()],
                                          help_text="To speed up the algorithm, accept only those bins with at least min_bin_freq points as seeds. If not defined, set to 1.",
                                          level=0,
                                          optionmenu_value="meanshift"))

# Birch options
clustering_section.add_entry(DefaultEntry(config_name="threshold",
                                          name="Threshold: ",
                                          default_values=[float()],
                                          help_text="The radius of the subcluster obtained by merging a new sample and the closest subcluster should be smaller than the threshold. Otherwise a new subcluster is started.",
                                          level=0,
                                          optionmenu_value="birch"))
clustering_section.add_entry(DefaultEntry(config_name="branching_factor",
                                          name="Branching factor: ",
                                          default_values=[int()],
                                          help_text="Maximum number of CF subclusters in each node.",
                                          level=0,
                                          optionmenu_value="birch"))
clustering_section.add_entry(DefaultEntry(config_name="n_clusters",
                                          name="Cluster number: ",
                                          default_values=[int()],
                                          help_text="Number of clusters after the final clustering step, which treats the subclusters from the leaves as new samples. By default, this final clustering step is not performed and the subclusters are returned as they are.",
                                          level=2,
                                          optionmenu_value="birch"))

# Kmeans options
clustering_section.add_entry(DefaultEntry(config_name="n_clusters",
                                          name="Cluster number: ",
                                          default_values=[int()],
                                          help_text="The number of clusters to form as well as the number of centroids to generate.",
                                          level=2,
                                          optionmenu_value="kmeans"))
clustering_section.add_entry(DefaultEntry(config_name="max_iter",
                                          name="Maximum number of iterations: ",
                                          default_values=[int()],
                                          help_text="Maximum number of iterations of the k-means algorithm for a single run.",
                                          level=0,
                                          optionmenu_value="kmeans"))
clustering_section.add_entry(DefaultEntry(config_name="n_init",
                                          name="n_init",
                                          default_values=[int()],
                                          help_text="Number of times the k-means algorithm will be run with different centroid seeds. The final results will be the best output of n_init consecutive runs in terms of inertia.",
                                          level=0,
                                          optionmenu_value="kmeans"))
clustering_section.add_entry(DefaultEntry(config_name="init",
                                          name="Init method: ",
                                          default_values=[("k-means++", "random")],
                                          help_text="Method for initialization, defaults to k-means++. Can be one of following: k-means++ or random.",
                                          level=0,
                                          optionmenu_value="kmeans"))
clustering_section.add_entry(DefaultEntry(config_name="tol",
                                          name="Tolerance: ",
                                          default_values=[float()],
                                          help_text="Relative tolerance with regards to inertia to declare convergence.",
                                          level=0,
                                          optionmenu_value="kmeans"))

clustering_section.add_entry(DefaultEntry(config_name="recursive_clustering",
                                          name="Recursive clustering: ",
                                          default_values=["clustering"],
                                          help_text="If it is set to name of some section that holds clustering method settings this method will be called in the next recursion of clusteriation. Default value for reclustering is None.",
                                          level=1))
clustering_section.add_entry(DefaultEntry(config_name="recursive_threshold",
                                          name="Recursive threshold: ",
                                          default_values=[str()],
                                          help_text="Allows to set threshold that excludes clusters of certain size from reclustering. Value of this option comprises of operator and value. Operator can be one of the following: >, >=, <=, <. Value have to be expressed as floating number and it have to be in the range of 0 to 1. One can use several definitions separated by a space character. Only clusters of size complying with all thresholds definitions are submitted to reclustering.",
                                          level=1))


DEFAULTS.append(clustering_section)

reclustering_section = DefaultSection(config_name="reclustering", name="Reclustering", level=0)
reclustering_section.add_entry(DefaultEntry(config_name="name",
                                            name="Name: ",
                                            default_values=["reclustering"],
                                            help_text="Used to refer other clustering method in \"Recursive clustering\" option",
                                            level=1))
reclustering_section.add_entry(DefaultEntry(config_name="method",
                                            name="Method: ",
                                            default_values=[
                                                ("barber", "dbscan", "affprop", "meanshift", "birch", "kmeans")],
                                            help_text="Name of clustering method. It has to be one of the following: barber, dbscan, affprop, meanshift, birch, kmeans. Default value depends whether it is clustering section (barber) or reclustering section (dbscan).",
                                            level=1))

# Barber options
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber",
                                            name="Auto barber: ",
                                            default_values=[str()],
                                            help_text="This option allows to select molecular entity used in Auto Barber procedure.",
                                            level=1,
                                            optionmenu_value="barber"))
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber_mincut",
                                            name="Auto Barber mincut: ",
                                            default_values=[str()],
                                            help_text="Minimal radius of spheres used in Auto Barber. If a sphere has radius smaller than this value, it is not used to cut. This option can be switched off by setting it to None.",
                                            level=0,
                                            optionmenu_value="barber"))
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut",
                                            name="Auto Barber maxcut: ",
                                            default_values=[str()],
                                            help_text="Maximal radius of spheres used in Auto Barber. If a sphere has radius greater than this value, it is not used to cut. This option can be switched off by setting it to None.",
                                            level=0,
                                            optionmenu_value="barber"))
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber_mincut_level",
                                            name="Auto Barber mincut level: ",
                                            default_values=[bool()],
                                            help_text="If set True, spheres of radius less than mincut are resized to mincut value.",
                                            level=0,
                                            optionmenu_value="barber"))
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber_maxcut_level",
                                            name="Auto Barber maxcut level:",
                                            default_values=[bool()],
                                            help_text="If set True, spheres of radius greater than maxcut are resized to maxcut value.",
                                            level=0,
                                            optionmenu_value="barber"))
reclustering_section.add_entry(DefaultEntry(config_name="auto_barber_tovdw",
                                            name="Auto Barber VdW: ",
                                            default_values=[bool()],
                                            help_text="If set True, cutting of spheres is corrected by decreasing its radius by VdW radius of the closest atom.",
                                            level=0,
                                            optionmenu_value="barber"))

# Dbscan options
reclustering_section.add_entry(DefaultEntry(config_name="eps",
                                            name="Maximum distance: ",
                                            default_values=[float()],
                                            help_text="The maximum distance between two samples for them to be considered as in the same neighborhood.",
                                            level=0,
                                            optionmenu_value="dbscan"))
reclustering_section.add_entry(DefaultEntry(config_name="min_samples",
                                            name="Maximum samples: ",
                                            default_values=[int()],
                                            help_text="The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.",
                                            level=0, optionmenu_value="dbscan"))
reclustering_section.add_entry(DefaultEntry(config_name="metric",
                                            name="Metric: ",
                                            default_values=[("euclidean", "cityblock", "cosine", "manhattan")],
                                            help_text="The metric to use when calculating distance between instances in a feature array.",
                                            level=0,
                                            optionmenu_value="dbscan"))
reclustering_section.add_entry(DefaultEntry(config_name="algorithm",
                                            name="Algorithm: ",
                                            default_values=[("auto", "ball_tree", "kd_tree", "brute")],
                                            help_text="The algorithm to be used by the NearestNeighbors module to compute pointwise distances and find nearest neighbors.",
                                            level=0,
                                            optionmenu_value="dbscan"))
reclustering_section.add_entry(DefaultEntry(config_name="leaf_size",
                                            name="Leaf size: ",
                                            default_values=[int()],
                                            help_text="Leaf size passed to BallTree or cKDTree.",
                                            level=0,
                                            optionmenu_value="dbscan"))

# Affprop options
reclustering_section.add_entry(DefaultEntry(config_name="damping",
                                            name="Damping factor: ",
                                            default_values=[float()],
                                            help_text="Damping factor between 0.5 and 1.",
                                            level=0,
                                            optionmenu_value="affprop"))
reclustering_section.add_entry(DefaultEntry(config_name="convergence_iter",
                                            name="Maximum no effect iterations: ",
                                            default_values=[int()],
                                            help_text="Number of iterations with no change in the number of estimated clusters that stops the convergence.",
                                            level=0,
                                            optionmenu_value="affprop"))
reclustering_section.add_entry(DefaultEntry(config_name="max_iter",
                                            name="Maximum number of iterations: ",
                                            default_values=[int()],
                                            help_text="Maximum number of iterations.",
                                            level=0,
                                            optionmenu_value="affprop"))
reclustering_section.add_entry(DefaultEntry(config_name="preference",
                                            name="Preference: ",
                                            default_values=[float()],
                                            help_text="Points with larger values of preferences are more likely to be chosen as exemplars.",
                                            level=0,
                                            optionmenu_value="affprop"))

# Meanshift options
reclustering_section.add_entry(DefaultEntry(config_name="bandwidth",
                                            name="Bandwidth: ",
                                            default_values=["Auto"],
                                            help_text="Bandwidth used in the RBF kernel. If Auto or None automatic method for bandwidth estimation is used.",
                                            level=1,
                                            optionmenu_value="meanshift"))
reclustering_section.add_entry(DefaultEntry(config_name="cluster_all",
                                            name="Cluster all points: ",
                                            default_values=[bool()],
                                            help_text="If true, then all points are clustered, even those orphans that are not within any kernel.",
                                            level=0,
                                            optionmenu_value="meanshift"))
reclustering_section.add_entry(DefaultEntry(config_name="bin_seeding",
                                            name="Bin seeding: ",
                                            default_values=[bool()],
                                            help_text="If true, initial kernel locations are not locations of all points, but rather the location of the discretized version of points, where points are binned onto a grid whose coarseness corresponds to the bandwidth.",
                                            level=0,
                                            optionmenu_value="meanshift"))
reclustering_section.add_entry(DefaultEntry(config_name="min_bin_freq",
                                            name="Minimum bin frequency",
                                            default_values=[int()],
                                            help_text="To speed up the algorithm, accept only those bins with at least min_bin_freq points as seeds. If not defined, set to 1.",
                                            level=0,
                                            optionmenu_value="meanshift"))

# Birch options
reclustering_section.add_entry(DefaultEntry(config_name="threshold",
                                            name="Threshold: ",
                                            default_values=[float()],
                                            help_text="The radius of the subcluster obtained by merging a new sample and the closest subcluster should be smaller than the threshold. Otherwise a new subcluster is started.",
                                            level=0,
                                            optionmenu_value="birch"))
reclustering_section.add_entry(DefaultEntry(config_name="branching_factor",
                                            name="Branching factor: ",
                                            default_values=[int()],
                                            help_text="Maximum number of CF subclusters in each node.",
                                            level=0,
                                            optionmenu_value="birch"))
reclustering_section.add_entry(DefaultEntry(config_name="n_clusters",
                                            name="Cluster number: ",
                                            default_values=[int()],
                                            help_text="Number of clusters after the final clustering step, which treats the subclusters from the leaves as new samples. By default, this final clustering step is not performed and the subclusters are returned as they are.",
                                            level=1,
                                            optionmenu_value="birch"))

# Kmeans options
reclustering_section.add_entry(DefaultEntry(config_name="n_clusters",
                                            name="Cluster number: ",
                                            default_values=[int()],
                                            help_text="The number of clusters to form as well as the number of centroids to generate.",
                                            level=1,
                                            optionmenu_value="kmeans"))
reclustering_section.add_entry(DefaultEntry(config_name="max_iter",
                                            name="Maximum number of iterations: ",
                                            default_values=[int()],
                                            help_text="Maximum number of iterations of the k-means algorithm for a single run.",
                                            level=0,
                                            optionmenu_value="kmeans"))
reclustering_section.add_entry(DefaultEntry(config_name="n_init",
                                            name="N Init: ",
                                            default_values=[int()],
                                            help_text="Number of times the k-means algorithm will be run with different centroid seeds. The final results will be the best output of n_init consecutive runs in terms of inertia.",
                                            level=0,
                                            optionmenu_value="kmeans"))
reclustering_section.add_entry(DefaultEntry(config_name="init",
                                            name="Init method: ",
                                            default_values=[("k-means++", "random")],
                                            help_text="Method for initialization, defaults to k-means++. Can be one of following: k-means++ or random.",
                                            level=0,
                                            optionmenu_value="kmeans"))
reclustering_section.add_entry(DefaultEntry(config_name="tol",
                                            name="Tolerance: ",
                                            default_values=[float()],
                                            help_text="Relative tolerance with regards to inertia to declare convergence.",
                                            level=0,
                                            optionmenu_value="kmeans"))

reclustering_section.add_entry(DefaultEntry(config_name="recursive_clustering",
                                            name="Recursive clustering: ",
                                            default_values=[str()],
                                            help_text="If it is set to name of some section that holds clustering method settings this method will be called in the next recursion of clusteriation. Default value for reclustering is None.",
                                            level=1))
reclustering_section.add_entry(DefaultEntry(config_name="recursive_threshold",
                                            name="Recursive threshold: ",
                                            default_values=[str()],
                                            help_text="Allows to set threshold that excludes clusters of certain size from reclustering. Value of this option comprises of operator and value. Operator can be one of the following: >, >=, <=, <. Value have to be expressed as floating number and it have to be in the range of 0 to 1. One can use several definitions separated by a space character. Only clusters of size complying with all thresholds definitions are submitted to reclustering.",
                                            level=1))
DEFAULTS.append(reclustering_section)

analysis_section = DefaultSection(config_name="analysis", name="Analysis", level=0)
analysis_section.add_entry(DefaultEntry(config_name="execute",
                                        name="Execute: ",
                                        default_values=[("run", "runonce", "skip")],
                                        help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run or runonce stage is executed and results is saved according to save option. If it is set to skip stage is skipped.",
                                        level=1))
analysis_section.add_entry(DefaultEntry(config_name="save",
                                        name="Save file: ",
                                        default_values=["5_analysis_results.txt"],
                                        help_text="File name for saving results.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="dump_config",
                                        name="Dump config: ",
                                        default_values=[True],
                                        help_text="If set to True configuration options, as seen by Valve, are added to the head of results.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="calculate_scope_object_size",
                                        name="Calculate scope and object size: ",
                                        default_values=[False],
                                        help_text="If set to True volumes and areas of object and scope approximated by convex hulls will be calculated for each analyzed frames and saved in output CSV file.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="scope_chull",
                                        name="Scope hull definition: ",
                                        default_values=[str()],
                                        help_text="Scope convex hull definition used in calculating volume and area.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="scope_chull_inflate",
                                        name="Scope convexhull inflate: ",
                                        default_values=[str()],
                                        help_text="Increase (or if negative - decrease) size of the scope convex hull.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="object_chull",
                                        name="Object convex hull: ",
                                        default_values=[str()],
                                        help_text="Object convex hull definition used in calculating volume and area.",
                                        level=0))
analysis_section.add_entry(DefaultEntry(config_name="cluster_area",
                                        name="Cluster area: ",
                                        default_values=[True],
                                        help_text="If set True, clusters’ areas are estimated with kernel density estimation method (KDE).",
                                        level=1))

DEFAULTS.append(analysis_section)

visualize_section = DefaultSection(config_name="visualize", name="Visualize", level=2)
visualize_section.add_entry(DefaultEntry(config_name="execute",
                                         name="Execute: ",
                                         default_values=[("run", "runonce", "skip")],
                                         help_text="Option controls stage execution. It can have one of three possible values: run, runonce, and skip. If it is set to run or runonce stage is executed and results is saved according to save option. If it is set to skip stage is skipped.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="save",
                                         name="Save file: ",
                                         default_values=["6_visualize_results.py"],
                                         help_text="File name for saving results.",
                                         level=0))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_raw",
                                         name="All paths raw: ",
                                         default_values=[True],
                                         help_text="If True produces one object in PyMOL that holds all paths visualized by raw coordinates.",
                                         level=2,
                                         group_label="Raw paths",
                                         warning_text=" "))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_smooth",
                                         name="All paths smooth: ",
                                         default_values=[True],
                                         help_text="If True produces one object in PyMOL that holds all paths visualized by smooth coordinates.",
                                         level=2,
                                         group_label="Smooth paths"))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_split",
                                         name="All split paths: ",
                                         default_values=[True],
                                         help_text="If is set True objects produced by all_paths_raw and all_paths_smooth are split into Incoming, Object, and Outgoing parts and visualized as three different objects.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_raw_io",
                                         name="All raw paths io: ",
                                         default_values=[True],
                                         help_text="If set True arrows pointing beginning and end of paths are displayed oriented accordingly to raw paths orientation.",
                                         level=1,
                                         group_label="Raw paths",
                                         warning_text=" "))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_smooth_io",
                                         name="All paths smooth io: ",
                                         default_values=[True],
                                         help_text="If set True arrows pointing beginning and end of paths are displayed oriented accordingly to smooth paths orientation.",
                                         level=1,
                                         group_label="Smooth paths"))
visualize_section.add_entry(DefaultEntry(config_name="all_paths_amount",
                                         name="Paths limit: ",
                                         default_values=[str()],
                                         help_text="Allows to limit number of visualised paths. If it is a number in range (0,1), then it is interpreted as a percent number of paths to be visualized. It is is a integer number >= 1 it is total number of all_paths visualized.",
                                         level=0))
visualize_section.add_entry(DefaultEntry(config_name="simply_smooths",
                                         name="Smoothing simplification: ",
                                         default_values=[("RecursiveVector", "HobbitVector", "OneWayVector",
                                                          "RecursiveTriangle", "HobbitTriangle", "OneWayTriangle"),
                                                         float()],
                                         help_text="Option indicates linear simplification method to be used in plotting smooth paths. Simplification removes points which do not (or almost do not) change the shape of smooth path. Optionally name of the method can be followed by a threshold value. For sane values of thresholds see appropriate documentation of each method. Default values work well. This option is not case sensitive. It is recommended to use default method or HobbitVector method.",
                                         level=0))
visualize_section.add_entry(DefaultEntry(config_name="paths_raw",
                                         name="Paths raw: ",
                                         default_values=[True],
                                         help_text="If set True raw paths are displayed as separate objects or as one object with states corresponding to number of path.",
                                         level=2,
                                         group_label="Raw paths",
                                         warning_text=" "))
visualize_section.add_entry(DefaultEntry(config_name="paths_smooth",
                                         name="Paths smooth: ",
                                         default_values=[True],
                                         help_text="If set True smooth paths are displayed as separate objects or as one object with states corresponding to number of path.",
                                         level=2,
                                         group_label="Smooth paths"))
visualize_section.add_entry(DefaultEntry(config_name="paths_raw_io",
                                         name="Paths raw io: ",
                                         default_values=[True],
                                         help_text="If set True arrows indicating beginning and end of paths, oriented accordingly to raw paths, are displayed as separate objects or as one object with states corresponding to number of paths.",
                                         level=1,
                                         group_label="Raw paths",
                                         warning_text=" "))
visualize_section.add_entry(DefaultEntry(config_name="paths_smooth_io",
                                         name="Paths smooth io: ",
                                         default_values=[True],
                                         help_text="If set True arrows indicating beginning and end of paths, oriented accordingly to smooth paths, are displayed as separate objects or as one object with states corresponding to number of paths.",
                                         level=1,
                                         group_label="Smooth paths"))
visualize_section.add_entry(DefaultEntry(config_name="paths_states",
                                         name="Paths states: ",
                                         default_values=[True],
                                         help_text="If True objects displayed by paths_raw, paths_smooth, paths_raw_io, and paths_smooth_io are displayed as one object with states corresponding to number of paths. Otherwise they are displayed as separate objects.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="ctypes_raw",
                                         name="Ctypes raw: ",
                                         default_values=[True],
                                         help_text="Displays raw paths in a similar manner as non split all_paths_raw but each cluster type is displayed in separate object.",
                                         level=2,
                                         group_label="Raw paths"))
visualize_section.add_entry(DefaultEntry(config_name="ctypes_smooth",
                                         name="Ctypes smooth: ",
                                         default_values=[True],
                                         help_text="Displays smooth paths in a similar manner as non split all_paths_smooth but each cluster type is displayed in separate object.",
                                         level=2,
                                         group_label="Smooth paths"))
visualize_section.add_entry(DefaultEntry(config_name="ctypes_amount",
                                         name="Ctypes limit: ",
                                         default_values=[float()],
                                         help_text="Allows to limit number of visualised ctypes. If it is a number in range (0,1), then it is interpreted as percent number of ctypes to be visualized. It is is a integer number >= 1, it is total number of visualized ctypes.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="split_by_type",
                                         name="Split by type: ",
                                         default_values=[False],
                                         help_text="If True, visualizations of all objects that correspond to particular types of traced molecules are splited. Appropriate molecule name is added to created objects.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="retain_all_types",
                                         name="Retain all types: ",
                                         default_values=[False],
                                         help_text="If True and *split_by_type* option is True as well, standard visualization for all types of traced moelcules will be generated as well.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="inlets_clusters",
                                         name="Visualize cluster of inlets: ",
                                         default_values=[True],
                                         help_text="If set True, clusters of inlets are visualized.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="inlets_clusters_amount",
                                         name="Inlets limit: ",
                                         default_values=[str()],
                                         help_text="Allows to limit number of visualised inlets. If it is a number in range (0,1) then it is interpreted as percent number of inlets to be visualized. It is is a integer number >= 1 it is total number of visualized inlets.",
                                         level=0,
                                         warning_text=" "))
visualize_section.add_entry(DefaultEntry(config_name="show_molecule",
                                         name="Show molecule: ",
                                         default_values=[["protein", "False"]],
                                         help_text="If is set to selection of some molecular object in the system, for example to protein, this object is displayed.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="show_molecule_frames",
                                         name="Molecule frames: ",
                                         default_values=[0],
                                         help_text="Allows to indicate which frames of object defined by show_molecule should be displayed. It is possible to set several frames. In that case frames would be displayed as states.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="show_scope_chull",
                                         name="Show scope convex hull: ",
                                         default_values=[["protein", str()]],
                                         help_text="If is set to selection of some molecular object in the system, for example to protein, convex hull of this object is displayed.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="show_scope_chull_inflate",
                                         name="Show scope convex hull inflate: ",
                                         default_values=[str()],
                                         help_text="Increase (or if negative decrease) size of the scope convex hull.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="show_scope_chull_frames",
                                         name="Scope convex hull frames: ",
                                         default_values=[0],
                                         help_text="Allows to indicate for which frames of object defined by \"Show scope convex hull\" should be displayed. It is possible to set several frames. In that case frames would be displayed as states.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="show_object_chull",
                                         name="Show object convex hull: ",
                                         default_values=[str()],
                                         help_text="If is set to selection of some molecular object in the system convex hull of this object is displayed. This works exacly the same way as show_chull but is meant to mark object shape. It can be achieved by using name * and molecular object definition plus some spatial constrains, for example those used in object definition.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="show_object_chull_frames",
                                         name="Object convex hull frames:",
                                         default_values=[0],
                                         help_text="Allows to indicate for which frames of object defined by show_object convex hull should be displayed. It is possible to set several frames. In that case frames would be displayed as states.",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="cluster_area",
                                         name="Cluster area: ",
                                         default_values=[False],
                                         help_text="If set True, clusters’ areas are estimated with kernel density estimation method (KDE) and plotted as countour.",
                                         level=2))
visualize_section.add_entry(DefaultEntry(config_name="cluster_area_precision",
                                         name="Cluster area precision: ",
                                         default_values=[20],
                                         help_text="Precision of KDE method in clusters’ areas estimation method. This options controls number of grid points per one square A as used in KDE. Higher values means better precision. Number of points can be calculated as P^(2/3).",
                                         level=1))
visualize_section.add_entry(DefaultEntry(config_name="cluster_area_expand",
                                         name="Expand cluster area: ",
                                         default_values=[2],
                                         help_text="Space occupied by clusters’ points can be expanded before KDE calculation. This option controls amount of A by which the cluster space is expanded. Average amount of expansion can be calcualted as E^(2/3).",
                                         level=1))
DEFAULTS.append(visualize_section)

smooth_section = DefaultSection(config_name="smooth", name="Smooth", level=1)
smooth_section.add_entry(DefaultEntry(config_name="method",
                                      name="Method: ",
                                      default_values=[("window", "mss", "window_mss", "awin", "awin_mss", "dwin",
                                                       "dwin_mss", "savgol")],
                                      help_text="Smoothing method.",
                                      level=0))
smooth_section.add_entry(DefaultEntry(config_name="recursive",
                                      name="Recursive runs: ",
                                      default_values=[int()],
                                      help_text="Number of recursive runs of smoothing method.",
                                      level=0))
smooth_section.add_entry(DefaultEntry(config_name="window",
                                      name="Window size: ",
                                      default_values=[float()],
                                      help_text="In window based method defines window size. In plain window it has to be int number. In savgol it has to be odd integer.",
                                      level=0))
smooth_section.add_entry(DefaultEntry(config_name="step",
                                      name="Step: ",
                                      default_values=[float()],
                                      help_text="In step based method defines size of the step.",
                                      level=0))
smooth_section.add_entry(DefaultEntry(config_name="function",
                                      name="Averaging function: ",
                                      default_values=[("mean", "median")],
                                      help_text="In window based methods defines averaging function.",
                                      level=0))
smooth_section.add_entry(DefaultEntry(config_name="polyorder",
                                      name="Polyorder: ",
                                      default_values=[int()],
                                      help_text="In savgol is polynomial order.",
                                      level=0))

separate_paths_section.add_entry(smooth_section)

# Added rest of easy clustering section entries
clust_easy_section.add_entry(copy.copy(clustering_section))

clust_easy_section.add_entry(DefaultEntry(config_name="join_clusters",
                                          name="Join clusters: ",
                                          default_values=[str()],
                                          help_text="This option allows to join selected clusters. Clusters’ IDs joined with + character lists clusters to be joined together. Several such blocks separated by space can be used. For example, if set to 1+3+4 5+6 clusters 1, 3, and 4 will be joined in one cluster and cluster 5, and 6 will be also joined in another cluster.",
                                          level=2,
                                          group_label="Post Clustering Optimization"))
clust_easy_section.add_entry(DefaultEntry(config_name="renumber_clusters",
                                          name="Renumber clusters: ",
                                          default_values=[False],
                                          help_text="If set True, clusters have consecutive numbers starting from 1 (or 0 if outliers are present) starting from the bigest cluster.",
                                          level=2,
                                          group_label="Post Clustering Optimization"))

VALVE_DEFAULTS = DefaultSection("", "", 0)
VALVE_DEFAULTS.add_entry(DefaultEntry(config_name="-c",
                                      name="Config absolute file path: ",
                                      default_values=[filetype()],
                                      help_text="Absolute path to config file.",
                                      level=None
                                      ))
VALVE_DEFAULTS.add_entry(DefaultEntry(config_name="-t",
                                      name="Max threads: ",
                                      default_values=[str()],
                                      help_text="Limit Aqua-Duct calculations to given number of threads.",
                                      level=None
                                      ))
# VALVE_DEFAULTS.add_entry(DefaultEntry(config_name="--force-save",
#                                       name="Force saving results: ",
#                                       default_values=[False],
#                                       help_text="Force saving results.",
#                                       level=None
#                                       ))
VALVE_DEFAULTS.add_entry(DefaultEntry(config_name="--debug",
                                      name="Debug mode: ",
                                      default_values=[False],
                                      help_text="Prints debug info.",
                                      level=None
                                      ))
VALVE_DEFAULTS.add_entry(DefaultEntry(config_name="--debug-file",
                                      name="Debug absolute file path: ",
                                      default_values=[filetype()],
                                      help_text="Absolute path to debug file.",
                                      level=None
                                      ))

POND_DEFAULTS = DefaultSection("", "", 0)
POND_DEFAULTS.add_entry(DefaultEntry(config_name="-c",
                                     name="Config filename: ",
                                     default_values=[str()],
                                     help_text="Config file filename.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="-t",
                                     name="Max threads: ",
                                     default_values=[str()],
                                     help_text="Limit Aqua-Duct calculations to given number of threads.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="-r",
                                     name="Results directory: ",
                                     default_values=[str()],
                                     help_text="Path to results directory.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--paths-types",
                                     name="Paths types: ",
                                     default_values=[str()],
                                     help_text="Limit calculations to given paths types, i.e. given molecules.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--raw",
                                     name="Use raw paths: ",
                                     default_values=[str()],
                                     help_text="Use raw data from paths instead of single paths.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--raw-master",
                                     name="Use raw master paths: ",
                                     default_values=[False],
                                     help_text="Use raw data from paths instead of single paths, only in master paths calculations.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--raw-discard-singletons",
                                     name="Discard raw singletons: ",
                                     default_values=[str()],
                                     help_text="Discard short scope only segments from raw data.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--window-full",
                                     name="Full window: ",
                                     default_values=[False],
                                     help_text="Return full window if windows is used.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--windows",
                                     name="Windows: ",
                                     default_values=[str()],
                                     help_text="Number of windows to calculate.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--wsize",
                                     name="Window size: ",
                                     default_values=[str()],
                                     help_text="Size of window in frames.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--reference",
                                     name="Reference: ",
                                     default_values=[str()],
                                     help_text="Selection of reference in the first frame of trajectory.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--reference-radius",
                                     name="Reference radius: ",
                                     default_values=[str()],
                                     help_text="Radius of reference.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--reference-mol",
                                     name="Molecule reference: ",
                                     default_values=[str()],
                                     help_text="Selection of reference molecules.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--temperature",
                                     name="Temperature: ",
                                     default_values=[str()],
                                     help_text="Simulation temperature.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--gsize",
                                     name="Grid cells size: ",
                                     default_values=[str()],
                                     help_text="Size of grid's cells.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--pockets",
                                     name="Pockets: ",
                                     default_values=[False],
                                     help_text="Calculate pockets.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--hotspots",
                                     name="Hotspots: ",
                                     default_values=[False],
                                     help_text="Calculates hotspots if pockets are calculated.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--master-radius",
                                     name="Master paths radius: ",
                                     default_values=[False],
                                     help_text="Calculate profiles for master paths with given radius.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--master-ctypes",
                                     name="Ctypes: ",
                                     default_values=[False],
                                     help_text="Limit calculations to given ctypes.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--debug",
                                     name="Debug mode: ",
                                     default_values=[False],
                                     help_text="Prints debug info.",
                                     level=None
                                     ))
POND_DEFAULTS.add_entry(DefaultEntry(config_name="--debug-file",
                                     name="Debug mode: ",
                                     default_values=[False],
                                     help_text="Debug log file.",
                                     level=None
                                     ))
# @formatter:on

MENUS = [
    # section:entry
    "clustering:method",
    "reclustering:method"
]

LEVELS = {
    "Easy": 2,
    "Normal": 1,
    "Expert": 0
}
