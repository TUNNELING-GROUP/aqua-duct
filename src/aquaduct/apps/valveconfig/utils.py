# -*- coding: utf-8 -*-

import Tkinter as tk
import os
import ttk
from tkFileDialog import askopenfile


class longstr(str):
    """
    Class used to specify type of default value.

    Represents Text field.
    """
    pass


class filetype(object):
    """
    Class used to specify type of default value.

    Represents Entry with file loading button.
    """
    pass


def get_widget_bg(widget):
    """
    Return background color of specified widget.

    :param widget: Ttk widget.
    :return: Background color.
    :rtype: str
    """
    try:
        return ttk.Style().lookup(widget["style"], "background")
    except tk.TclError:
        # FIXME: Appending T may cause errors
        return ttk.Style().lookup("T" + widget.winfo_class(), "background")


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


def get_default_section(section_name):
    """
    Return default informations.

    :param section_name: Name of section which informations are demaned.
    :return: Default section informations.
    :rtype: dict
    """
    for section in DEFAULTS:
        if section_name == section["config_name"]:
            return section


def get_default_entry(section, option):
    """
    Return default entries.

    :param section: Name of section where option is located.
    :param option: Option name which default values are demanded.
    :return: tuple
    """
    for entry in get_default_section(section)["entries"]:
        if option == entry[0]:
            return entry


def widget_factory(parent, default, state=tk.NORMAL):
    """
    Creates widget depending on default argument.

    :param parent: Parent of new widget.
    :param default: Default widget value.
    :param state: State of widget.
    :return: Widget and variable attached to it.
    :rtype: tuple
    """
    if isinstance(default, longstr):  # Due to inheritance from str, longstr must be checked first
        v = tk.StringVar()
        v.set(default)

        w = Text(parent, textvariable=v, wrap=tk.NONE, state=state, width=34, height=5)

    elif isinstance(default, str):
        v = tk.StringVar()
        v.set(default)

        w = ttk.Entry(parent, textvariable=v, width=30, state=state, background=get_widget_bg(parent))
    elif isinstance(default, bool):
        v = tk.BooleanVar()
        v.set(default)

        w = ttk.Checkbutton(parent, variable=v, offvalue=False, onvalue=True, state=state)
    elif isinstance(default, int):
        v = tk.IntVar()
        v.set(default)

        w = ttk.Entry(parent, textvariable=v, width=5, state=state, background=get_widget_bg(parent))
    elif isinstance(default, float):
        v = tk.DoubleVar()
        v.set(default)

        w = ttk.Entry(parent, textvariable=v, width=5, state=state, background=get_widget_bg(parent))
    elif isinstance(default, tuple):
        v = tk.StringVar()

        w = ttk.OptionMenu(parent, v, default[0], *default)
    elif isinstance(default, list):
        raise TypeError("Use tuple for option menu, not list")
    else:
        raise TypeError("There is no specified behaviour for {} type".format(type(default)))

    return w, v


def entry_factory(parent, row, entry_name_long, default, help, state=tk.NORMAL):
    """
    Determines which class is used to handle specified default value.

    :param parent: Parent of widget.
    :param row: Row number where first Entry will be grided.
    :param entry_name_long: Readable entry name.
    :param default: Default values of entry.
    :param help: Text which will be displayed in tooltip.
    :param state: State of widget.
    :return: Entry based on default value.
    """
    if isinstance(default, list):
        if len(default) > 2:
            raise RuntimeError("There can be only two values in config defaults({})".format(entry_name_long))

        input_default, control_default = default

        if isinstance(control_default, bool):
            return BoolEntry(parent, row, entry_name_long, input_default, control_default, help)
        elif isinstance(control_default, filetype):
            if isinstance(input_default, longstr):
                return ManyFileEntry(parent, row, entry_name_long, input_default, help)
            elif isinstance(input_default, str):
                return FileEntry(parent, row, entry_name_long, input_default, help)
            else:
                raise TypeError(
                    "File can be loaded only into str or longstr widget type in {} option({})".format(entry_name_long,
                                                                                                      type(default)))
        elif isinstance(control_default, float):
            return ParenthesedEntry(parent, row, entry_name_long, input_default, control_default, help)
        else:
            raise TypeError("There is no specified behaviour for {} type(for {} option). "
                            "First must be input widget, then control widget"
                            .format(type(control_default), entry_name_long))
    else:
        return StandardEntry(parent, row, entry_name_long, default, help, state)


class Text(tk.Text, object):
    def __init__(self, parent, textvariable, **kwargs):
        """
         Text widget with ability to assign content to variable.

        :param parent: Parent of widget.
        :param textvariable: String variable to which will contain Text content.
        :param kwargs: Arguments which will be passed to original tk.Text widget.
        """
        super(Text, self).__init__(parent, **kwargs)

        self.var = textvariable

        # Block performing callback recursivelly
        self.block = False

        if self.var:
            # Every written char assign content to var
            self.bind("<Key>", self._on_update)

        if self.var:
            # Change Text content when variable will change
            self.var.trace("w", self._on_var_update)

    def _on_update(self, e):
        """
        Text widget callback which assign its content to variable.

        :param e: Event informations.
        """
        if not self.block:
            self.block = True

            self.var.set(self.get("1.0", tk.END))

            self.block = False

    def _on_var_update(self, *args):
        """
        Variable callback which assign its content to Text widget.

        :param args: Event informations.
        """
        if not self.block:
            self.block = True

            self.delete('1.0', tk.END)
            self.insert("1.0", self.var.get())

            self.block = False


class Entry(object):
    """
    Abstract class for various Entries that manage different type of default values.

    Represents single row of configuration option with Label and all input widgets.
    """
    def __init__(self):
        self.input_var = None
        self.control_var = None

        self.frame_pady = 5
        self.frame_sticky = "w"

        self.label_sticky = "E"

    def get(self):
        """
        Gets Entry value.

        :return: Entry value.
        """
        raise NotImplementedError()

    def set(self, value):
        """
        Sets Entry value.

        :param value: New value of Entry.
        """
        raise NotImplementedError()


class StandardEntry(Entry):
    def __init__(self, parent, row, entry_name_long, default, help, state):
        """
        Entry with standard widget.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(StandardEntry, self).__init__()

        ttk.Label(parent, text=entry_name_long, background=get_widget_bg(parent)).grid(sticky=self.label_sticky,
                                                                                       row=row, column=0)

        input_frame = tk.Frame(parent)
        input_frame.grid(row=row, column=1, sticky=self.frame_sticky, pady=self.frame_pady, padx=7)

        widget, self.input_var = widget_factory(input_frame, default, state)
        widget.pack()

        ToolTip.create(widget, help)

    def get(self):
        """
        Gets Entry value.

        :return: Entry value.
        """
        return self.input_var.get()

    def set(self, value):
        """
        Sets Entry value.

        :param value: New value of Entry.
        """
        self.input_var.set(value)


class BoolEntry(Entry):
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help):
        """
        Entry with Checkbox and Entry or text widget.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(BoolEntry, self).__init__()

        self.entry_name_long = entry_name_long

        ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_frame = ttk.Frame(parent)
        input_frame.grid(row=row, column=1, pady=self.frame_pady, sticky=self.frame_sticky)

        input_widget, self.input_var = widget_factory(input_frame, input_default)
        input_widget.pack(side=tk.RIGHT)

        ToolTip.create(input_widget, help)

        control_widget, self.control_var = widget_factory(input_frame, control_default)
        control_widget.pack(side=tk.LEFT, padx=6)

        ToolTip.create(control_widget, help)

    def get(self):
        """
        Get Entry value.

        :return: If Checkbox is checked it return input widget value, otherwise False.
        """
        if self.control_var.get():
            return self.input_var.get()
        else:
            return False

    def set(self, value):
        """
        Set Entry value.

        :param value: Value that input widget will be set too.
        If False it will set only Checkbox value to False and input widget to "".
        If value have str type it will set Checkbox to True and input widget to that value.
        """
        if value == "False" or value == False:  # If value = 0 expression should be true
            self.control_var.set(False)
            self.input_var.set("")
        else:
            self.control_var.set(True)
            self.input_var.set(value)


class FileEntry(Entry):
    def __init__(self, parent, row, entry_name_long, default, help):
        """
        Entry with Entry widget and button to load and append file name to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(FileEntry, self).__init__()

        ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_frame = ttk.Frame(parent)
        input_frame.grid(row=row, column=1, sticky=self.frame_sticky, pady=self.frame_pady)

        self.input_widget, self.input_var = widget_factory(input_frame, default)
        self.input_widget.pack(side=tk.LEFT, padx=7)

        ToolTip.create(self.input_widget, help)

        load_file_button = ttk.Button(input_frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.RIGHT)

        load_file_button.bind("<Button-1>", self.callback_load_file)

    def callback_load_file(self, e):
        """
        Callback for selecting file.

        Sets widget content to loaded file name.
        """
        try:
            with askopenfile("r") as f:
                self.input_var.set(f.name)
        except AttributeError:  # In case of cancel selecting file
            pass

    def get(self):
        """
        Gets Entry value.

        :return: Entry value.
        """
        return self.input_var.get()

    def set(self, value):
        """
        Sets Entry value.

        :param value: New value of Entry.
        """
        self.input_var.set(value)


class ManyFileEntry(Entry):
    def __init__(self, parent, row, entry_name_long, default, help):
        """
        Entry with Text widget and button to load and append file names to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(ManyFileEntry, self).__init__()

        ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_frame = ttk.Frame(parent)
        input_frame.grid(row=row, column=1, sticky=self.frame_sticky, pady=self.frame_pady)

        self.input_widget, self.input_var = widget_factory(input_frame, default)
        self.input_widget.pack(side=tk.LEFT, padx=7)

        ToolTip.create(self.input_widget, help)

        load_file_button = ttk.Button(input_frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.RIGHT)

        load_file_button.bind("<Button-1>", self.callback_load_file)

    def callback_load_file(self, e):
        """
        Callback for selecting file.

        Appends loaded file name at the end of Text widget.
        """
        try:
            with askopenfile("r") as f:
                self.input_widget.mark_set("insert", tk.END)
                self.input_var.set(self.input_var.get() + f.name + "\n")
        except AttributeError:  # In case of cancel selecting file
            pass

    def get(self):
        """
        Gets Entry value.

        :return: Entry value.
        """
        # TODO: validate path separator
        value = self.input_var.get()
        if self.input_var.get().endswith("\n"):
            value = value[:-1]

        return value.replace("\n", os.pathsep)

    def set(self, value):
        """
        Sets Entry value.

        :param value: New value. It can be single path or paths separated by os.pathsep.
        """
        self.input_var.set(value.replace(os.pathsep, "\n"))


class ParenthesedEntry(Entry):
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help):
        """
        Entry with Text widget and button to load and append file names to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(ParenthesedEntry, self).__init__()

        ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_frame = ttk.Frame(parent)
        input_frame.grid(row=row, column=1, pady=self.frame_pady, sticky=self.frame_sticky)

        input_widget, self.input_var = widget_factory(input_frame, input_default)
        input_widget.pack(side=tk.LEFT, padx=6)

        ToolTip.create(input_widget, help)

        control_widget, self.control_var = widget_factory(input_frame, control_default)
        control_widget.pack(side=tk.RIGHT, padx=6)

        ToolTip.create(control_widget, help)

    def get(self):
        """
        Gets Entry value.

        :return: Value of Entry widget and value of second Entry in parentheses, eg. Value1(Value2).
        """
        formatter_string = "{}"
        if self.control_var.get():
            formatter_string = "{}({})"

        return formatter_string.format(self.input_var.get(), self.control_var.get())

    def set(self, value):
        """
        Sets Entry value.

        :param value: First value and second value in parentheses or without second value and parentheses.
        """
        try:
            b_pos1 = value.index("(")
            b_pos2 = value.index(")")
            self.input_var.set(value[0:b_pos1])
            self.control_var.set(value[b_pos1 + 1:b_pos2])
        except ValueError:  # When there is no threshold
            self.input_var.set(value)
            self.control_var.set(0.0)


class HidingFrame(ttk.Frame, object):
    def __init__(self, parent, row, text, **kwargs):
        """
        Frame that remembers inner row for griding new widgets.
        Used to keep methods that depends on option menu value.

        :param parent: Parent of widget.
        :param row: Row where widgets will be grided.
        :param text: Title of the Frame.
        :param kwargs: Arguments which will be passed to original ttk.Frame widget.
        """
        super(HidingFrame, self).__init__(parent, **kwargs)

        ttk.Label(self, text=text, style="HF.TFrame.Label").grid(
            sticky="W", row=0, column=0, columnspan=2, padx=5, pady=5)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        self.row = row
        self.inner_row = 1

    def show(self):
        """ Method to grid Frame with predefinied configuration. """
        self.grid(sticky="EW", row=self.row, column=0, columnspan=2, padx=10, pady=10, ipady=5)


class CallbackWrapper(object):
    def __init__(self, callback, *args):
        """ Allow to use callbacks with predefined list of arguments. """
        self.callback = callback
        self.args = args

    def __call__(self, *args, **kwargs):
        self.callback(*self.args)


class ToolTip(object):
    """
    Code found: http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml
    """

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    @staticmethod
    def create(widget, text):
        toolTip = ToolTip(widget)

        def enter(event):
            toolTip.showtip(text)

        def leave(event):
            toolTip.hidetip()

        widget.bind('<Enter>', enter)
        widget.bind('<Leave>', leave)

    def showtip(self, text):
        # Display text in tooltip window
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 15
        y = y + cy + self.widget.winfo_rooty() + 15
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        try:
            # For Mac OS
            tw.tk.call("::tk::unsupported::MacWindowStyle",
                       "style", tw._w,
                       "help", "noActivates")
        except tk.TclError:
            pass
        label = ttk.Label(tw, text=self.text, justify=tk.LEFT,
                          background="#ffffe0", relief=tk.SOLID, borderwidth=1,
                          font=("TkDefaultFont", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


# http://tkinter.unpythonic.net/wiki/VerticalScrolledFrame
class VerticalScrolledFrame(tk.Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame
    * Construct and pack/place/grid normally
    * This frame only allows vertical scrolling
    """

    def __init__(self, parent, *args, **kw):
        tk.Frame.__init__(self, parent, *args, **kw)

        # create a canvas object and a vertical scrollbar for scrolling it
        vscrollbar = tk.Scrollbar(self, orient=tk.VERTICAL)
        vscrollbar.pack(fill=tk.Y, side=tk.RIGHT, expand=tk.FALSE)
        canvas = tk.Canvas(self, bd=0, highlightthickness=0,
                           yscrollcommand=vscrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.TRUE)
        vscrollbar.config(command=canvas.yview)

        # reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = tk.Frame(canvas)
        interior_id = canvas.create_window(0, 0, window=interior,
                                           anchor="nw")

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(event):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                canvas.config(width=interior.winfo_reqwidth())

        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())

        canvas.bind('<Configure>', _configure_canvas)


# @formatter:off
DEFAULTS = [
    {
        "config_name": "global",
        "name": "General",
        "name_long": "General options",
        "entries": [
            # Name used directly in code
            ("top", "Topology file: ", [str(), filetype()],
             "Path to topology file.\nAqua-Duct supports PDB, PRMTOP, PFS topology files.", None, 1),
            # Name used directly in code
            ("trj", "Trajectory file: ", [longstr(), filetype()],
             "Path to trajectory file.\nAqua-Duct supports NC and DCD trajectory files.", None, 1)
        ]
    },
    {
        "config_name": "traceable_residues",
        "name": "Stage I",
        "name_long": "Traceable residues",
        "entries": [
            ("execute", "Execute: ", ("runonce", "run", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run calculations are always performed and if dump is set dump file is saved.\nIf it is set to runonce calculations are performed if there is no dump file specified by dump option.\nIf it is present calculations are skipped and data is loaded from the file,\nIf it is set to skip calculations are skip and if dump is set data is loaded from the file.",
             None, 0),
            ("dump", "Dump file: ", "1_traceable_residues_data.dump",
             "File name of dump data.\nIt is used to save results of calculations or to load previously calculated data - this depends on execute option.",
             None, 0),
            ("scope", "Scope: ", str(), "Definition of Scope of interest.\nSee also Scope definition.", None, 1),
            ("scope_convexhull", "Convex hull scope: ", True,
             "Flag to set if Scope is direct or convex hull definition.", None, 1),
            ("scope_everyframe", "Everyframe scope", False,
             "Flag to set Scope evaluation mode.\nIf set True Scope is evaluated in every frame.\nThis make sense if the definition is complex and depends on distances between molecular entities.",
             None, 0),
            ("scope_convexhull_inflate", "scope_convexhull_inflate", "None", "", None, 0),
            ("object", "object", "None", "", None, 1),
            ("add_passing", "add_passing", "None", "", None, 0)
        ]
    },
    {
        "config_name": "raw_paths",
        "name": "Stage II",
        "name_long": "Raw paths",
        "entries": [
            ("execute", "Execute: ", ("runonce", "run", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run calculations are always performed and if dump is set dump file is saved.\nIf it is set to runonce calculations are performed if there is no dump file specified by dump option.\nIf it is present calculations are skipped and data is loaded from the file,\nIf it is set to skip calculations are skip and if dump is set data is loaded from the file.",
             None, 0),
            ("dump", "Dump file: ", "2_raw_paths_data.dump",
             "File name of dump data.\nIt is used to save results of calculations or to load previously calculated data - this depends on execute option.",
             None, 0),
            ("scope", "Scope: ", str(),
             "Definition of Scope of interest.\nSee also Scope definition.\nIf None value form previous stage is used.",
             None, 0),
            ("scope_convexhull", "Convex hull scope: ", "None",
             "Flag to set if the Scope is direct or convex hull definition.", None, 0),
            ("scope_everyframe", "Everyframe scope", False,
             "Flag to set Scope evaluation mode.\nIf set True Scope is evaluated in every frame.\nThis make sense if the definition is complex and depends on distances between molecular entities.\nIf None value from previous stage is used.",
             None, 0),
            ("scope_convexhull_inflate", "scope_convexhull_inflate", "None", "", None, 0),
            ("object", "Object: ", str(),
             "Definition of Object of interest.\nSee also Object definition.\nIf None value from the previous stage is used.",
             None, 0),
            ("clear_in_object_info", "Recalculate coś tam: ", False,
             "If it is set to True information on occupation of Object site by traceable residues calculated in the previous stage is cleared and have to be recalculated.\nThis is useful if definition of Object was changed.",
             None, 0),
            ("discard_singletons", "discard_singletons", 1, "", None, 0),
            ("discard_empty_paths", "discard_empty_paths", True, "", None, 0)
        ]
    },
    {
        "config_name": "separate_paths",
        "name": "Stage III",
        "name_long": "Separate paths",
        "entries": [
            ("execute", "Execute: ", ("runonce", "run", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run calculations are always performed and if dump is set dump file is saved.\nIf it is set to runonce calculations are performed if there is no dump file specified by dump option.\nIf it is present calculations are skipped and data is loaded from the file.\nIf it is set to skip calculations are skip and if dump is set data is loaded from the file.",
             None, 0),
            ("dump", "Dump file: ", "3_separate_paths_data.dump",
             "File name of dump data.\nIt is used to save results of calculations or to load previously calculated data - this depends on execute option.",
             None, 0),
            (
            "discard_empty_paths", "Discard empty paths: ", True, "If set to True empty paths are discarded.", None, 0),
            ("sort_by_id", "Sort by ID: ", True,
             "If set to True separate paths are sorted by ID.\nOtherwise they are sorted in order of appearance.", None,
             0),
            ("discard_short_paths", "Discard short paths: ", 20,
             "This option allows to discard paths which are shorter than the threshold which is defined as total number of frames.",
             None, 0),
            ("discard_short_object", "discard_short_object", 2.0,
             "This option allows to discard paths which objects are shorter than the threshold which is defined as total length in metric units.",
             None, 0),
            ("discard_short_logic", "discard_short_logic", "or",
             "If both discard_short_paths and discard_short_object options are used, this option allows to set combination logic.\nIf it is set or a path is discarded if any of discard criterion is met.\nIf it is set and both criteria have to be met to discard path.",
             None, 0),
            ("auto_barber", "auto_barber", str(),
             "This option allows to select molecular entity used in Auto Barber procedure.\nSee also Auto Barber and barber_with_spheres().",
             None, 1),
            ("auto_barber_mincut", "auto_barber_mincut", str(),
             "Minimal radius of spheres used in Auto Barber.\nIf a sphere has radius smaller then this value it is not used in AutoBarber procedure.\nThis option can be switched off by setting it to None.",
             None, 0),
            ("auto_barber_maxcut", "auto_barber_maxcut", 2.8,
             "Maximal radius of spheres used in Auto Barber.\nIf a sphere has radius greater then this value it is not used in AutoBarber procedure.\nThis option can be switched off by setting it to None.",
             None, 0),
            ("auto_barber_mincut_level", "auto_barber_mincut_level", True,
             "If set True spheres of radius smaller than mincut are resized to mincut value.", None, 0),
            ("auto_barber_maxcut_level", "auto_barber_maxcut_level", True,
             "If set True spheres of radius greater than maxcut are resized to maxcut value.", None, 0),
            ("auto_barber_tovdw", "auto_barber_tovdw", True,
             "Correct cutting sphere by decreasing its radius by VdW radius of the closest atom.", None, 0),
            ("allow_passing_paths", "allow_passing_paths", False,
             "If set True paths that do not enter the object are detected and added to the rest of paths as ‘passing’ paths.",
             None, 0),
        ]
    },
    {
        "config_name": "inlets_clusterization",
        "name": "Stage IV",
        "name_long": "Inlets clusterization",
        "entries": [

            ("execute", "Execute", ("runonce", "run", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run calculations are always performed and if dump is set dump file is saved.\nIf it is set to runonce calculations are performed if there is no dump file specified by dump option.\nIf it is present calculations are skipped and data is loaded from the file.\nIf it is set to skip calculations are skip and if dump is set data is loaded from the file.",
             None, 0),
            ("dump", "Dump file: ", "4_inlets_clusterization_data.dump",
             "File name of dump data.\nIt is used to save results of calculations or to load previously calculated data - this depends on execute option.",
             None, 0),
            ("recluster_outliers", "recluster_outliers", False,
             "If set to True reclusterization of outliers is executed according to the method defined in reclusterization section.",
             None, 0),
            ("detect_outliers", "detect_outliers", [str(), False],
             "If set, detection of outliers is executed.\nIt could be set as a floating point distance threshold or set to Auto.\nSee Clusterization of inlets for more details",
             None, 1),
            ("singletons_outliers", "singletons_outliers", ["False", False],
             "Maximal size of cluster to be considered as outliers.\nIf set to number > 0 clusters of that size are removed and their objects are moved to outliers.\nSee Clusterization of inlets for more details.",
             None, 1),
            # If 0 default clusterization section will be performed
            ("max_level", "max_level", 0, "Maximal number of recursive clusterization levels.", None, 0),
            ("create_master_paths", "create_master_paths", False,
             "If set to True master paths are created (fast CPU and big RAM recommended; 50k frames long simulation may need ca 20GB of memory)",
             None, 0),
            ("exclude_passing_in_clusterization", "exclude_passing_in_clusterization", True,
             "If set to True passing paths are not clustered with normal paths.", None, 0),
            ("add_passing_to_clusters", "add_passing_to_clusters", str(),
             "Allows to run procedure for adding passing paths inlets to clusters with Auto Barber method.\nTo enable this the option should be set to molecular entity that will be used by Auto Barber.",
             None, 0),
            ("renumber_clusters", "renumber_clusters", False, "", None, 1),
            ("join_clusters", "join_clusters", "None", "", None, 1)
        ]
    },
    {
        "config_name": "analysis",
        "name": "Stage V",
        "name_long": "Analysis",
        "entries": [
            ("execute", "Execute", ("run", "runonce", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run or runonce stage is executed and results is saved according to save option.\nIf it is set to skip stage is skipped.",
             None, 0),
            ("save", "Save file: ", "???", "File name for saving results.", None, 0),
            ("dump_config", "dump_config", True,
             "If set to True configuration options, as seen by Valve, are added to the head of results.", None, 0),
            ("calculate_scope_object_size", "calculate_scope_object_size", False,
             "If set to True volumes and areas of object and scope approximated by convex hulls will be calculated for each analyzed frames and saved in output CSV file.",
             None, 0),
            ("scope_chull", "scope_chull", str(), "Scope convex hull definition used in calculating volume and area.",
             None, 0),
            ("scope_chull_inflate", "scope_chull_inflate", str(),
             "Scope convex hull definition used in calculating volume and area.", None, 0),
            (
            "object_chull", "object_chull", str(), "Object convex hull definition used in calculating volume and area.",
            None, 0),
        ]
    },
    {
        "config_name": "visualize",
        "name": "Stage VI",
        "name_long": "Visualize",
        "entries": [
            ("execute", "Execute: ", ("run", "runonce", "skip"),
             "Option controls stage execution.\nIt can have one of three possible values: run, runonce, and skip.\nIf it is set to run or runonce stage is executed and results is saved according to save option.\nIf it is set to skip stage is skipped.",
             None, 0),
            ("save", "Save file: ", "???", "File name for saving results.", None, 0),
            ("all_paths_raw", "all_paths_raw", False,
             "If True produces one object in PyMOL that holds all paths visualized by raw coordinates.", None, 1),
            ("all_paths_smooth", "all_paths_smooth", False,
             "If True produces one object in PyMOL that holds all paths visualized by smooth coordinates.", None, 1),
            ("all_paths_split", "all_paths_split", False,
             "If is set True objects produced by all_paths_raw and all_paths_smooth are split into Incoming, Object, and Outgoing parts and visualized as three different objects.",
             None, 1),
            ("all_paths_raw_io", "all_paths_raw_io", False,
             "If set True arrows pointing beginning and end of paths are displayed oriented accordingly to raw paths orientation.",
             None, 1),
            ("all_paths_smooth_io", "all_paths_smooth_io", False,
             "If set True arrows pointing beginning and end of paths are displayed oriented accordingly to smooth paths orientation.",
             None, 1),
            ("all_paths_amount", "all_paths_amount", "None", "", None, 0),
            ("simply_smooths", "simply_smooths", [("RecursiveVector", "HobbitVector", "OneWayVector",
                                                   "RecursiveTriangle", "HobbitTriangle", "OneWayTriangle"), float()],
             "Option indicates linear simplification method to be used in plotting smooth paths.\nSimplification removes points which do not (or almost do not) change the shape of smooth path.\nPossible choices are: RecursiveVector (LinearizeRecursiveVector), HobbitVector (LinearizeHobbitVector), OneWayVector (LinearizeOneWayVector), RecursiveTriangle (LinearizeRecursiveTriangle), HobbitTriangle (LinearizeHobbitTriangle), OneWayTriangle (LinearizeOneWayTriangle).\nOptionally name of the method can be followed by a threshold value in parentheses, i.e.\nRecursiveVector(0.05).\nFor sane values of thresholds see appropriate documentation of each method.\nDefault values work well.\nThis option is not case sensitive.\nIt is recommended to use default method or HobbitVector method.",
             None, 0),
            ("paths_raw", "paths_raw", False,
             "If set True raw paths are displayed as separate objects or as one object with states corresponding to number of path.",
             None, 1),
            ("paths_smooth", "paths_smooth", False,
             "If set True smooth paths are displayed as separate objects or as one object with states corresponding to number of path.",
             None, 1),
            ("paths_raw_io", "paths_raw_io", False,
             "If set True arrows indicating beginning and end of paths, oriented accordingly to raw paths, are displayed as separate objects or as one object with states corresponding to number of paths.",
             None, 1),
            ("paths_smooth_io", "paths_smooth_io", False,
             "If set True arrows indicating beginning and end of paths, oriented accordingly to smooth paths, are displayed as separate objects or as one object with states corresponding to number of paths.",
             None, 1),
            ("paths_states", "paths_states", False,
             "If True objects displayed by paths_raw, paths_smooth, paths_raw_io, and paths_smooth_io are displayed as one object with states corresponding to number of paths.\nOtherwise they are displayed as separate objects.",
             None, 1),
            ("ctypes_raw", "ctypes_raw", False,
             "Displays raw paths in a similar manner as non split all_paths_raw but each cluster type is displayed in separate object.",
             None, 1),
            ("ctypes_smooth", "ctypes_smooth", False,
             "Displays smooth paths in a similar manner as non split all_paths_smooth but each cluster type is displayed in separate object.",
             None, 1),
            ("ctypes_amount", "ctypes_amount", "None", "", None, 1),
            ("inlets_clusters", "inlets_clusters", False, "", None, 1),
            ("inlets_clusters_amount", "inlets_clusters_amount", "None", "", None, 0),
            ("show_molecule", "show_molecule", [str(), False],
             "If is set to selection of some molecular object in the system, for example to protein, this object is displayed.",
             None, 1),
            ("show_molecule_frames", "show_molecule_frames", 0,
             "Allows to indicate which frames of object defined by show_molecule should be displayed.\nIt is possible to set several frames.\nIn that case frames would be displayed as states.",
             None, 1),
            ("show_scope_chull", "show_scope_chull", [str(), False],
             "If is set to selection of some molecular object in the system, for example to protein, convex hull of this object is displayed.",
             None, 1),
            ("show_scope_chull_inflate", "show_scope_chull_inflate", "None", "", None, 1),
            ("show_scope_chull_frames", "show_scope_chull_frames", 0,
             "Allows to indicate for which frames of object defined by show_chull convex hull should be displayed.\nIt is possible to set several frames.\nIn that case frames would be displayed as states.",
             None, 1),
            ("show_object_chull", "show_object_chull", [str(), False],
             "If is set to selection of some molecular object in the system convex hull of this object is displayed.\nThis works exacly the same way as show_chull but is meant to mark object shape.\nIt can be achieved by using name * and molecular object definition plus some spatial constrains, for example those used in object definition.",
             None, 1),
            ("show_object_chull_frames", "show_object_chull_frames", 0,
             "Allows to indicate for which frames of object defined by show_object convex hull should be displayed.\nIt is possible to set several frames.\nIn that case frames would be displayed as states.",
             None, 1),
        ]
    },
    {
        "config_name": "smooth",
        "name": "Smooth",
        "name_long": "Smooth",
        "entries": [
            ("method", "method", ("window", "mss", "window_mss", "awin", "awin_mss", "dwin", "dwin_mss", "savgol"),
             "Smoothing method.", None, 0),
            ("recursive", "recursive", "???", "Number of recursive runs of smoothing method.", None, 0),
            ("window", "window", "???",
             "In window based method defines window size.\nIn plain window it has to be int number.\nIn savgol it has to be odd integer.",
             None, 0),
            ("step", "step", "???", "In step based method defines size of the step.", None, 0),
            ("function", "function", "???",
             "In window based methods defines averaging function.\nCan be mean or median.", None, 0),
            ("polyorder", "polyorder", "???", "In savgol is polynomial order.", None, 0),
        ]
    },
    {
        "config_name": "clusterization",
        "name": "Clusterization",
        "name_long": "Clusterization",
        "entries": [
            ("name", "Name", "clusterization", "clusterization name", None, 1),
            ("method", "method", ("barber", "dbscan", "affprop", "meanshift", "birch", "kmeans"),
             "Name of clusterization method.\nIt has to be one of the following: barber, dbscan, affprop, meanshift, birch, kmeans.\nDefault value depends whether it is clusterization section (barber) or reclusterization section (dbscan).",
             None, 1),
            ("recursive_clusterization", "recursive_clusterization", "clusterization",
             "If it is set to name of some section that holds clusterization method settings this method will be called in the next recursion of clusteriation.\nDefault value for reclusterization is None.",
             None, 1),
            ("recursive_threshold", "recursive_threshold", str(),
             "Allows to set threshold that excludes clusters of certain size from reclusterization.\nValue of this option comprises of operator and value.\nOperator can be one of the following: >, >=, <=, <.\nValue have to be expressed as floating number and it have to be in the range of 0 to 1.\nOne can use several definitions separated by a space character.\nOnly clusters of size complying with all thresholds definitions are submitted to reclusterization.",
             None, 1),

            # This is exceptional tuple, 5th parameter indicates to which method it belong
            # so there is possibility to keep visible only that one which belong to selected method

            # Barber options
            ("auto_barber", "auto_barber", str(), str(), "barber", 1),
            ("auto_barber_mincut", "auto_barber_mincut", str(), "", "barber", 0),
            ("auto_barber_maxcut", "auto_barber_maxcut", 2.8, "", "barber", 0),
            ("auto_barber_mincut_level", "auto_barber_mincut_level", True, "", "barber", 0),
            ("auto_barber_maxcut_level", "auto_barber_maxcut_level", True, "", "barber", 0),
            ("auto_barber_tovdw", "auto_barber_tovdw", True, "", "barber", 0),

            # Dbscan options
            ("eps", "eps", float(), "", "dbscan", 0),
            ("min_samples", "min_samples", int(), "", "dbscan", 0),
            ("metric", "metric", ("euclidean", "cityblock", "cosine", "manhattan"), "", "dbscan", 0),
            ("algorithm", "algorithm", ("auto", "ball_tree", "kd_tree", "brute"), "", "dbscan", 0),
            ("leaf_size", "leaf_size", int(), "", "dbscan", 0),

            # Affprop options
            ("damping", "damping", float(), "", "affprop", 0),
            ("convergence_iter", "convergence_iter", int(), "", "affprop", 0),
            ("max_iter", "max_iter", int(), "", "affprop", 0),
            ("preference", "preference", float(), "", "affprop", 0),

            # Meanshift options
            ("bandwidth", "bandwidth", "Auto", "", "meanshift", 1),
            ("cluster_all", "cluster_all", bool(), "", "meanshift", 0),
            ("bin_seeding", "bin_seeding", bool(), "", "meanshift", 0),
            ("min_bin_freq", "min_bin_freq", int(), "", "meanshift", 0),

            # Bircz options
            ("threshold", "threshold", float(), "", "birch", 0),
            ("branching_factor", "branching_factor", int(), "", "birch", 0),
            ("n_clusters", "n_clusters", int(), "", "birch", 1),

            # Kmeans options
            ("n_clusters", "n_clusters", int(), "", "kmeans", 1),
            ("max_iter", "max_iter", int(), "", "kmeans", 0),
            ("n_init", "n_init", int(), "", "kmeans", 0),
            ("init", "init", str(), "", "kmeans", 0),
            ("tol", "tol", float(), "", "kmeans", 0),
        ]
    },
    {
        "config_name": "reclusterization",
        "name": "Reclusterization",
        "name_long": "Reclusterization",
        "entries": [
            ("name", "Name", "reclusterization", "clusterization name", None, 1),
            ("method", "method", ("barber", "dbscan", "affprop", "meanshift", "birch", "kmeans"),
             "Name of clusterization method.\nIt has to be one of the following: barber, dbscan, affprop, meanshift, birch, kmeans.\nDefault value depends whether it is clusterization section (barber) or reclusterization section (dbscan).",
             None, 1),
            ("recursive_clusterization", "recursive_clusterization", "clusterization",
             "If it is set to name of some section that holds clusterization method settings this method will be called in the next recursion of clusteriation.\nDefault value for reclusterization is None.",
             None, 1),
            ("recursive_threshold", "recursive_threshold", str(),
             "Allows to set threshold that excludes clusters of certain size from reclusterization.\nValue of this option comprises of operator and value.\nOperator can be one of the following: >, >=, <=, <.\nValue have to be expressed as floating number and it have to be in the range of 0 to 1.\nOne can use several definitions separated by a space character.\nOnly clusters of size complying with all thresholds definitions are submitted to reclusterization.",
             None, 1),

            # This is exceptional tuple, 5th parameter indicates to which method it belong
            # so there is possibility to keep visible only that one which belong to selected method

            # Barber options
            ("auto_barber", "auto_barber", str(), str(), "barber", 1),
            ("auto_barber_mincut", "auto_barber_mincut", str(), "", "barber", 0),
            ("auto_barber_maxcut", "auto_barber_maxcut", 2.8, "", "barber", 0),
            ("auto_barber_mincut_level", "auto_barber_mincut_level", True, "", "barber", 0),
            ("auto_barber_maxcut_level", "auto_barber_maxcut_level", True, "", "barber", 0),
            ("auto_barber_tovdw", "auto_barber_tovdw", True, "", "barber", 0),

            # Dbscan options
            ("eps", "eps", float(), "", "dbscan", 0),
            ("min_samples", "min_samples", int(), "", "dbscan", 0),
            ("metric", "metric", ("euclidean", "cityblock", "cosine", "manhattan"), "", "dbscan", 0),
            ("algorithm", "algorithm", ("auto", "ball_tree", "kd_tree", "brute"), "", "dbscan", 0),
            ("leaf_size", "leaf_size", int(), "", "dbscan", 0),

            # Affprop options
            ("damping", "damping", float(), "", "affprop", 0),
            ("convergence_iter", "convergence_iter", int(), "", "affprop", 0),
            ("max_iter", "max_iter", int(), "", "affprop", 0),
            ("preference", "preference", float(), "", "affprop", 0),

            # Meanshift options
            ("bandwidth", "bandwidth", "Auto", "", "meanshift", 1),
            ("cluster_all", "cluster_all", bool(), "", "meanshift", 0),
            ("bin_seeding", "bin_seeding", bool(), "", "meanshift", 0),
            ("min_bin_freq", "min_bin_freq", int(), "", "meanshift", 0),

            # Bircz options
            ("threshold", "threshold", float(), "", "birch", 0),
            ("branching_factor", "branching_factor", int(), "", "birch", 0),
            ("n_clusters", "n_clusters", int(), "", "birch", 1),

            # Kmeans options
            ("n_clusters", "n_clusters", int(), "", "kmeans", 1),
            ("max_iter", "max_iter", int(), "", "kmeans", 0),
            ("n_init", "n_init", int(), "", "kmeans", 0),
            ("init", "init", str(), "", "kmeans", 0),
            ("tol", "tol", float(), "", "kmeans", 0),
        ]
    },
]
# @formatter:on

# Set option menu as "menu"
# Selecting option of option menu will show only entries which 4th element in DEFAULTSs is the same as option name
MENUS = [
    # section:entry
    "clusterization:method",
    "reclusterization:method"
]

LEVELS = {
    "Easy": 1,
    "Normal": 0
}

LOGO_ENCODED = """
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
"""
