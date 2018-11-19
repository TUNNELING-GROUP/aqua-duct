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
import os
import ttk
from tkFileDialog import askopenfile

import defaults


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


def widget_factory(parent, default, state=tk.NORMAL):
    """
    Creates widget depending on default argument.

    :param parent: Parent of new widget.
    :param default: Default widget value.
    :param state: State of widget.
    :return: Widget and variable attached to it.
    :rtype: tuple
    """
    if isinstance(default, defaults.longstr):  # Due to inheritance from str, longstr must be checked first
        v = tk.StringVar()
        v.set(default)

        w = Text(parent, textvariable=v, wrap=tk.WORD, state=state, width=34, height=5)
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
        v = tk.StringVar()

        w = ttk.Combobox(parent, textvariable=v)
        w["values"] = default
        w.current(0)
    else:
        raise TypeError("There is no specified behaviour for {} type".format(type(default)))

    return w, v


def entry_factory(parent, row, entry_name, default, help, state=tk.NORMAL, info_text=None, warning_text=None):
    """
    Determines which class is used to handle specified default value.

    :param parent: Parent of widget.
    :param row: Row number where first Entry will be grided.
    :param entry_name: Readable entry name.
    :param default: Default values of entry.
    :param help: Text which will be displayed in tooltip.
    :param state: State of widget.
    :return: Entry based on default value.
    """
    if len(default) > 2:
        raise RuntimeError("There can be only two values in config defaults({})".format(entry_name))
    elif len(default) == 2:
        input_default, control_default = default

        if isinstance(control_default, bool):
            return BoolEntry(parent, row, entry_name, input_default, control_default, help, info_text, warning_text)
        elif isinstance(control_default, defaults.filetype):
            if isinstance(input_default, str):
                return FileEntry(parent, row, entry_name, input_default, help, info_text, warning_text)
            else:
                raise TypeError(
                    "File can be loaded only into str widget type in {} option({})".format(entry_name,
                                                                                           type(default)))
        elif isinstance(control_default, defaults.manyfiletype):
            if isinstance(input_default, str):
                return ManyFileEntry(parent, row, entry_name, input_default, help, info_text, warning_text)
            else:
                raise TypeError(
                    "File can be loaded only into str widget type in {} option({})".format(entry_name,
                                                                                           type(default)))
        elif isinstance(control_default, float):
            return ParenthesedEntry(parent, row, entry_name, input_default, control_default, help, info_text,
                                    warning_text)
        else:
            raise TypeError("There is no specified behaviour for {} type(for {} option). "
                            "First must be input widget, then control widget"
                            .format(type(control_default), entry_name))
    else:
        return StandardEntry(parent, row, entry_name, default[0], help, state, info_text, warning_text)


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
    def __init__(self, parent, row, entry_name_long, default, help, state, info_text=None, warning_text=None):
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
        widget.pack(side=tk.LEFT)

        if info_text:
            InfoIconWidget(input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help, info_text, warning_text):
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

        if info_text:
            InfoIconWidget(input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, default, help, info_text, warning_text):
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
        load_file_button.pack(side=tk.LEFT)

        load_file_button.bind("<Button-1>", self.callback_load_file)

        if info_text:
            InfoIconWidget(input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, default, help, info_text, warning_text):
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

        self.parent = parent
        self.row = row
        self.entry_name_long = entry_name_long
        self.default = default
        self.help = help

        self.input_vars = []
        self.frames = []  # Store single row of widgets

        ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        self.input_frame = ttk.Frame(parent)
        self.input_frame.grid(row=row, column=1, sticky=self.frame_sticky, pady=self.frame_pady)

        self.append_entry()

        if info_text:
            InfoIconWidget(self.input_frame, info_text).grid(row=0, column=2)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).grid(row=0, column=2)

    def append_entry(self):
        """ Creates new entry with input widget and load button """
        frame = ttk.Frame(self.input_frame)
        frame.pack()
        self.frames.append(frame)

        input_widget, input_var = widget_factory(frame, self.default)

        input_widget.pack(side=tk.LEFT, padx=7)

        self.input_vars.append(input_var)

        ToolTip.create(input_widget, self.help)

        load_file_button = ttk.Button(frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.LEFT)

        callback = CallbackWrapper(self.callback_load_file, self.input_vars.index(input_var))
        load_file_button.bind("<Button-1>", callback)

    def callback_load_file(self, index):
        """
        Callback for selecting file.

        Appends loaded file name at the end of Text widget.

        :param index: Index of variable in self.input_vars
        """
        try:
            with askopenfile("r") as f:
                self.input_vars[index].set(f.name)

                # Append new entry only if file is loaded to last entry
                if len(self.input_vars) - 1 <= index:
                    self.append_entry()
        except AttributeError:  # In case of cancel selecting file
            pass

    def get(self):
        """
        Gets Entry value.

        :return: Entry value.
        """
        # TODO: validate path separator
        return os.pathsep.join([var.get() for var in self.input_vars if var.get()])

    def set(self, value):
        """
        Sets Entry value.

        If value is set to "" it deletes all input widgets, except first and sets it to ""

        :param value: New value. It can be single path or paths separated by os.pathsep.
        """
        if value == "":
            self.input_vars[0].set("")

            for i in reversed(range(1, len(self.frames))):
                self.frames[i].pack_forget()
                del self.frames[i]
                del self.input_vars[i]

            return

        for i, path in enumerate(value.split(os.pathsep)):
            self.input_vars[i].set(path)
            self.append_entry()


class ParenthesedEntry(Entry):
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help, info_text, warning_text):
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
        control_widget.pack(side=tk.LEFT, padx=6)

        if info_text:
            InfoIconWidget(input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(input_frame, warning_text).pack(side=tk.LEFT)

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


class WarningIconWidget(ttk.Label, object):
    def __init__(self, parent, text):
        """
        Widget with waring icon and Tooltip information

        :param parent: Parent of widget.
        :param text: Content of tooltip
        """
        self.image = tk.PhotoImage(data=WARNING_ICON)
        super(WarningIconWidget, self).__init__(parent, image=self.image, padding=0)

        ToolTip.create(self, text)


class InfoIconWidget(ttk.Label, object):
    def __init__(self, parent, text):
        """
        Widget with info icon and Tooltip information

        :param parent: Parent of widget.
        :param text: Content of tooltip
        """
        self.image = tk.PhotoImage(data=INFO_ICON)
        super(InfoIconWidget, self).__init__(parent, image=self.image, padding=0)

        ToolTip.create(self, text)


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
                          font=("TkDefaultFont", "8", "normal"), wraplength=400)
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


INFO_ICON = """
iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gsRFS
M6amXeigAAABl0RVh0Q29tbWVudABDcmVhdGVkIHdpdGggR0lNUFeBDhcAAAAzSURBVDjLY9Sdd+8/AwWACZ/kpURFhkuJiuQbQAxgwSepN/8+ZV4Y
GmEwasCwMGDgkzIA61QLdLWHcd4AAAAASUVORK5CYII=
"""

WARNING_ICON = """
iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gsRFS
IvHqMLIAAAABl0RVh0Q29tbWVudABDcmVhdGVkIHdpdGggR0lNUFeBDhcAAAA+SURBVDjLY/z/n+E/AwWAiYFCgNcARkYIHjgXDA0DWPBJsrFRaICg
IOFYwGuAsDBhFzCOpkT8gfj/Px3CAABswAh9ETyJyQAAAABJRU5ErkJggg==
"""

LOGO_ENCODED = """
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
"""
