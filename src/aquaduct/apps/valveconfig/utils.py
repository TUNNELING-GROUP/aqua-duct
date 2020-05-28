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
import os
import tkinter.ttk
from tkinter.filedialog import askopenfile, askdirectory

from . import defaults
from aquaduct.apps.valveconfig import get_img


def get_widget_bg(widget):
    """
    Return background color of specified widget.

    :param widget: Ttk widget.
    :return: Background color.
    :rtype: str
    """
    try:
        return tkinter.ttk.Style().lookup(widget["style"], "background")
    except tk.TclError:
        # FIXME: Appending T may cause errors
        return tkinter.ttk.Style().lookup("T" + widget.winfo_class(), "background")


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

        w = tkinter.ttk.Entry(parent, textvariable=v, width=30, state=state, background=get_widget_bg(parent))
    elif isinstance(default, bool):
        v = tk.BooleanVar()
        v.set(default)

        w = tkinter.ttk.Checkbutton(parent, variable=v, offvalue=False, onvalue=True, state=state)
    elif isinstance(default, int):
        v = tk.IntVar()
        v.set(default)

        w = tkinter.ttk.Entry(parent, textvariable=v, width=5, state=state, background=get_widget_bg(parent))
    elif isinstance(default, float):
        v = tk.DoubleVar()
        v.set(default)

        w = tkinter.ttk.Entry(parent, textvariable=v, width=5, state=state, background=get_widget_bg(parent))
    elif isinstance(default, tuple):
        v = tk.StringVar()

        w = tkinter.ttk.OptionMenu(parent, v, default[0], *default)
    elif isinstance(default, list):
        v = tk.StringVar()

        w = tkinter.ttk.Combobox(parent, textvariable=v)
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
        elif isinstance(control_default, float):
            return ParenthesedEntry(parent, row, entry_name, input_default, control_default, help, info_text,
                                    warning_text)
        else:
            raise TypeError("There is no specified behaviour for {} type(for {} option). "
                            "First must be input widget, then control widget"
                            .format(type(control_default), entry_name))
    elif isinstance(default[0], defaults.filetype):
        return FileEntry(parent, row, entry_name, str(), help, info_text, warning_text)
    elif isinstance(default[0], defaults.manyfiletype):
        return ManyFileEntry(parent, row, entry_name, str(), help, info_text, warning_text)
    elif isinstance(default[0], defaults.dirtype):
        return DirEntry(parent, row, entry_name, default[0], help, info_text, warning_text)

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

    def __init__(self, parent, row):
        self.input_var = None
        self.control_var = None

        self.input_frame = tk.Frame(parent)
        self.input_frame.grid(row=row, column=1, sticky="w", pady=5)

        self.hightlight_color = "OrangeRed2"
        self.default_background = None

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

    def highlight(self):
        """
        Changes color of input frame.

        Used to highlight, which required entry is unfilled.
        """
        self.default_background = self.input_frame.cget("background")
        self.input_frame.config(background=self.hightlight_color)

    def unhighlight(self):
        """
        Sets entry to default color.
        """
        self.input_frame.configure(bg=self.default_background)


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
        super(StandardEntry, self).__init__(parent, row)

        tkinter.ttk.Label(parent, text=entry_name_long, background=get_widget_bg(parent)).grid(sticky=self.label_sticky,
                                                                                       row=row, column=0)

        widget, self.input_var = widget_factory(self.input_frame, default, state)
        widget.pack(side=tk.LEFT, padx=5, pady=5)

        if info_text:
            InfoIconWidget(self.input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help, info_text=None,
                 warning_text=None):
        """
        Entry with Checkbox and Entry or text widget.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(BoolEntry, self).__init__(parent, row)

        self.entry_name_long = entry_name_long

        tkinter.ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_widget, self.input_var = widget_factory(self.input_frame, input_default)
        input_widget.pack(side=tk.RIGHT)

        ToolTip.create(input_widget, help)

        control_widget, self.control_var = widget_factory(self.input_frame, control_default)
        control_widget.pack(side=tk.LEFT)

        if info_text:
            InfoIconWidget(self.input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, default, help, info_text=None, warning_text=None):
        """
        Entry with Entry widget and button to load and append file name to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(FileEntry, self).__init__(parent, row)

        tkinter.ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        self.input_widget, self.input_var = widget_factory(self.input_frame, default)
        self.input_widget.pack(side=tk.LEFT, padx=5, pady=5)

        ToolTip.create(self.input_widget, help)

        load_file_button = tkinter.ttk.Button(self.input_frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.LEFT, padx=5)

        load_file_button.bind("<Button-1>", self.callback_load_file)

        if info_text:
            InfoIconWidget(self.input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).pack(side=tk.LEFT)

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
    def __init__(self, parent, row, entry_name_long, default, help, info_text=None, warning_text=None):
        """
        Entry with Text widget and button to load and append file names to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(ManyFileEntry, self).__init__(parent, row)

        self.parent = parent
        self.row = row
        self.entry_name_long = entry_name_long
        self.default = default
        self.help = help

        self.input_vars = []
        self.frames = []  # Store single row of widgets

        tkinter.ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        self.append_entry()

        if info_text:
            InfoIconWidget(self.input_frame, info_text).grid(row=0, column=2)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).grid(row=0, column=2)

    def append_entry(self):
        """ Creates new entry with input widget and load button """
        frame = tk.Frame(self.input_frame)
        frame.pack(anchor="nw", ipady=5)
        self.frames.append(frame)

        input_widget, input_var = widget_factory(frame, self.default)
        input_widget.pack(side=tk.LEFT, padx=5)

        self.input_vars.append(input_var)

        ToolTip.create(input_widget, self.help)

        load_file_button = tkinter.ttk.Button(frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.LEFT, padx=5)

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

            for i in reversed(list(range(1, len(self.frames)))):
                self.frames[i].pack_forget()
                del self.frames[i]
                del self.input_vars[i]

            return

        for i, path in enumerate(value.split(os.pathsep)):
            self.input_vars[i].set(path)
            self.append_entry()

    def highlight(self):
        self.default_background = self.input_frame.cget("background")
        for frame in self.frames:
            frame.config(background=self.hightlight_color)

    def unhighlight(self):
        for frame in self.frames:
            frame.config(background=self.default_background)


class DirEntry(Entry):
    def __init__(self, parent, row, entry_name_long, default, help, info_text=None, warning_text=None):
        """
        Entry with Entry widget and button to load and append file name to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(DirEntry, self).__init__(parent, row)

        tkinter.ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        self.input_widget, self.input_var = widget_factory(self.input_frame, default)
        self.input_widget.pack(side=tk.LEFT, padx=5, pady=5)

        ToolTip.create(self.input_widget, help)

        load_file_button = tkinter.ttk.Button(self.input_frame, text="Load", style="File.TButton")
        load_file_button.pack(side=tk.LEFT, padx=5)

        load_file_button.bind("<Button-1>", self.callback_load_dir)

        if info_text:
            InfoIconWidget(self.input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).pack(side=tk.LEFT)

    def callback_load_dir(self, e):
        """
        Callback for selecting dir.

        Sets widget content to loaded dir name.
        """
        try:
            selected_dir = askdirectory()
            self.input_var.set(selected_dir)
        except AttributeError:  # In case of cancel selecting dir
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


class ParenthesedEntry(Entry):
    def __init__(self, parent, row, entry_name_long, input_default, control_default, help, info_text=None,
                 warning_text=None):
        """
        Entry with Text widget and button to load and append file names to it.

        :param parent: Parent of widgets.
        :param row: Row where widgets will be grided.
        :param entry_name_long: Readable entry name.
        :param default: Default values of entry.
        :param help: Text which will be displayed in tooltip.
        :param state: State of widget.
        """
        super(ParenthesedEntry, self).__init__(parent, row)

        tkinter.ttk.Label(parent, text=entry_name_long).grid(sticky=self.label_sticky, row=row, column=0)

        input_widget, self.input_var = widget_factory(self.input_frame, input_default)
        input_widget.pack(side=tk.LEFT)

        ToolTip.create(input_widget, help)

        control_widget, self.control_var = widget_factory(self.input_frame, control_default)
        control_widget.pack(side=tk.LEFT)

        if info_text:
            InfoIconWidget(self.input_frame, info_text).pack(side=tk.LEFT)
        elif warning_text:
            WarningIconWidget(self.input_frame, warning_text).pack(side=tk.LEFT)

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


class WarningIconWidget(tkinter.ttk.Label, object):
    def __init__(self, parent, text):
        """
        Widget with waring icon and Tooltip information

        :param parent: Parent of widget.
        :param text: Content of tooltip
        """
        self.image = tk.PhotoImage(file=get_img("warning.gif"))
        super(WarningIconWidget, self).__init__(parent, image=self.image, padding=0)

        ToolTip.create(self, text)


class InfoIconWidget(tkinter.ttk.Label, object):
    def __init__(self, parent, text):
        """
        Widget with info icon and Tooltip information

        :param parent: Parent of widget.
        :param text: Content of tooltip
        """
        self.image = tk.PhotoImage(file=get_img("info.gif"))
        super(InfoIconWidget, self).__init__(parent, image=self.image, padding=0)

        ToolTip.create(self, text)


class HidingFrame(tkinter.ttk.Frame, object):
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

        tkinter.ttk.Label(self, text=text, style="HF.TFrame.Label").grid(
            sticky="W", row=0, column=0, columnspan=2, padx=5, pady=5)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

        self.row = row

    def show(self):
        """ Method to grid Frame with predefinied configuration. """
        self.grid(row=self.row, column=0, columnspan=2, padx=10, pady=10, ipady=5, ipadx=10)


class CallbackWrapper(object):
    def __init__(self, callback, *args, **kwargs):
        """ Allow to use callbacks with predefined list of arguments. """
        self.callback = callback
        self.args = args
        self.kwargs = kwargs

    def __call__(self, *args, **kwargs):
        self.callback(*self.args, **self.kwargs)


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
        label = tkinter.ttk.Label(tw, text=self.text, justify=tk.LEFT,
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
        self.canvas = tk.Canvas(self, bd=0, highlightthickness=0,
                                yscrollcommand=vscrollbar.set)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=tk.TRUE)
        vscrollbar.config(command=self.canvas.yview)

        # reset the view
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)

        # create a frame inside the canvas which will be scrolled with it
        self.interior = interior = tk.Frame(self.canvas)
        interior_id = self.canvas.create_window(0, 0, window=interior, anchor="nw")

        # track changes to the canvas and frame width and sync them,
        # also updating the scrollbar
        def _configure_interior(e):
            # update the scrollbars to match the size of the inner frame
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            self.canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != self.canvas.winfo_width():
                # update the canvas's width to fit the inner frame
                self.canvas.config(width=interior.winfo_reqwidth())

        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(e):
            if interior.winfo_reqwidth() != self.canvas.winfo_width():
                # update the inner frame's width to fill the canvas
                self.canvas.itemconfigure(interior_id, width=self.canvas.winfo_width())

        self.canvas.bind('<Configure>', _configure_canvas)

        self.canvas.bind('<Enter>', self._bind_mousewheel)
        self.canvas.bind('<Leave>', self._unbind_mousewheel)

    def _bind_mousewheel(self, e):
        def _mousewheel_handler(e, canvas):
            if e.num == 5:
                canvas.yview_scroll(1, "units")
            if e.num == 4:
                canvas.yview_scroll(-1, "units")

        self.canvas.bind_all("<MouseWheel>", lambda e: _mousewheel_handler(e, self.canvas))
        self.canvas.bind_all("<Button-4>", lambda e: _mousewheel_handler(e, self.canvas))
        self.canvas.bind_all("<Button-5>", lambda e: _mousewheel_handler(e, self.canvas))

    def _unbind_mousewheel(self, e):
        self.canvas.unbind("<MouseWheel>")
