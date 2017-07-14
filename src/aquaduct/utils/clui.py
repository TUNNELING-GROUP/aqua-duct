# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>
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

"""
Module comprises convieniences functions and definitios for different operations related to command line user interface.
"""

import logging
logger = logging.getLogger(__name__)

from collections import OrderedDict
from aquaduct import logger as root_logger
from aquaduct import __mail__ as mail
import datetime
import time
from os import linesep
from sys import stderr
from os import linesep
from functools import partial

from multiprocessing import Queue, Manager, Lock, Value, Process
from collections import OrderedDict

def emit_message_to_file_in_root_logger(mess):
    # emits message to the file used by file handler in the root logger
    # assumes there is only one file handler
    if logging.FileHandler in map(type, root_logger.handlers):
        fh = root_logger.handlers[map(type, root_logger.handlers).index(logging.FileHandler)]
        with fh.lock:
            with open(fh.baseFilename, 'a') as logfile:
                logfile.write(mess)


def message_special(mess):
    emit_message_to_file_in_root_logger(mess)


def message(mess, cont=False):
    """
    Prints message to standard error.
    If FileHandler is present in the :class:`root_logger` the same message is appended to the log file.

    :param str mess: message to print
    :param bool cont: if set True no new line is printed
    """
    if cont:
        mess = mess + ' '
    else:
        mess = mess + linesep
    emit_message_to_file_in_root_logger(mess)
    stderr.write(mess)


class fbm(object):
    # feedback message
    def __init__(self, info, cont=True):
        self.__cont = cont
        self.__info = info
        message(info + '...', cont=self.__cont)

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        if typ is None:
            if self.__cont:
                message("OK.")
            else:
                message(self.__info + ": DONE.")

    def __call__(self, info):
        if self.__cont:
            message(linesep + '\t' + info + '...', cont=True)
        else:
            message(info + '...', cont=False)


class tictoc(object):
    def __init__(self, mess):
        self.__tic = 0
        self.__toc = 0
        self.__mess = mess

    def __enter__(self):
        self.__tic = time.time()
        return self

    def __exit__(self, typ, value, traceback):
        self.__toc = time.time()
        if typ is None:
            #logger.debug('Execution time of [%s] %s', self.__mess, smart_time_string(self.__toc - self.__tic))
            logger.debug('Execution time of [%s] %f', self.__mess, (self.__toc - self.__tic))


gregorian_year_in_days = 365.2425
'''Length of Gregorian year in days. Average value. Source: https://en.wikipedia.org/wiki/Year'''


def smart_time_string(s, rl=0, t=1.1, maximal_length=None, maximal_units=5):
    """
    Function transforms time in seconds to nicely formatted string of
    length defined by :attr:`maximal_length`. Depending on number of seconds time is represented with
    one or more of the following units:

    ========= =================
    Unit name Unit abbreviation
    ========= =================
    seconds   s
    minutes   m
    hours     h
    days      d
    years     y
    ========= =================

    Maximal number of units used in time string can be set with :attr:`maximal_units`.

    :param int s: Input time in seconds.
    :param int rl: Number of units already used for representing time.
    :param float t: Exces above standard number of current time units.
    :param int maximal_length: Maximal length of the output string. Must be greater then 0.
    :param int maximal_units: Maximal number of units used in the output string. Must be greater then 0 and lower then 6.

    :return: string of nicely formated time
    :rtype: str
    """
    # assert isinstance(maximal_length, (int, long))
    # assert maximal_length > 0
    assert isinstance(maximal_units, (int, long))
    assert maximal_units > 0
    assert maximal_units < 6

    output = ''
    rl += 1
    if rl > maximal_units:
        output = ''
    else:
        # seconds
        if s < t * 60:
            output = "%2.2d s" % s
        elif s < t * 3600:
            output = ("%2.2d m" % (int(s) / 60)) + ' ' + smart_time_string(int(s) % 60, rl).strip()
        elif s < t * 3600 * 24:
            output = ("%2.2d h" % (int(s) / 3600)) + ' ' + smart_time_string(int(s) % 3600, rl).strip()
        elif s < t * 3600 * 24 * gregorian_year_in_days:
            output = ("%d d" % (int(s) / (3600 * 24))) + ' ' + smart_time_string(int(s) % (3600 * 24), rl).strip()
        elif True:
            output = ("%d y" % (int(s) / (3600 * 24 * gregorian_year_in_days))) + ' ' + smart_time_string(
                int(s) % (3600 * 24 * gregorian_year_in_days), rl).strip()
    if maximal_length:
        return (output + " " * maximal_length)[:maximal_length]
    return output


###################################
# vis separators

def gsep(sep='-', times=72, length=None):
    """
    Generic separator.

    :param str sep: Element(s) of separator.
    :param int times: Number of times :attr:`sep` is printed.
    :param int length: Optional maximal length of output.
    :return: String separator.
    :rtype: str
    """
    return (sep * times)[:length]


def tsep(line):
    """
    :param str line: Input line.
    :return: Returns default :func:`gsep` of length of :attr:`line`.
    """
    return gsep(sep='-', times=len(line))


def underline(line):
    """
    :param str line: Input line.
    :return: String made by concatenation of :attr:`line`, :mod:`os.linesep`, and output of :func:`tsep` called with :attr:`line`.
    :rtype: str
    """
    uline = line
    uline += linesep
    uline += tsep(line)
    return uline


def thead(line):
    """
    :param str line: Input line.
    :return: String made by concatenation of output of :func:`tsep` called with :attr:`line`, :attr:`line`, :mod:`os.linesep`, and again output of :func:`tsep` called with :attr:`line`.
    :rtype: str
    """
    header = tsep(line)
    header += linesep
    header += line
    header += linesep
    header += tsep(line)
    return header


###################################


class SimpleProgressBar(object):
    """
    Simple progress bar displaying progress with percent indicator, progress bar and ETA.
    Progress is measured by iterations.

    :cvar str rotate: String comprising characters with frames of a rotating toy.
    :cvar int barlenght: Length of progress bar.
    :ivar int maxval: maximal number of iterations
    :ivar int current: current number of iterations
    :ivar bool overrun_notice: if True, overrun above :attr:`maxval` iterations causes insert of newline
    :ivar bool overrun: flag of overrun
    :ivar int begin: time in seconds at the initialization of the :class:`SimpleProgressBar` class.
    :ivar int tcurrent: time in seconds of current iteration
    """

    rotate = '\\|/-'
    # rotate = '<^>v'
    # rotate = '.:|:.'
    # rotate = 'x+'

    barlenght = 24

    def __init__(self, maxval=None, mess=None):
        """
        :param int maxval: Maximal number of iterations stored to :attr:`maxval`.
        :param str mess: Optional message displayed at progress bar initialization.
        """

        #self.lock = Manager().Lock()

        assert isinstance(maxval,
                          (int, long)), 'Parameter maxval should be of int or long type, %r given instead.' % type(
            maxval)
        if maxval < 1:
            self.maxval = 1
        else:
            self.maxval = maxval

        self.tens = []

        self.current = 0

        self.overrun_notice = True
        self.overrun = False

        self.begin = time.time()
        self.tcurrent = self.begin

        self.last_rotate_time = self.begin
        self.last_rotate_idx = 0

        if mess is not None:
            message(mess)
        self.show()

    def bar(self):
        barval = int(self.percent() / 100 * self.barlenght)
        if barval > self.barlenght:
            barval = self.barlenght
        bar = '#' * barval
        if self.current:
            if self.tcurrent - self.last_rotate_time > 1. / 4:  # FIXME: magic constant, remove it!
                self.last_rotate_idx += 1
                self.last_rotate_time = self.tcurrent
            if self.last_rotate_idx > len(self.rotate) - 1:
                self.last_rotate_idx = 0
            bar += self.rotate[self.last_rotate_idx]
        bar += ' ' * self.barlenght
        return '[%s]' % bar[:self.barlenght]

    def ETA(self):
        """
        Returns ETA calculated on the basis of current number of iterations
        :attr:`current` and current time :attr:`tcurrent`. If number of
        iterations is 0 returns ``?``.
        Time is formated wiht :func:`smart_time_string`.

        :return: ETA as string.
        :rtype: str
        """
        if self.current == 0:
            return '?'
        diff = self.tcurrent - self.begin
        periteration = diff / self.current
        expected = periteration * self.maxval
        eta = periteration * (self.maxval - self.current)

        return smart_time_string(eta)

    def percent(self):
        """
        Returns float number of precent progress calculated in the basis
        of current number of iterations :attr:`current`. Should return
        number between 0 and 100.

        :returns: percent progress number
        :rtype: float
        """
        percent = float(self.current) / float(self.maxval) * 100
        return percent

    def show(self):
        """
        Shows current progress.

        If value returned by :meth:`percent` is =< 100 then progres is printed as
        percent indicator leaded by ETA calculated by :meth:`ETA`.

        If value returned by :meth:`percent` is > 100 then progress is printed as
        number of iterations and total time.

        Progress bar is writen to standard error.
        """
        percent = self.percent()
        # TODO: create some unittests for pbar
        mess = ''
        mess_spec = ''
        if percent > 100 or self.overrun:
            if self.overrun_notice:
                stderr.write(linesep)
                self.overrun_notice = False
                self.overrun = True
            mess_spec = "%d iterations out of %d. Total time: %s" % (self.current, self.maxval, self.ttime())
            mess = "\r" + mess_spec
        elif not self.overrun:
            mess_spec = "%3d%% %s ETA: %s" % (self.percent(), self.bar(), self.ETA())
            mess = "\r" + mess_spec + "\033[K"  # FIXME: magic constant!

        stderr.write(mess)
        # TODO: do not use last_rotate_time here, use separate marker, last_rotate_time can be used in over run notice
        if int(percent) / 10 not in self.tens:  # FIXME: magic constant!
            if percent > 100:
                if self.tcurrent - self.last_rotate_time > 60.:  # FIXME: magic constant, remove it!
                    message_special(mess_spec + linesep)
                    self.tens.append(int(percent) / 10)
                    self.last_rotate_time = self.tcurrent
            else:
                message_special(mess_spec + linesep)
                self.tens.append(int(percent) / 10)

    def heartbeat(self):
        #with self.lock:
        if time.time() - self.last_rotate_time > 2.:  # FIXME: magic constant, remove it!
            self.tcurrent = time.time()
            self.show()

    def next(self):
        return self.update(self.current+1)

    def update(self, step):
        """
        Updates number of current iterations :obj:`current` by one if :obj:`step` is > 0.
        Otherwise number of current iterations is not updated.
        In boths cases time of current iteration :obj:`tcurrent` is updated and
        :meth:`show` is called.

        :param int step: update step
        """
        # TODO: change logic of step == 1 vs step > 1 - add or set?
        #with self.lock:
        if step > 0:
            self.current = step
        self.tcurrent = time.time()
        if (step == self.maxval) or (
                        self.tcurrent - self.last_rotate_time > 1. / 4):  # FIXME: magic constant, remove it!
            # TODO: check for last_rotate_time is done twice, SimpleProgressBar code needs revision
            self.show()

    def ttime(self):
        """
        Calculates and returns total time string formated with :func:`smart_time_string`.

        :return: string of total time
        :rtype: str
        """
        return smart_time_string(self.tcurrent - self.begin)

    def finish(self):
        """
        Finishes progress bar. First, :meth:`update` is called with :obj:`step` = 0. Next message of total time
        is writen to standard error.
        """
        if self.current < self.maxval:
            self.update(self.maxval)
        else:
            self.update(0)
            self.show()
        stderr.write(linesep)
        message("Total time: %s" % self.ttime())
        # stderr.write(linesep)


pbar = SimpleProgressBar  # default progress bar


def get_str_timestamp():
    # returns time stamp as string
    return str(datetime.datetime(*tuple(time.localtime())[:6]))



class SimpleTree(object):
    def __init__(self,name=None,message=None):
        self.name = name
        self.message = []
        self.add_message(message)
        self.branches = []

    def __repr__(self):
        return "%s {%s} %s" % (str(self.name), "; ".join(self.message),str(self.branches))

    def is_leaf(self):
        return len(self.branches)==0

    @property
    def leafs_names(self):
        return [leaf.name for leaf in self.branches]

    def get_leaf(self,name):
        assert name in self.leafs_names
        return [leaf for leaf in self.branches if name == leaf.name][0]


    def add_message(self,message=None,toleaf=None,replace=False):
        if toleaf is not None:
            return self.add_message_to_leaf(message=message,toleaf=toleaf,replace=replace)
        if message is not None:
            if isinstance(message,list):
                if replace:
                    self.message = message
                else:
                    self.message += message
            else:
                if replace:
                    self.message = [message]
                else:
                    self.message += [message]

    def add_message_to_leaf(self,message=None,toleaf=None,replace=False):
        if toleaf in self.leafs_names:
            leaf = self.get_leaf(toleaf)
            return leaf.add_message(message,replace=replace)
        else:
            for leaf in self.branches:
                leaf.add_message_to_leaf(message=message,toleaf=toleaf,replace=replace)

    def add_leaf(self,name=None,message=None,toleaf=None):
        if toleaf is not None:
            return self.add_leaf_to_leaf(name=name,message=message,toleaf=toleaf)
        leaf = SimpleTree(name=name,message=message)
        self.branches.append(leaf)

    def add_leaf_to_leaf(self,name=None,message=None,toleaf=None):
        if toleaf in self.leafs_names:
            leaf = self.get_leaf(toleaf)
            return leaf.add_leaf(name=name,message=message)
        else:
            for leaf in self.branches:
                leaf.add_leaf_to_leaf(name=name,message=message,toleaf=toleaf)


def print_simple_tree(st,prefix=None,multiple=False,concise=True):
    _l = '|'
    _t = ' '
    _c = '-+'+_t

    def name_str(name):
        if name is None:
            return ''
        return str(name)
    def name_len(name):
        return len(name_str(name))

    def message_str(message):
        if len(message):
            return '{%s}' % '; '.join(message)
        return ''

    prefix_ = ''
    if prefix is not None:
        prefix_ += prefix
    out = ''
    # name
    out += prefix_
    out += name_str(st.name)
    if not st.is_leaf():
        if name_len(st.name):
            out += _c
        out += message_str(st.message)
        out += linesep
        if multiple:
            new_prefix = prefix_+_l+(_t*(name_len(st.name)-1))
        else:
            new_prefix = prefix_+(_t*(name_len(st.name)))
        new_prefix += _t
        if not concise:
            out += new_prefix+_l+linesep
        out_rec = []
        for nr,branch in enumerate(st.branches):
            out_rec.append(partial(print_simple_tree,prefix=new_prefix, multiple=len(st.branches) - nr > 1)(branch))
        out += ''.join(out_rec)
        if st.branches[-1].is_leaf() and not concise:
            out += new_prefix.rstrip(_t) + linesep
    else:
        out += _t
        out += message_str(st.message)
        out += linesep
    return out
