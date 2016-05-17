"""
Module comprises convieniences functions and definitios for logging
purposes including progress bar helpers.
"""

# __all__ = ['debug','info','warning','error','critical','pbar']

import logging
from os import linesep
from sys import stderr
import time
import datetime

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
    :ivar bool overrun_notice: if True, overrun above :ivar:`maxval` iterations causes insert of newline
    :ivar bool overrun: flag of overrun
    :ivar int begin: time in seconds at the initialization of the :class:`SimpleProgressBar` class.
    :ivar int tcurrent: time in seconds of current iteration
    """

    rotate = '\\|/-'
    # rotate = '<^>v'
    # rotate = '.:|:.'
    # rotate = 'x+'

    barlenght = 24

    def __init__(self, mess=None, maxval=None):
        """
        :param int maxval: Maximal number of iterations stored to :attr:`maxval`.
        :param str mess: Optional message displayed at progress bar initialization.
        """

        assert isinstance(maxval, (int, long))
        if maxval < 1:
            self.maxval = 1
        else:
            self.maxval = maxval

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
        if percent > 100:
            if self.overrun_notice:
                stderr.write(linesep)
                self.overrun_notice = False
                self.overrun = True
            stderr.write("\r%d iterations out of %d. Total time: %s" % (self.current, self.maxval, self.ttime()))
        elif not self.overrun:
            stderr.write(
                "\r%3d%% %s ETA: %s" % (self.percent(), self.bar(), self.ETA()) + "\033[K")  # FIXME: magic constant!

    def update(self, step):
        """
        Updates number of current iterations :obj:`current` by one if :obj:`step` is > 0.
        Otherwise number of current iterations is not updated.
        In boths cases time of current iteration :obj:`tcurrent` is updated and
        :meth:`show` is called.

        :param int step: update step
        """
        if step > 0:
            if step == 1:
                self.current += 1
            else:
                self.current = step
        self.tcurrent = time.time()
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
        stderr.write(linesep)
        stderr.write("Total time: %s" % self.ttime())
        stderr.write(linesep)


# TODO: try to add tqdm as an alternative for progressbarr which does not work very well in ipython


class pbar(object):
    """
    Progress bar wrapper class.
    It can use several types of progress bars, including :class:`SimpleProgressBar`. Additionaly, it can
    handle progress bars with following packages (must be installed separately):

    * progressbar
    * tqdm
    * pyprind

    :ivar int __maxval: maximal number of iterations
    :ivar str __kind: type of progress bar
    :ivar int __curval: current number of iterations
    :ivar __pbar: progress bar child object

    .. warning::

        :class:`SimpleProgressBar` is the only one kind of progress bar recommended.

    """

    def __init__(self, maxval=100, kind='simple'):
        """
        :param int maxval: maximal number of iterations stored to :ivar:`__maxval` and passed child progress bar object
        :param str kind: type of progress bar, available types: simple, progressbar, tqdm, pyprind
        """

        assert kind in ['simple', 'progressbar', 'tqdm', 'pyprind']

        self.__kind = kind
        self.__maxval = maxval
        self.__curval = 0

        if self.__kind == "simple":
            self.__pbar = SimpleProgressBar(maxval=maxval)
        if self.__kind == "progressbar":
            import progressbar as pb
            widgets = [pb.Percentage(), ' ',
                       pb.Bar(), ' ',
                       pb.ETA()]
            self.__pbar = pb.ProgressBar(widgets=widgets, maxval=maxval).start()
        elif self.__kind == "tqdm":
            # it does not work as expected
            from tqdm import tqdm
            self.__pbar = tqdm(leave=True, ascii=True, total=maxval)
        elif self.__kind == "pyprind":
            import pyprind
            self.__pbar = pyprind.ProgBar(maxval)

    def update(self, val):
        """
        Updates progress bar with value of :obj:`val` parameter. Exact behavior depends
        on the type of progress bar.

        :param int val: value used to update progress bar
        """
        if self.__kind in ["tqdm", "pyprind"]:
            if val < 1:
                return
        if self.__kind in ["tqdm", "pyprind"]:
            self.__pbar.update()
        else:
            self.__pbar.update(val)

    def finish(self):
        """
        Finishes progress bar. It cals appropriate, if exist, method of child progress bar object.
        is writen to standard error.
        """
        if self.__kind in ["progressbar", "simple"]:
            self.__pbar.finish()
        if self.__kind == "tqdm":
            self.__pbar.close()


def get_str_timestamp():
    return str(datetime.datetime(*tuple(time.localtime())[:6]))


level = logging.WARNING
# format = linesep+'AQUARIUM:%(levelname)1.1s:[%(module)s|%(funcName)s@s%(lineno)d]:'+linesep+'%(message)s'
log_format = linesep + 'AQUARIUM:%(levelname)s:[%(module)s|%(funcName)s@s%(lineno)d]:' + linesep + '%(message)s'

logging.basicConfig(format=log_format, level=level)

debug = logging.debug
info = logging.info
warning = logging.warning
error = logging.error
critical = logging.critical


class fbm(object):
    # feedback message
    def __init__(self, info):
        message(info + '...', cont=True)

    def __enter__(self):
        return self

    def __exit__(self, typ, value, traceback):
        if typ is None:
            message("OK.")
        else:
            raise typ(value)

    def __call__(self, info):
        message(linesep + '\t' + info + '...', cont=True)


def message(mess, cont=False):
    """
    Prints message to standard error.

    :param str mess: message to print
    :param bool cont: if set True no new line is printed
    """
    if cont:
        print >> stderr, mess,
    else:
        print >> stderr, mess
