import logging
from aquaduct import logger, logger_name

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)  # default level is WARNING
logger.addHandler(ch)

import cPickle as pickle
import gzip
from aquaduct.utils import clui
from aquaduct import version as aquaduct_version
from collections import namedtuple, OrderedDict  # TODO: check if OrderedDict is REALLY used


# load_dump functions&firends will be replaced by imports

def version():
    return 7, 7, 7


class LoadDumpWrapper(object):
    """This is wrapper for pickled data that provides compatibility
    with earlier versions of Aqua-Duct.

    Conversions in use:

    1) replace 'aquaduct.' by 'aquaduct.'

    """

    def __init__(self, filehandle):
        self.fh = filehandle

    def convert(self, s):
        new_s = s
        new_s = new_s.replace('aqueduct.', 'aquaduct.')
        new_s = new_s.replace('aqueduct_version', 'aquaduct_version')
        return new_s

    def read(self, *args, **kwargs):
        return self.convert(self.fh.read(*args, **kwargs))

    def readline(self, *args, **kwargs):
        return self.convert(self.fh.readline(*args, **kwargs))


def load_dump(filename):
    with clui.fbm('Loading data dump from %s file' % filename):
        with gzip.open(filename, mode='r') as protof:
            f = LoadDumpWrapper(protof)
            # version!
            loaded_data = pickle.load(f)
            check_versions(loaded_data)

            # loaded data!
            loaded_data = pickle.load(f)
            # enything else?
            try:
                other_data = pickle.load(f)
                loaded_data.update({'other_data': other_data})
            except:
                pass
        return loaded_data
        # loaded_data_nt = namedtuple('LoadedData', loaded_data.keys())
        # return loaded_data_nt(**loaded_data)


def check_version_compliance(current, loaded, what):
    if current[0] > loaded[0]:
        logger.error('Loaded data has %s major version lower then the application.' % what)
    if current[0] < loaded[0]:
        logger.error('Loaded data has %s major version higher then the application.' % what)
    if current[0] != loaded[0]:
        logger.error('Possible problems with API compliance.')
    if current[1] > loaded[1]:
        logger.warning('Loaded data has %s minor version lower then the application.' % what)
    if current[1] < loaded[1]:
        logger.warning('Loaded data has %s minor version higher then the application.' % what)
    if current[1] != loaded[1]:
        logger.warning('Possible problems with API compliance.')


def check_versions(version_dict):
    assert isinstance(version_dict, (dict, OrderedDict)), "File is corrupted, cannot read version data."
    assert 'version' in version_dict, "File is corrupted, cannot read version data."
    assert 'aquaduct_version' in version_dict, "File is corrupted, cannot read version data."
    check_version_compliance(aquaduct_version(), version_dict['aquaduct_version'], 'Aqua-Duct')
    check_version_compliance(version(), version_dict['version'], 'Valve')
#testowy komentarz
