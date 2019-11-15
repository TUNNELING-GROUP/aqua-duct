# -*- coding: utf-8 -*-

# Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
# Copyright (C) 2018-2019  Michał Banas <info@aquaduct.pl>
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

"""Module used for caching results of specific function using database to store data"""

from hashlib import md5
from aquaduct import logger
import numpy as np


class DBCache(object):
    def __init__(self, db_dir):
        self.db_dir = db_dir

    def __call__(self, func, *args, **kwargs):
        def wrapper(*args, **kwargs):
            key = md5(','.join(map(str, args)) + '&' + ','.join(
                [':'.join(map(str, kv)) for kv in iter(kwargs.items())])).hexdigest()

            logger.debug('Looking for cache key {}'.format(key))

            try:
                coords = np.load(self.db_dir + key + ".npy")
                print("Data fetched")
                return coords
            except IOError:
                coords = func(*args, **kwargs)

                np.save(self.db_dir + key, coords)

                logger.debug("New key added to cache.")

                return coords

        return wrapper
