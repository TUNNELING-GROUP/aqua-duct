#!/usr/bin/env python2.7
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

import gzip
import json
import random
import string

from numpy import log
from pymol import cmd, stored


def _sele_exists(sele):
    return sele in cmd.get_names("all")


def _random_string(n=10):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))


def hs_resize(meta_file, selection):
    """
    Changes size of hotspots depending on normalized partial_charge values.

    USAGE:
    hs_resize /path/to/pond_meta.json hotspot_selection
    """
    if not _sele_exists(selection):
        raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

    # Find free sele name
    temp_sele = _random_string()
    while _sele_exists(temp_sele):
        temp_sele = _random_string()

    states = cmd.count_states(selection=selection)

    with gzip.open(meta_file) as f:
        ref = float(json.load(f)["reference_density_correction"])

    for state in range(1, states + 1):
        stored.info = []
        cmd.iterate_state(state, selection, "stored.info.append((ID, partial_charge))")

        for id_, partial_charge in stored.info:
            size = log(partial_charge / ref * 1. + 1)

            cmd.select(temp_sele, "{} and id {}".format(selection, id_), state=state)
            cmd.set("sphere_scale", value=size, selection=temp_sele)
            cmd.alter(temp_sele, "b={}".format(partial_charge))

    cmd.delete(temp_sele)


cmd.extend("hs_resize", hs_resize)
