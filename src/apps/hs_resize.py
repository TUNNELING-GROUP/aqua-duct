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

from pymol import cmd, stored
import random
import string


def _sele_exists(sele):
	return sele in cmd.get_names("all")


def _random_string(n = 10):
	return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))


def hs_resize(selection):
	"""
Changes size of hotspots depending on normalized partial_charge values.

USAGE:
hs_resize selection
	"""
	if not _sele_exists(selection):
		raise RuntimeError("Selection \"{}\" does not exists.".format(selection))

	# Find free sele name
	temp_sele = _random_string()
	while _sele_exists(temp_sele):
		temp_sele = _random_string()

	states = cmd.count_states(selection=selection)

	# Find min/max
	min_ = 0
	max_ = 0
	for state in range(1, states+1):
		stored.info = []
		cmd.iterate_state(state, "all", "stored.info.append((partial_charge))")

		state_min = min(stored.info)
		state_max = max(stored.info)

		if state_min < min_:
			min_ = state_min

		if state_max > max_:
			max_ = state_max

	for state in range(1, states+1):
		stored.info = []
		cmd.iterate_state(state, "all", "stored.info.append((partial_charge))")

		min_ = min(partial_charge)
		max_ = max(partial_charge)

		for id_, partial_charge in stored.info:
			size = (partial_charge - min_) * (1. - 0.1) / (max_ - min_) + 0.1

			cmd.select(temp_sele, "id {}".format(id_), state=state)
			cmd.set("sphere_scale", value=size, selection=temp_sele)
			cmd.alter(temp_sele, "b={}".format(partial_charge))

	cmd.delete(temp_sele)


cmd.extend("hs_resize", hs_resize)
