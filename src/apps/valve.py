#!/usr/bin/env python2.7
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
This is a driver for Aqua-Duct.
"""

################################################################################
# reuse AQ logger

import logging
from aquaduct import logger, logger_name

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING) # default level is WARNING
logger.addHandler(ch)

################################################################################


from aquaduct.apps.valvecore import *


################################################################################


if __name__ == "__main__":

    with clui.tictoc('Aqua-Duct calculations'):

        ############################################################################
        # argument parsing
        import argparse

        description_version = '''Aquaduct library version %s''' % (aquaduct_version_nice(), )
        description = '''Valve, Aquaduct driver'''

        parser = argparse.ArgumentParser(description=description,
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument("--debug", action="store_true", dest="debug", required=False, help="Prints debug info.")
        parser.add_argument("--debug-file", action="store", dest="debug_file", required=False, help="Debug log file.")
        parser.add_argument("--dump-template-config", action="store_true", dest="dump_template_conf", required=False,
                            help="Dumps template config file. Suppress all other output or actions.")
        parser.add_argument("-t", action="store", dest="threads", required=False, default=None,
                            help="Limit Aqua-Duct calculations to given number of threads.")
        parser.add_argument("-c", action="store", dest="config_file", required=False, help="Config file filename.")
        parser.add_argument("--sps", action="store_true", dest="sps", required=False, help="Use single precision to store data.")
        parser.add_argument("--max-frame", action="store", dest="max_frame", type=int, required=False, help="Maximal number of frame.")
        parser.add_argument("--min-frame", action="store", dest="min_frame", type=int, required=False, help="Minimal number of frame.")
        parser.add_argument("--step-frame", action="store", dest="step_frame", type=int, required=False, help="Frames step.")
        parser.add_argument("--version", action="store_true", dest="print_version", required=False,
                            help="Prints versions and exits.")
        parser.add_argument("--license", action="store_true", dest="print_license", required=False,
                            help="Prints short license info and exits.")

        args = parser.parse_args()

        ############################################################################
        # debug
        # at this stage logger is the AQ root logger
        if args.debug:
            logger.removeHandler(ch)  # remove old ch handlers
            ch = logging.StreamHandler()
            ch.setFormatter(formatter)
            ch.setLevel(logging.DEBUG)
            logger.addHandler(ch)
        if args.debug_file:
            formatter = logging.Formatter('%(asctime)s: ' + formatter_string)
            fh = logging.FileHandler(args.debug_file)
            fh.setFormatter(formatter)
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)
        # finally, get valve logger
        logger = logging.getLogger(logger_name + '.valve')
        logger.info('Initialization of Valve logging done.')

        ############################################################################
        # single precision storage
        if args.sps:
            logger.info('Single precision data storage activated.')
            from aquaduct.utils.maths import defaults
            from numpy import float32, int32
            defaults.float_default = float32
            defaults.int_default = int32

        ############################################################################
        # special option for dumping template config
        config = ValveConfig()  # config template
        if args.dump_template_conf:
            #import cStringIO as StringIO
            #config_dump = StringIO.StringIO()
            #config.save_config_stream(config_dump)
            #print config_dump.getvalue()
            import os
            print os.linesep.join(config.dump_config(dump_template=True))
            exit(0)
        # special case of version
        if args.print_version:
            print description
            print description_version
            exit(0)
        # special case of license
        if args.print_license:
            valve_begin()
            print "Licensed under GNU GPL v3. Full text of the license is distributed"
            print "with installation package and is also available at"
            print "https://www.gnu.org/licenses/gpl-3.0.txt"
            print ""
            print '''Aqua-Duct, a tool facilitating analysis of the flow of solvent molecules in molecular dynamic simulations
    Copyright (C) 2016-2017  Tomasz Magdziarz, Alicja Płuciennik, Michał Stolarczyk <info@aquaduct.pl>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    '''
            exit(0)


        ############################################################################
        # begin!

        valve_begin()
        valve_load_config(args.config_file, config)

        # get global options
        goptions = config.get_global_options()
        # pbar_name = goptions.pbar

        if args.threads is None:
            optimal_threads.threads_count = optimal_threads.cpu_count + 1
        else:
            optimal_threads.threads_count = int(args.threads)
        clui.message("Number of threads Valve is allowed to use: %d" % optimal_threads.threads_count)
        if (1 < optimal_threads.threads_count < 3) or (optimal_threads.threads_count - 1 > optimal_threads.cpu_count):
            clui.message("Number of threads is not optimal; CPU count reported by system: %d" % optimal_threads.cpu_count)
        # because it is used by mp.Pool it should be -1???
        if optimal_threads.threads_count > 1:
            optimal_threads.threads_count -= 1
            clui.message("Main process would use 1 thread.")
            clui.message("Concurent calculations would use %d threads." % optimal_threads.threads_count)

        # At this point calculations starts. All options are read.
        # Options:
        # goptions - global options
        # config - configuration file object

        ############################################################################
        # STAGE 0

        # Maximal frame checks
        frames_window = slice(args.min_frame,args.max_frame,args.step_frame)

        # TODO: Is it always required?
        reader = valve_read_trajectory(goptions.top, goptions.trj, frames_window=frames_window) # trajectory reader

        with reader.get() as traj_reader:
            clui.message("Frames window: %d:%d step %d" % (traj_reader.get_start_frame(),
                                                           traj_reader.get_stop_frame(),
                                                           traj_reader.get_step_frame()))

        ## TODO: Is it reported correctly?
        #clui.message("Using %d of %d available frames." % (max_frame + 1, reader.max_frame + 1))

        # container for collecting whether particular stage was executed
        run_status = {}

        # STAGE I
        result1 = valve_exec_stage(0, config, stage_I_run,
                                   run_status=run_status,
                                   reader=reader)

        # STAGE II
        result2 = valve_exec_stage(1, config, stage_II_run,
                                   run_status=run_status,
                                   reader=reader,
                                   **result1)

        # STAGE III
        result3 = valve_exec_stage(2, config, stage_III_run,
                                   run_status=run_status,
                                   reader=reader,
                                   **result2)

        # STAGE IV
        result4 = valve_exec_stage(3, config, stage_IV_run,
                                   run_status=run_status,
                                   reader=reader,
                                   **result3)

        # STAGE V
        results = {}
        for result in (result2, result3, result4):
            results.update(result)

        result5 = valve_exec_stage(4, config, stage_V_run,
                                   run_status=run_status,
                                   no_io=True,
                                   **results)

        # STAGE VI
        results = {}
        for result in (result3, result4):
            results.update(result)

        result6 = valve_exec_stage(5, config, stage_VI_run,
                                   run_status=run_status,
                                   no_io=True,
                                   reader=reader,
                                   **results)
        ############################################################################
        # end!

        valve_end()
        logger.info('Valve calulations finished.')
