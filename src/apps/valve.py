

'''
This is driver for aqueduct.
'''

import ConfigParser

def get_stage_names(nr=None):
    if nr is None:
        return [get_stage_names(nr) for nr in range(5)]
    else:
        if nr == 0:
            return 'find traceable residues'
            pass
        elif nr == 1:
            return 'find raw paths'
        elif nr == 2:
            return 'create separate frames'
        elif nr == 3:
            pass
        elif nr == 4:
            pass


def get_default_config():

    config = ConfigParser.RawConfigParser()

    def common(section):
        config.set(section, 'exec', 'skip') # what to do: skip, run
        config.set(section, 'load') # load previous results form file name
        config.set(section, 'save') # save results to file name

    def common_traj_data(section):
        config.set(section, 'top') # top file name
        config.set(section, 'nc') # netcdf file name
        config.set(section, 'scope') # scope definition
        config.set(section, 'scope_convexhull','true') # take convex hull of scope
        config.set(section, 'object') # object definition

    ################
    # stage I
    section = 'find traceable residues'
    config.add_section(section)

    common(section)
    common_traj_data(section)

    ################
    # stage II
    section = 'find raw paths'
    config.add_section(section)

    common(section)
    common_traj_data(section)

    config.set(section,'discard_empty_paths','true')

    ################
    # stage III
    section = 'create separate frames'

    common(section)

    config.set(section,'discard_empty_paths','true')
