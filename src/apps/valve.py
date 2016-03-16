'''
This is driver for aqueduct.
'''

import ConfigParser


class ValveConfig(object):
    def __init__(self):
        self.config = self.get_default_config()

    def common_config_names(self):
        # exec - what to do: skip, run
        # load - load previous results form file name
        # save - save results to file name
        return 'exec load save'.split()

    def common_traj_data_config_names(self):
        # top - top file name
        # nc - netcdf file name
        # scope - scope definition
        # scope_convexhull - take convex hull of scope, true of false
        # object - object definition
        return 'top nc scope scope_convexhull object'.split()

    def stage_names(self, nr=None):
        if nr is None:
            return [self.stage_names(nr) for nr in range(5)]
        else:
            if nr == 0:
                return 'find traceable residues'
            elif nr == 1:
                return 'find raw paths'
            elif nr == 2:
                return 'create separate frames'
            elif nr == 3:
                return 'inlets clusterisation'
            elif nr == 4:
                return 'analysis'
        raise NotImplementedError('Stage %r is not implemented.' % nr)

    def get_common_traj_data(self, stage):
        assert isinstance(stage, int)
        options = {name: None for name in self.common_traj_data_config_names()}
        for nr in range(stage + 1)[::-1]:
            section = self.stage_names(nr)
            for name in self.common_traj_data_config_names():
                if self.config.has_option(section, name):
                    value = self.config.get(section, name)
                    if (value is not None) and (options[name] is None):
                        options.update({name: value})
        return options

    def get_options(self, stage):
        assert isinstance(stage, int)
        stage_name = self.stage_names(stage)
        names = self.config.options(stage_name)
        options = {name: self.config.get(stage_name, name) for name in names}
        options.update(self.get_common_traj_data(stage))
        return options

    def get_default_config(self):
        config = ConfigParser.RawConfigParser()

        def common(section):
            for setting in self.common_config_names():
                config.set(section, setting)

        def common_traj_data(section):
            for setting in self.common_traj_data_config_names():
                config.set(section, setting)

        ################
        # stage I
        # find traceable residues
        section = self.stage_names(0)
        config.add_section(section)

        common(section)
        common_traj_data(section)

        ################
        # stage II
        # find raw paths
        section = self.stage_names(1)
        config.add_section(section)

        common(section)
        common_traj_data(section)


        ################
        # stage III
        # create separate frames
        section = self.stage_names(2)
        config.add_section(section)

        common(section)

        config.set(section, 'discard_empty_paths', 'true')

        ################
        # stage IV
        # inlets clusterisation
        section = self.stage_names(3)
        config.add_section(section)

        common(section)

        ################
        # stage V
        # analysis
        section = self.stage_names(4)
        config.add_section(section)

        common(section)

        return config

    def load_config(self, filename):
        self.config.read(filename)

    def save_config(self, filename):
        with open(filename, 'w') as f:
            self.config.write(f)

