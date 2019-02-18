# -*- coding: utf-8 -*-

import ConfigParser
import os

# TODO: Defaults should be owned by configcheck, not valveconfig
import aquaduct.apps.valveconfig.defaults as defaults


class NonExistingValue(Exception):
    pass


class ConfigCheck(object):
    def __init__(self, config_filename):
        self.config = ConfigParser.RawConfigParser()

        with open(config_filename, "r") as config_file:
            self.config.readfp(config_file)
    def check(self):
        """
        Check if all config options are valid.

        :return: None if options are correct. Otherwise raise proper exception.
        """
        for section_name in self.config.sections():
            for option_name in self.config.options(section_name):
                # it should use DEFAULTS values to determine proper method to get option value
                value = self.config.get(section_name, option_name)
                self.valid(section_name, option_name, value)

    # value is casted to type of target_value
    @staticmethod
    def dynamic_cast(value, target_value):
        try:
            result = type(target_value)(value)
            return result
        except (TypeError, ValueError) as e:
            pass  # types not compatible, effectively returning None

    @staticmethod
    def valid(section_name, option_name, value):
        # TODO: value must be casted to the right type
        # strings to be directly compared, empty string is wildcard

        option_info = defaults.get_default_entry(section_name, option_name)

        if len(option_info.default_values) == 1:
            if isinstance(option_info.default_values[0], tuple):
                if value not in option_info.default_values[0]:

                    raise NonExistingValue("")
            elif isinstance(option_info.default_values[0], list):
                # Combobox always contain string, because it's universal and can keep bool/str/int
                pass
            elif isinstance(option_info.default_values[0], defaults.filetype):
                if not os.path.isfile(value):
                    raise NonExistingValue("File does not exist")
            elif isinstance(option_info.default_values[0], defaults.manyfiletype):
                files = value.split(os.pathsep)
                for elem in files:
                    if not os.path.isfile(elem):
                        raise NonExistingValue("File does not exist")
            else:
                if not isinstance(ConfigCheck.dynamic_cast(value, option_info.default_values[0]), type(option_info.default_values[0])):
                    raise TypeError()

        return True

if __name__ == '__main__':
    confCheck = ConfigCheck("/home/kara/Desktop/aq dev files/config.txt")
    confCheck.check()
