# -*- coding: utf-8 -*-

import ConfigParser

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
                # TODO: To get value in right type,
                # it should use DEFAULTS values to determine proper method to get option value
                value = self.config.get(section_name, option_name)
                self.valid(section_name, option_name, value)

    @staticmethod
    def valid(section_name, option_name, value):
        # TODO: value must be casted to the right type
        option_info = defaults.get_default_entry(section_name, option_name)

        if len(option_info.default_values) == 1:
            if isinstance(option_info.default_values[0], tuple):
                if value not in option_info.default_values[0]:
                    raise NonExistingValue()
            elif isinstance(option_info.default_values[0], list):
                # Combobox always contain string, because it's universal and can keep bool/str/int
                pass
            else:
                if not isinstance(value, option_info.default_values[0]):
                    raise TypeError()
        # There is 2 options only when there is file/many files loading and
        # checkbox is used(which was replaced by combobox and not used for now).
        elif len(option_info.default_values) == 2:
            if isinstance(option_info.default_values[0], tuple):
                if value not in option_info.default_values[0]:
                    raise NonExistingValue()
            elif isinstance(option_info.default_values[0], list):
                # Combobox again
                pass
            elif isinstance(option_info.default_values[1], bool):
                # When checkbox is used, boolean value can be False only, otherwise it have type of first element in list
                if not isinstance(value, option_info.default_values[0]) or (isinstance(value, bool) and value):
                    raise TypeError()
            else:
                if not isinstance(value, option_info.default_values[0]):
                    raise TypeError()

        return True
