# -*- coding: utf-8 -*-

import logging
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

from aquaduct import logger
from aquaduct.apps.valvecore import load_stage_dump

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)  # default level is WARNING
logger.addHandler(ch)

import argparse

# load_dump functions&firends will be replaced by imports


if __name__ == "__main__":
        description = '''Mario, a tool for trajectories' properties visualization'''# TODO: przemyslec
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-i", action="store", dest="filename", required=True, help="Input file with paths to visualize.")
        parser.add_argument("-n", action="store", dest="num", required=True, help="ID of molecule to visualize.")
        parser.add_argument("-p", action="store", dest="plot_type", required=False, help="Type of plot(s) to be displayed. Separated by comma.")
        #Dodac nieobowiazkowe parametry polyorder i window_length
        #Kolory (midpoints)
        #create help!
        #0. Create list^ of plots available, 2. Read the string with names from -p args, 2. Divide by comma (list* will be created) 3. Check elements of the list* for any misamtched with list^. 4. Think about automatization of plotting process
        args = parser.parse_args()
        # Argument checking section
        try:
            args.num = int(args.num)
        except ValueError:
            logger.error("Oops!  That was no valid number. Integer required.")
            raise ValueError("Oops!  That was no valid number. Integer required.")

        loaded_data = load_stage_dump(args.filename)
        assert ((args.num >= 0) & (args.num <= len(loaded_data['spaths']))), "Path number should be in range 0 - " + str(len(loaded_data['spaths'])) #Check if path ID is accessible

        divided_paths = map(len, loaded_data['spaths'][args.num].paths)
        velocity = loaded_data['spaths'][args.num].get_velocity_cont() #Getting velocity for each frame the path with ID = args.num was detected
        SG_filtered = savgol_filter(x=velocity, window_length=75, polyorder=2)
        plt.plot()
        plt.plot(velocity)
        plt.show()

        distance = loaded_data['spaths'][args.num].get_distance_cont()  #Getting distance travelled for each frame the path with ID = args.num was detected
        plt.plot(distance)


        plt.plot(distance,velocity)
        plt.show()
