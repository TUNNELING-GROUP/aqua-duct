# -*- coding: utf-8 -*-

import logging
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import argparse
from aquaduct import logger
from aquaduct.apps.valvecore import load_stage_dump

formatter_string = '%(name)s:%(levelname)s:[%(module)s|%(funcName)s@%(lineno)d]: %(message)s'
# create and add console handler with WARNING level to the AQ logger
formatter = logging.Formatter(formatter_string)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)  # default level is WARNING
logger.addHandler(ch)

if __name__ == "__main__":
        description = '''Mario, a tool for trajectories' properties visualization'''# TODO: przemyslec
        parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-i", action="store", dest="filename", required=True, help="Input file with paths to visualize.")
        parser.add_argument("-n", action="store", dest="path_number", required=True, help="ID of path to visualize.")
        parser.add_argument("-p", action="store", dest="plot_type", required=True, help="Type of plot(s) to be displayed. Separated by comma. Plots available: velocity, acceleration", default=None)
        parser.add_argument("-o", action="store", dest="poly_order", required=False, help="Order of polynomial to be used in Savitzky-Golay filter.", default=2)
        parser.add_argument("-w", action="store", dest="window_size", required=False, help="Window length to be used in Savitzky-Golay filter.", default=51)

        plots_avail = ['velocity', 'acceleration']
        args = parser.parse_args()

        # Argument checking section
        try:
            args.path_number = int(args.path_number)
        except ValueError:
            logger.error("Oops!  That was not a valid path number. Integer required.")
            raise ValueError("Oops!  That was no valid path number. Integer required.")

        try:
            args.poly_order = int(args.poly_order)
        except ValueError:
            logger.error("Oops!  That was not a valid polynomial order. Integer required.")
            raise ValueError("Oops!  That was no valid polynomial order. Integer required.")

        try:
            args.window_size = int(args.window_size)
        except ValueError:
            logger.error("Oops!  That was not a valid window size. Integer required.")
            raise ValueError("Oops!  That was no valid window size. Integer required.")

        loaded_data = load_stage_dump(args.filename)
        assert ((args.path_number >= 0) & (args.path_number <= len(loaded_data['spaths']))), "Path number should be in range 0 - " + str(len(loaded_data['spaths'])) #Check if path ID is accessible
        assert (len(loaded_data['spaths'][args.path_number].get_velocity_cont()) >= 3), "Path is too short to be visualized."

        plots_user = args.plot_type.split(',')
        plots = list(set(plots_avail).intersection(plots_user))
        assert (len(plots) > 0), "Oops!  Invalid name(s) of plots provided."

        # Plotting

        divided_paths = map(len, loaded_data['spaths'][args.path_number].paths)
        print divided_paths
        if 0 not in divided_paths:
            indices_incoming = range(divided_paths[0])
            indices_incoming.append(indices_incoming[-1]+1)
            indices_inside = range(divided_paths[0], divided_paths[1]+divided_paths[0])
            indices_inside.append(indices_inside[-1] + 1)
            indices_outgoing = range(divided_paths[1]+divided_paths[0], divided_paths[0]+divided_paths[1]+divided_paths[2])
        else:
            indices_incoming = range(divided_paths[0])
            indices_inside = range(divided_paths[0], divided_paths[1] + divided_paths[0])
            indices_outgoing = range(divided_paths[1] + divided_paths[0], divided_paths[0] + divided_paths[1] + divided_paths[2])

        if 'velocity' in plots:
            velocity = loaded_data['spaths'][args.path_number].get_velocity_cont() #Getting velocity for each frame the path with ID = args.num was detected
            distance = loaded_data['spaths'][args.path_number].get_distance_cont()  # Getting distance travelled for each frame the path with ID = args.num was detected
            try:
                SG_filtered = savgol_filter(x=velocity, window_length=args.window_size, polyorder=args.poly_order)
            except IOError:
                logger.error("Oops!  This path is too short to be visualized.")
                raise ValueError("Oops!  This path is too short to be visualized.")

            plt.plot(distance, velocity, 'k')
            plt.plot(distance[indices_incoming], SG_filtered[indices_incoming], 'r')
            plt.plot(distance[indices_inside], SG_filtered[indices_inside], 'g')
            plt.plot(distance[indices_outgoing], SG_filtered[indices_outgoing], 'b')
            plt.title("Velocity plot")
            plt.xlabel("Distance")
            plt.ylabel("Velocity")
            plt.show()

        if 'acceleration' in plots:
            acceleration = loaded_data['spaths'][args.path_number].get_acceleration_cont() #Getting acceleration for each frame the path with ID = args.num was detected
            distance = loaded_data['spaths'][args.path_number].get_distance_cont()  # Getting distance travelled for each frame the path with ID = args.num was detected
            try:
                SG_filtered = savgol_filter(x=acceleration, window_length=args.window_size, polyorder=args.poly_order)
            except IOError:
                logger.error("Oops!  This path is too short to be visualized.")
                raise ValueError("Oops!  This path is too short to be visualized.")

            plt.plot(distance, acceleration, 'k')
            plt.plot(distance[indices_incoming], SG_filtered[indices_incoming], 'r')
            plt.plot(distance[indices_inside], SG_filtered[indices_inside], 'g')
            plt.plot(distance[indices_outgoing], SG_filtered[indices_outgoing], 'b')
            plt.title("Acceleration plot")
            plt.xlabel("Distance")
            plt.ylabel("Velocity")
            plt.show()


