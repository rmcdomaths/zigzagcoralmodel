import numpy as np
import sys
import math
import pandas as pd
from multiprocessing import Pool
from functools import partial
import TDAtools as TDAtl

# Parameters

number_of_processors = int(sys.argv[1])
number_of_simulations = int(sys.argv[2])

coral_percent = float(sys.argv[3]) / 100
macroalgae_percent = float(sys.argv[4]) / 100
turf_percent = round(1 - coral_percent - macroalgae_percent, 2)

grid_option = int(sys.argv[5])
number_of_rows, number_of_columns = int(sys.argv[6]), int(sys.argv[7])
neighborhood_threshold = float(sys.argv[8])
record_rate = int(sys.argv[9])
image_return = str(sys.argv[10])
image_record_rate = int(sys.argv[11])
image_counter = image_record_rate
r, d, a, g, y, dt, tf = [float(sys.argv[n]) for n in range(12, 19)]

# Extra parameters
organism = int(sys.argv[19])
group = int(sys.argv[20])
denoise1 = int(sys.argv[21])
denoise2 = int(sys.argv[22])
isUnion = int(sys.argv[23])
tims = int(sys.argv[24])
number_of_timesteps = int(tf / dt)

# Main code

path = './output/' + str(number_of_rows) + 'x' + str(number_of_columns) + '/grid' + str(grid_option) + '/grazing' + str(
    round(100 * g)) + '/threshold' + str(int(neighborhood_threshold * 100)) + '/'

summary = pd.read_csv(
    path + 'coral' + str(int(coral_percent * 100)) + '-macro' + str(int(macroalgae_percent * 100)) + '-r' + str(
        int(r * 10)) + '-d' + str(int(d * 100)) + '-a' + str(int(a * 100)) + '-y' + str(int(y * 100)) + '-time' + str(
        number_of_timesteps) + '-rec' + str(record_rate) + '-nsim' + str(number_of_simulations) + '.csv', sep=',')

bars_images = math.floor(tims / (record_rate * dt))

# Extract images
thisImages = TDAtl.extract_images(number_of_simulations, coral_percent, macroalgae_percent, grid_option, number_of_rows,
                                  number_of_columns, neighborhood_threshold, record_rate, r, d, a, g, y, dt, tf, tims)

# Calculate zig-zag bars
hdim = 0  # can adjust
func = partial(TDAtl.images_to_zz_bars, thisImages, record_rate, dt, hdim, organism, isUnion, group, denoise1, denoise2)
with Pool(number_of_processors) as p:
    ZZ_bars = p.map(func, np.arange(number_of_simulations))

bar_list = []
land_list = []
path = './output/' + str(number_of_rows) + 'x' + str(number_of_columns) + '/grid' + str(grid_option) + '/grazing' + str(
    round(100 * g)) + '/threshold' + str(int(neighborhood_threshold * 100)) + '/'

for ii in range(number_of_simulations):

    psRbars = ZZ_bars[ii]

    np.savetxt(path + 'bars' + '/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
               'all_bars_sim' + str(ii) + '_hdim' + str(hdim) + '.txt', psRbars, delimiter=" ", fmt='%s', header='')
    bar_list.append(path + 'bars' + '/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
                    'all_bars_sim' + str(ii) + '_hdim' + str(hdim) + '.txt')
    land_list.append(path + 'bars' + '/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
                     'all_bars_sim' + str(ii) + '_hdim' + str(hdim) + '.txt.land')

np.savetxt(path + 'bars' + '/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'all_bars_list_hdim' + str(hdim) + ".txt", bar_list, delimiter=" ", fmt='%s', header='')

np.savetxt(path + 'bars' + '/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'all_lands_list_hdim' + str(hdim) + ".txt", land_list, delimiter=" ", fmt='%s', header='')
