import numpy as np
import sys
import math
import pandas as pd
import TDAtools as TDAtl
from multiprocessing import Pool
from functools import partial

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

# Work out which species dominated
number_of_timesteps = int(tf / dt)
Ntime = math.floor(tf / (image_record_rate * dt))

path = './output/' + str(number_of_rows) + 'x' + str(number_of_columns) + '/grid' + str(grid_option) + '/grazing' + str(
    round(g * 100)) + '/threshold' + str(round(neighborhood_threshold * 100)) + '/'

summary = pd.read_csv(path + 'coral' + str(int(coral_percent * 100)) + '-macro' + str(int(macroalgae_percent * 100)) +
                      '-r' + str(int(r * 10)) + '-d' + str(int(d * 100)) + '-a' + str(int(a * 100)) + '-y' +
                      str(int(y * 100)) + '-time' + str(number_of_timesteps) + '-rec' + str(record_rate) + '-nsim' +
                      str(number_of_simulations) + '.csv', sep=',')

summary = pd.DataFrame(summary).to_numpy()
coralDom = []
macroDom = []
for ii in range(number_of_simulations):

    finalCoral = summary[-1 + Ntime * (ii + 1)][2]
    finalMacro = summary[-1 + Ntime * (ii + 1)][4]

    if finalMacro <= finalCoral:
        coralDom.append(int(ii))
        print('c', ii, finalCoral, finalMacro)
    elif finalMacro > finalCoral:
        macroDom.append(int(ii))
        print('m', ii, finalCoral, finalMacro)

# Extract images
thisImages = TDAtl.extract_images(number_of_simulations, coral_percent, macroalgae_percent, grid_option, number_of_rows,
                                  number_of_columns, neighborhood_threshold, image_record_rate, r, d, a, g, y, dt, tf,
                                  tims)

# Calculate zig-zag
hdim = 0  # can adjust
func = partial(TDAtl.images_to_zz_bars, thisImages, record_rate, dt, hdim, organism, isUnion, group, denoise1, denoise2)
with Pool(number_of_processors) as p:
    ZZ_bars = p.map(func, np.arange(number_of_simulations))

coral_success_list = []
coral_success_list_land = []
macro_success_list = []
macro_success_list_land = []
print(coralDom)
print(macroDom)

for ii in range(number_of_simulations):
    psRbars = ZZ_bars[ii]

    if ii in coralDom:
        # If coral dominated, save the bars in a .txt file, and save the file path for the bars and the landscape 
        np.savetxt(path + 'bars' +
                   '/coral_success/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
                   'coral_success_sim' + str(ii) + '_hdim' + str(hdim) + ".txt", psRbars, delimiter=" ", fmt='%s',
                   header='')
        coral_success_list.append(path + 'bars' + '/coral_success/c' + str(int(coral_percent * 100)) + 'm' +
                                  str(int(macroalgae_percent * 100)) + 'coral_success_sim' + str(ii) + '_hdim' +
                                  str(hdim) + '.txt')
        coral_success_list_land.append(path + 'bars' + '/coral_success/c' + str(int(coral_percent * 100)) + 'm' +
                                       str(int(macroalgae_percent * 100)) +
                                       'coral_success_sim' + str(ii) + '_hdim' + str(hdim) + '.txt.land')
    elif ii in macroDom:
        np.savetxt(path + 'bars' + '/macro_success/c'
                   + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
                   'macro_success_sim' + str(ii) + '_hdim' + str(hdim) + ".txt", psRbars, delimiter=" ", fmt='%s',
                   header='')
        macro_success_list.append(path + 'bars' + '/macro_success/c' + str(int(coral_percent * 100)) + 'm' +
                                  str(int(macroalgae_percent * 100)) + 'macro_success_sim' + str(ii) + '_hdim' +
                                  str(hdim) + '.txt')
        macro_success_list_land.append(path + 'bars' + '/macro_success/c' + str(int(coral_percent * 100)) + 'm' +
                                       str(int(macroalgae_percent * 100)) +
                                       'macro_success_sim' + str(ii) + '_hdim' + str(hdim) + '.txt.land')

# Save .txt files with the list of paths to the bars and landscapes
np.savetxt(path + 'bars' +
           '/coral_success/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'coral_success_all_hdim' + str(hdim) + ".txt", coral_success_list, delimiter=" ", fmt='%s', header='')

np.savetxt(path + 'bars' +
           '/coral_success/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'coral_success_all_land_hdim' + str(hdim) + ".txt", coral_success_list_land, delimiter=" ", fmt='%s',
           header='')

np.savetxt(path + 'bars' +
           '/macro_success/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'macro_success_all_hdim' + str(hdim) + ".txt", macro_success_list, delimiter=" ", fmt='%s', header='')

np.savetxt(path + 'bars' +
           '/macro_success/c' + str(int(coral_percent * 100)) + 'm' + str(int(macroalgae_percent * 100)) +
           'macro_success_all_land_hdim' + str(hdim) + ".txt", macro_success_list_land, delimiter=" ", fmt='%s',
           header='')
