import math
from math import ceil
import bats
import numpy as np
import pandas as pd


def find_all_neighbors(node, grid_shape, threshold):
    """
    Function to find all the neighbors of a point in a finite grid
    Args:
        node ([int, int]): point on a grid
        grid_shape ([int, int]): number of nodes in each dimension of the grid
        threshold (float): neighborhood radius

    Returns:
        all_neighbors (list): list of all neighbors in the grid

    """

    # Initialise list of neighbors
    all_neighbors = []

    # Calculate the max distance between the node and a neighbor
    max_distance = int(np.round(threshold ** 2))

    # Find all neighbors. Note that the node is not its own neighbor
    for ii in range(1 - max_distance, max_distance):
        for jj in range(1 - max_distance, max_distance):
            i = node[0] + ii
            j = node[1] + jj
            if 0 <= i < grid_shape[0] and 0 <= j < grid_shape[1] and \
                    ii ** 2 + jj ** 2 <= max_distance and \
                    not (ii == 0 and jj == 0):
                all_neighbors.append([i, j])

    return all_neighbors


def coral_to_bool(coral_image, organism, group, denoise1=0, denoise2=99):
    """
    A function taking a ternary matrix (with entries 0, 1 or 2) and returning a boolean matrix with the locations of a
    single organism of interest (0, 1 or 2)
    Args:
        coral_image (np.array): ternary two dimensional np array
        organism (int): 0, 1 or 2, denoting which organism is of interest
        group (int): option to group turf (1) with either coral (0) or macroalgae depending on which it is mixed with
        denoise1 (int): option to delete any nodes with fewer than denoise1 neighbors
        denoise2 (int): option to delete any nodes with more than denoise2 neighbors

    Returns:
        coral_bool (np.array): binary two dimensional np array
    """

    shape_image = np.shape(coral_image)

    # Record the species not being selected (macroalgae or coral)
    if organism == 0:
        other = 2
    else:
        other = 0

    # If the coral_image represents a stationary state (either coral or macroalgae extinct), then the boolean matrix is
    # all true
    if group == 1 and sum(sum(1 * np.logical_or(coral_image == organism, coral_image == 1))) == shape_image[0] * \
            shape_image[1]:
        coral_bool = np.full((np.shape(coral_image)), True)

    else:
        coral_bool = np.zeros(shape_image)
        # Group turf nodes with coral/macro if the majority of their neighbours are coral/macro
        if group == 1:
            # Consider each node in turn
            for i in range(shape_image[0]):
                for j in range(shape_image[1]):
                    # If there is turf at that node
                    if coral_image[i, j] == 1:
                        # Find all immediate neighbors of that node
                        neighbour_locations = [[ii, jj] for ii in [i - 1, i, i + 1] for jj in [j - 1, j, j + 1] if
                                               0 <= ii < shape_image[0] and 0 <= jj < shape_image[1] and (
                                                       ii != i or jj != j)]
                        # Record the species at each of the neighboring nodes
                        neighbour_species = [coral_image[neighbour_locations[n][0], neighbour_locations[n][1]] for n in
                                             range(len(neighbour_locations))]
                        # Assign the turf node to organism if there are more neighbors of that type
                        if neighbour_species.count(organism) >= neighbour_species.count(other) and \
                                neighbour_species.count(organism) > 0:
                            coral_bool[i, j] = True

        # Combine reassigned turf nodes with nodes of the organism type
        coral_bool = np.logical_or(coral_bool, coral_image == organism)

        # If denoising is selected, remove some nodes
        if denoise1 >= 1 or denoise2 <= 99:
            conc_profile = binary_to_profile(coral_bool)
            coral_bool = np.logical_and(conc_profile >= denoise1, conc_profile <= denoise2)

    return coral_bool


def bool_to_cubical(coral_bool):
    """
    Function to create a cubical complex from a boolean matrix
    Args:
        coral_bool (np.array): array of boolean (True/False) entries representing coral positions

    Returns:
        coral_complex (bats.CubicalComplex()): cubical complex representing the boolean array

    """

    # Initialise coral complex
    coral_complex = bats.CubicalComplex(2)

    # Cycle through all entries in matrix
    shape_bool = coral_bool.shape
    for i in range(shape_bool[0]):
        for j in range(shape_bool[1]):
            if coral_bool[i][j]:
                # Find which neighbors of the entry are Trues
                adj = np.array([j < shape_bool[1] - 1 and coral_bool[i][j + 1],
                                j < shape_bool[1] - 1 and i < shape_bool[0] - 1 and coral_bool[i + 1, j + 1],
                                i < shape_bool[0] - 1 and coral_bool[i + 1, j]])
                numAdj = np.sum(adj)

                # Add a square if node has three neighbors (left, down, diagdown)
                if numAdj == 3:
                    coral_complex.add_recursive([i, i + 1, j, j + 1])

                else:
                    # Add vertex
                    coral_complex.add_recursive([i, i, j, j])
                    # Add intervals
                    if adj[0]:
                        coral_complex.add_recursive([i, i, j, j + 1])

                    if adj[2]:
                        coral_complex.add_recursive([i, i + 1, j, j])

    return coral_complex


def coral_union(coral_bool_1, coral_bool_2):
    """
    Args:
        coral_bool_1 (np.array): array of boolean values
        coral_bool_2 (np.array): array of boolean values

    Returns:
        (np.array): array of boolean values
    """

    # Returns the union of two binary matrices
    return np.logical_or(coral_bool_1, coral_bool_2)


def coral_intersection(coral_bool_1, coral_bool_2):
    """
    Args:
        coral_bool_1 (np.array): array of boolean values
        coral_bool_2 (np.array): array of boolean values

    Returns:
        (np.array): array of boolean values
    """
    if np.shape(coral_bool_1) != np.shape(coral_bool_2):
        raise Exception("Not same shape")
    # Return the intersection of two binary matrices
    return np.logical_and(coral_bool_1, coral_bool_2)


def binary_to_profile(m, threshold=1.45):
    """
    Function to create a profile matrix based on the number of neighbors each entry of a binary matrix has
    Args:
        m (np.array): binary np.array
        threshold (float): number defining the neighborhood size

    Returns:
        profile (np.array): profile array of number of neighbors
    """

    shape = np.shape(m)
    profile = np.zeros(shape)

    # Cycle through all nodes in the matrix
    for i in range(shape[0]):
        for j in range(shape[1]):

            # If the node contains the species of interest
            if m[i][j]:

                # Cycle through all neighbors of the node
                all_neighbors = find_all_neighbors([i, j], shape, threshold)

                for neighbor in all_neighbors:
                    profile[i][j] += m[neighbor[0], neighbor[1]]

    return profile


def conc_to_ph_diagram(profile, max_conc):
    """

    Args:
        profile: a neighbor profile
        max_conc: the maximum concentration to consider

    Returns:
        D (bats.CubicalComplexDiagram): the PH diagram built from filtering the profile
    """

    D = bats.CubicalComplexDiagram()
    all_concs = np.arange(start=max_conc, stop=-1, step=-1)

    # Add complexes for each concentration step (from 8 to 0)
    for i in range(len(all_concs)):
        this_conc = all_concs[i]
        super_conc = np.copy(profile)
        super_conc[profile < this_conc] = -1
        super_conc = super_conc != -1
        D.add_node(bool_to_cubical(super_conc))

    # Add in all the maps
    for i in range(max_conc):
        D.add_edge(i, i + 1, bats.CubicalMap(D.node_data(i), D.node_data(i + 1)))

    return D


def ph_bars(coral_bool, max_dim=1, max_conc=8):
    """

    Args:
        coral_bool (np.array): a np.array of boolean values representing a coral iterate
        max_dim (int): the maximum dimension of homology to compute
        max_conc (int): the maximum number of neighbors to consider

    Returns:

    """

    # Compute the total number of coral nodes present
    total_benthic = sum((sum(coral_bool)))

    #  If there are no coral nodes, return no bars
    if total_benthic == 0:
        return [[[-1, -1]], [[-1, -1]]]

    # else if there are coral nodes
    else:

        # Get the concentration profile
        conc_profile = binary_to_profile(coral_bool)
        # Get the diagram of complexes
        diagram = conc_to_ph_diagram(2 * conc_profile, 2 * max_conc)

        # Create the chain complex (coeffs in F2)
        chain = bats.Chain(diagram, bats.F2())

        # Initialise the list of bars
        bars = []

        # Compute PH in each dimension
        for hdim in range(max_dim + 1):

            ps = []

            # Compute hdim'th homology
            persistence_module = bats.Hom(chain, hdim)

            # Extract generators (bars)
            ps.extend(bats.barcode_sparse(persistence_module, hdim))

            # Create a list of bars
            bars_hdim = np.zeros((len(ps), 2))

            # Copy bars into the list
            for ii in range(len(ps)):
                if not ps[ii].birth() == ps[ii].death():
                    bars_hdim[ii, 0] = ps[ii].birth()
                    bars_hdim[ii, 1] = ps[ii].death()

            # Sort the list of bars by their persistence (death minus birth)
            bars_hdim = sorted(bars_hdim / 2,
                               key=lambda sub: -(sub[1] - sub[0]))  # sort barcodes from longest to shortest

            # Create a modified list of bars (for visual purposes)
            modified_bars_hdim = []

            # Cycle through the list of (unmodified) bars
            for ii in range(len(bars_hdim)):
                # invert the direction of bars so that they go from large conc to low
                modified_bars_hdim.append([max_conc - bars_hdim[ii][0], max_conc - bars_hdim[ii][1]])

            # Round up the lower end of each bar
            for ii in range(len(modified_bars_hdim)):
                modified_bars_hdim[ii][1] = ceil(modified_bars_hdim[ii][1])

            bars.append(modified_bars_hdim)

    return bars


def ph_barcode(ax, image, fontdict,
               organism=0, group=1, denoise1=1, denoise2=99, max_conc=8, show_legend=False, thickness=2):
    """

    Args:
        ax (matplotlib.axis): existing axis to add barcode to
        image (np.array): ternary matrix representing ABM iterate
        fontdict (dict): parameters for the plot
        organism (int): organism of interest (usually coral)
        group (int): option to group turf nodes with their surrounding nodes
        denoise1 (int): option to ignore nodes with fewer neighbors
        denoise2 (int): option to ignore nodes with more neighbors
        max_conc (int): maximum conc for the filtration
        show_legend (bool): option to include legend below plot
        thickness (int): option to set the thickness of the bars on the plot

    Returns:
        ax with barcode plotted

    """
    xmax = max_conc  # max x tick
    over_tol = 0.15  # amount which the bars overshoot the integer
    bars = ph_bars(coral_to_bool(image, organism, group, denoise1, denoise2), 1, max_conc)
    lstyle = ['solid', 'dotted']
    num_bars = sum([len(bars[i]) for i in range(len(bars))]) + 1

    b = 1
    # cycle through each dimension
    for i in range(len(bars)):

        # cycle through each bar
        for ii in range(len(bars[i])):

            # flip bars (high conc to low conc)
            if bars[i][ii][0] == -1 and bars[i][ii][1] == -1:
                modified_bar = bars[i][ii]
            else:
                modified_bar = [min([max_conc, bars[i][ii][0] + over_tol]),
                                max([0, math.ceil(bars[i][ii][1]) - over_tol])]

            # plot bars
            xxx = np.linspace(modified_bar[0], modified_bar[1], 30)

            if b == 1:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1) / num_bars) + 0 * xxx, linewidth=thickness, color='tab:gray',
                        linestyle=lstyle[i], label="$H_0$")
            elif i == 1 and ii == 0:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1 + 1 * (i == 1)) / num_bars) + 0 * xxx, linewidth=thickness,
                        color='tab:gray', linestyle=lstyle[i], label="$H_1$")
            else:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1 + 1 * (i == 1)) / num_bars) + 0 * xxx, linewidth=thickness,
                        color='tab:gray', linestyle=lstyle[i])
            b = b + 1

    # modify plot
    ax.set_xlim([xmax + 0.25, -0.5])
    ax.set_ylim([0, 1])
    ax.set_xticks(np.linspace(max_conc, 0, max_conc + 1))
    ax.set_xticklabels([str(int(lab)) for lab in np.linspace(max_conc, 0, max_conc + 1)], fontdict=fontdict)

    # option to include legend
    if show_legend:
        ax.legend(bbox_to_anchor=(1, -.1), fontsize=fontdict['fontsize'])
        ax.set_xlabel("direct coral \n neighbors", fontdict=fontdict)
        ax.xaxis.set_label_coords(.275, -.15)

    else:
        ax.set_xlabel("direct coral neighbors", fontdict=fontdict)

    # modify plot
    ax.set_yticks([])
    #ax.set_aspect(aspect=xmax - 0.5, adjustable='box')

    return ax


def extract_images(number_of_simulations, coral_percent, macroalgae_percent, grid_option, number_of_rows,
                   number_of_columns, neighborhood_threshold, record_rate, r, d, a, g, y, dt, tf, tims):
    """
    function to extract images from file
    Args:
        number_of_simulations (int): number of simulations of the ABM
        coral_percent (float): initial coral percentage used in simulations
        macroalgae_percent (float): initial macroalgae percentage used in simulations
        grid_option (int): initial configuration chosen for simulations
        number_of_rows (int): size of grid
        number_of_columns (int): size of grid
        neighborhood_threshold (float): definition of node neighborhood
        record_rate (float): number of timesteps before an image is recorded
        r, d, a, g, y, dt, tf (float): parameters from coral ABM
        tf (int): final time of simulation
        tims (int): final time to pull images (min 0, max tf)

    Returns:
        list of lists of images

    """
    number_of_timesteps = int(round(tf / dt))
    image_data = []

    # read from file
    df = pd.read_csv('./output/' + str(number_of_rows) + 'x' + str(number_of_columns) + '/grid' + str(
        grid_option) + '/grazing' + str(round(g * 100)) + '/threshold' + str(
        int(neighborhood_threshold * 100)) + '/coral' + str(int(coral_percent * 100)) + '-macro' + str(
        int(macroalgae_percent * 100)) + '-r' + str(int(r * 10)) + '-d' + str(int(d * 100)) + '-a' + str(
        int(a * 100)) + '-y' + str(int(y * 100)) + '-time' + str(number_of_timesteps) + '-rec' + str(
        record_rate) + '-nsim' + str(number_of_simulations) + '.csv', sep=',')

    # pull images
    images_df = df['image']

    # distiguish between the total number of images and the images we are interested in analysing
    total_images = math.floor(tf / (record_rate * dt))
    bars_images = math.floor(tims / (record_rate * dt))
    # save each image
    for i in range(number_of_simulations):
        images = []
        for ii in range(bars_images):
            image = images_df[i * total_images + ii]
            image = image[1:-1]  # remove square brackets
            image = np.fromstring(image, sep=" ")
            image = np.reshape(image, (number_of_rows, number_of_columns))
            images.append(image)

        image_data.append(images)

    return image_data


def images_to_zz_bars(simulation_images, record_rate, dt, hdim, species, is_union, group_turf, denoise1, denoise2,
                      simulation=0):
    """
    function that takes images (from many simulations) to extract zigzag bars
    Args:
        simulation_images (list): list of simulations (each entry is a list of images)
        record_rate, dt (int): parameters to work out the timescale
        hdim (int): the dimension of z-z persistence to compute
        species (int): the species of interest (usually coral)
        is_union (int): option to insert the union of iterates (1) or the intersection (0)
        group_turf (int): option to group turf with its surrounding nodes
        denoise1 (int): option to ignore nodes with fewer neighbors
        denoise2 (int): option to ignore nodes with more neighbors
        simulation (int): which simulation to compute bars from

    Returns:

    """
    print('Generating ZZD ', simulation + 1, ', hdim ', hdim)

    images = simulation_images[simulation]

    max_time = len(images)

    d = bats.CubicalComplexDiagram()

    if is_union == 1:

        for i in range(max_time):

            this_rbool = coral_to_bool(images[i], species, group_turf, denoise1, denoise2)

            this_r = bool_to_cubical(this_rbool)

            ii_d = d.add_node(this_r)

            if i != 0:
                # add unions
                d.add_edge(ii_d, ii_d - 1, bats.CubicalMap(this_r, union_r))

            if i != max_time - 1:
                other_rbool = coral_to_bool(images[i + 1], species, group_turf, denoise1, denoise2)

                # add union
                union_rbool = coral_union(this_rbool, other_rbool)

                union_r = bool_to_cubical(union_rbool)

                d.add_node(union_r)

                d.add_edge(ii_d, ii_d + 1, bats.CubicalMap(d.node_data(ii_d), d.node_data(ii_d + 1)))

        f = bats.Chain(d, bats.F2())

    else:

        for i in range(max_time):

            this_rbool = coral_to_bool(images[i], species, group_turf, denoise1, denoise2)

            this_r = bool_to_cubical(this_rbool)

            ii_d = d.add_node(this_r)

            if i != 0:
                # add intersection
                d.add_edge(ii_d - 1, ii_d, bats.CubicalMap(intersection_r, this_r))

            if i != max_time - 1:
                other_rbool = coral_to_bool(images[i + 1], species, group_turf, denoise1, denoise2)

                # add intersection
                intersection_rbool = coral_intersection(this_rbool, other_rbool)

                intersection_r = bool_to_cubical(intersection_rbool)

                d.add_node(intersection_r)

                d.add_edge(ii_d + 1, ii_d, bats.CubicalMap(d.node_data(ii_d + 1), d.node_data(ii_d)))

        f = bats.Chain(d, bats.F2())

        # compute homology
        h = bats.Hom(f, hdim)

        # save bars
        ps = []
        ps.extend(bats.barcode_sparse(h, hdim))
        ps_bars = get_zz_bars(ps, record_rate, dt)

        return ps_bars


def get_zz_bars(ps, record_rate, dt):
    """
    function to extract zig-zag bars from the bats module
    Args:
        ps: persistence object from the bats module
        record_rate, dt (int): parameters to work out the time scale of the zigzag bars

    Returns:
        list of zigzag bars

    """
    bars = []
    for i in range(len(ps)):
        birth = ps[i].birth()
        death = ps[i].death()
        if not birth == death:
            bars.append(np.array([birth, death]) * (record_rate * dt / 2))
    return bars


def zz_barcode(ax, images, fontdict, record_rate, dt, is_union=1, organism=0, group=1, denoise1=1, denoise2=99):
    # Set some parameters
    linewidth = 2  # thickness of bars
    lstyle = ['solid', 'dotted']

    # Extract bars of dim 0 and dim 1
    bars_0 = images_to_zz_bars([images], record_rate, dt, 0, organism, is_union, group, denoise1, denoise2)
    bars_1 = images_to_zz_bars([images], record_rate, dt, 1, organism, is_union, group, denoise1, denoise2)
    bars = [bars_0, bars_1]
    num_bars = sum([len(bars[i]) for i in range(len(bars))]) + 1

    b = 1
    # cycle through bars of all dimension
    for i in range(len(bars)):

        # cycle through all bars
        for ii in range(len(bars[i])):

            # plot bars
            xxx = np.linspace(bars[i][ii][0], bars[i][ii][1], 30)

            if b == 1:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1) / num_bars) + 0 * xxx, linewidth=linewidth, color='r',
                        linestyle=lstyle[i], label="$H_0$")
            elif i == 1 and ii == 0:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1 + 1 * (i == 1)) / num_bars) + 0 * xxx, linewidth=linewidth,
                        color='r', linestyle=lstyle[i], label="$H_1$")
            else:
                ax.plot(xxx, 0.95 - 0.95 * ((b - 1 + 1 * (i == 1)) / num_bars) + 0 * xxx, linewidth=linewidth,
                        color='r', linestyle=lstyle[i])
            b = b + 1

    # modify axes
    ax.set_xlim([0, (-1 + len(images)) * record_rate * dt])
    ax.set_ylim([0, 1])
    ax.set_xticks(np.linspace(0, (-1 + len(images)) * record_rate * dt, len(images)))
    ax.legend(loc='lower right', fontsize=fontdict['fontsize'])
    ax.set_yticks([])
    ax.set_xlabel("$t$", fontsize=fontdict['fontsize'])

    return ax


def get_landscape_critical_points(landscape_filepath):
    """
    function to extract critical points (the peaks) from a landscape file
    Args:
        landscape_filepath (str): path on local drive to landscape file

    Returns:
        list of landscape critical points

    """
    lan_file = open(landscape_filepath, 'r')
    landscapes = []
    final_landscape = False
    LS = 0
    while not final_landscape:
        line = lan_file.readline()
        if len(line) == 0:
            final_landscape = True
            continue
        if LS == 0:
            LS = LS + 1
            this_landscape = []
        elif line[0] == "#":
            landscapes.append(this_landscape)
            LS = LS + 1
            this_landscape = []
        else:
            this_landscape.append(np.fromstring(line, dtype=float, sep='  '))

    return landscapes


def cp_landscape_to_pd(landscape_cps):
    """
    function which takes a list of landscape critical points and converts this to an equivalent persistence diagram
    Args:
        landscape_cps (list): list of landscape critical points

    Returns:
        list of persistence diagram points

    """
    PD_CPs = []
    for landscape_index in landscape_cps:
        landscape_i = []
        for landscape in landscape_index:
            landscape_i.append([landscape[0] - landscape[1], landscape[0] + landscape[1]])
        PD_CPs.append(landscape_i)

    return PD_CPs


def integrate_landscape(landscape_cps):
    """
    function to integrate a landscape based on its critical points
    Args:
        landscape_cps (list): list of landscape critical points

    Returns:
        float of integral of landscape

    """
    integral = 0
    for cp in range(len(landscape_cps) - 1):
        cp1 = landscape_cps[cp]
        cp2 = landscape_cps[cp + 1]
        # integrate using the trapezium rule (noting that first and last critical points are always (0, 0))
        integral = integral + ((cp1[1] + cp2[1]) / 2) * (cp2[0] - cp1[0])

    return integral
