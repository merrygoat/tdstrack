import numpy as np
from glob import glob
from pandas import read_csv
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def slice_correlation(data, num_slices):

    matches = []

    plt.figure(figsize=(10, 10))

    # Build neighbour matrix
    for framenumber in range(0, 1):
        for slice1 in xrange(0, num_slices):
            for slice2 in range(slice1 + 1, num_slices):
                matches.append((cdist(data[slice1][data[slice1][:, 5] == framenumber, 0:2],
                                      data[slice2][data[slice2][:, 5] == framenumber, 0:2]) < 1))
                print framenumber, slice1, slice2

                plt.imshow(matches[len(matches) - 1])
                plt.savefig(str(slice1) + "_" + str(slice2) + ".png", dpi=100)

    return num_slices

def readinputfile():
    data = []

    # Add slice number to end of frame.
    for slice_number, inputfile in enumerate(glob("slices/slice_[1-7]_raw_coords.txt")):
        frame = read_csv(inputfile, sep="\t", header=None).values
        slice_array = np.full(frame.shape[0], slice_number)[:, np.newaxis]
        frame = np.append(frame, slice_array, axis=1)

        data.append(frame)

    # Need to split out sperate frames into highest level list.

    data = np.concatenate(data[:])



    return data

def initialise_cell_list(data, num_slices):

    max_x = 0
    cell_width = 2

    #sort data by x-coordinate and find max x coodinate
    for z_slice in range(0, num_slices):
        data[z_slice] = data[z_slice][np.argsort(data[z_slice][:, 0]), :]
        newmax_x = np.max(data[z_slice][:, 0])
        if newmax_x > max_x:
            max_x = newmax_x

    n_cells = np.ceil(max_x / cell_width)

    for z_slice in range(0, num_slices):
        histo, binedges = np.histogram(data[z_slice][:, 0], bins=n_cells, range=[0,n_cells*cell_width])
        # Need to convert from frequency to cumulative frequency for using in np.split
        histo = np.cumsum(histo)
        data[z_slice] = np.split(data[z_slice], histo)

    return data, n_cells

def particle_correlation(data, n_slices, n_cells):

    pass

def main():

    data = readinputfile()
    n_slices = len(data)

    data, n_cells = initialise_cell_list(data, n_slices)

    # slice_correlation(data, num_slices)

    print "Fin."


main()
