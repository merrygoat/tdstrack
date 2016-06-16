import numpy as np
from glob import iglob
from pandas import read_csv
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def slicecorrelation(data):

    matches = []
    num_slices = len(data)

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

    for inputfile in iglob("slices/slice_[1-7]_raw_coords.txt"):
        data.append(read_csv(inputfile, sep="\t", header=None).values)
    # Once appending is done then condense python list into np array


    return data

def initialise_cell_list(data, num_slices):

    max_x = 0

    #sort data by x-coordinate
    for z_slice in range(0, num_slices):
        data = data[z_slice][np.lexsort(data[:][:, 0]), :]
        newmax_x = np.max(data[:][:, 0])-np.min(data[:][:,0])
        if newmax_x > max_x:
            max_x = newmax_x

    return data


def main():

    data = readinputfile()
    num_slices = slicecorrelation(data)

    data = initialise_cell_list(data, num_slices)

    print "Fin."


main()
