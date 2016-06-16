import numpy as np
from glob import iglob
from pandas import read_csv
from scipy.spatial.distance import cdist
from PIL import Image
import matplotlib.pyplot as plt


def main():
    data = []
    matches = []

    for inputfile in iglob(
            "/Volumes/WIN_DATA/Confocal/STED/Hard spheres/15-12-21/slice_[1-3]_raw_coords.txt"):
        data.append(read_csv(inputfile, sep="\t", header=None).values)

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

    particlepositions = []
    # Get particle positions from neighbour matrix
    for slice1 in range(0, num_slices):
        for particle in range(0, matches[slice].shape[0]):

    # Loop through particles
    #    for slice in range(0, num_slices):
    #       for particle in range(0, matches[slice].shape[0]):
    #            for secondslice in range(slice + 1, num_slices):


    print "Fin."


main()
