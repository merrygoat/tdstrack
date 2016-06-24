import numpy as np
from glob import glob
from pandas import read_csv
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def particle_correlation(data, num_slices):

    distance_cutoff = 2
    depth_cutoff = 3
    # Processes only one timestep
    plot, ax = plt.subplots(figsize=(8, 10))
    print max(data[:, 0])
    print max(data[:, 1])
    # Build neighbour matrix
    correlations = cdist(data[:, 0:2], data[:, 0:2]) < distance_cutoff

    correlation_sum = (np.sum(correlations, axis=0) > depth_cutoff)
    particles = data[correlation_sum, :]

    patches = []
    for particle in range(particles.shape[0]):
        patches.append(mpatches.CirclePolygon((data[particle, 0], data[particle, 1]), radius=5))
    p = PatchCollection(patches, alpha=0.1)

    ax.add_collection(p)
    plt.xlim(0, 400)
    plt.ylim(0, 500)
    #plt.hist(histo, bins=np.arange(10))
    #plt.ylabel("Frequency")
    #plt.xlabel("Number of slices")
    plt.show()

    #plt.imshow(correlations, interpolation="None")
    #plt.savefig("correlations.png", dpi=1000)


def readinputfile():
    data = []

    # Add slice number to end of data.
    # Data is formatted: x, y, brightness, radius, eccentricity, frame number, particle number, slice number
    for slice_number, inputfile in enumerate(glob("slices/slice_[1-7]_raw_coords.txt")):
        frame = read_csv(inputfile, sep="\t", header=None).values
        slice_array = np.full(frame.shape[0], slice_number, dtype=np.int8)[:, np.newaxis]
        frame = np.append(frame, slice_array, axis=1)
        data.append(frame)

    n_frames = int(max(data[0][:, 5]))
    n_slices = len(data)
    sorted_data = [[] for i in range(n_frames)]

    # Data is currently sorted into slices at highest level.
    # Need to reorder to split out seperate frames into highest level list.
    for frame_number in range(n_frames):
        for slice_number in range(n_slices):
            mask = data[slice_number][:, 5] == frame_number
            sorted_data[frame_number].append(data[slice_number][mask, :])

    for frame_number in range(n_frames):
        sorted_data[frame_number] = np.concatenate(sorted_data[frame_number])

    return sorted_data, n_frames, n_slices

def main():

    data, n_frames, n_slices = readinputfile()

    for frametime in range(1):
        particle_correlation(data[frametime], n_frames)

    print "Fin."


main()
