import numpy as np
from glob import glob
from pandas import read_csv
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def particle_correlation(data, distance_cutoff=2, depth_cutoff=3, startslice=0, endslice=1, rawstub=""):
    """Take particle positions from multiple slices and return the coordinates of particles common to several slices."""
    particlelist = []

    # Build neighbour matrix
    correlations = cdist(data[:, 0:2], data[:, 0:2]) < distance_cutoff

    # Find particles present in more than x layers
    correlation_sum = (np.sum(correlations, axis=0) > depth_cutoff)
    filteredparticles = data[correlation_sum, 0:2]
    linkagearray = linkage(filteredparticles)
    # plotclusterdistances(linkagearray)
    particleclusters = fcluster(linkagearray, t=distance_cutoff+1, criterion="distance")

    # Collect together particles in each cluster and average their coordinates to get final result
    for i in range(1, max(particleclusters)):
        particlelist.append(np.mean(filteredparticles[particleclusters == i], axis=0))

    if rawstub != "":
        plotresult(filteredparticles, np.array(particlelist), start=startslice, stop=endslice, stub=rawstub)

    return particlelist


def plotclusterdistances(data):
    """Plot the cophenetic distances from a linkage array as a histogram."""
    plt.figure(figsize=(25, 10))
    plt.hist(data[:, 2], bins=100)
    plt.show()


def plotresult(allparticles, finalparticles, plotraw=True, plotallparticles=False, plotselected=True, start=0, stop=1, stub="", filename="output"):
    """Plot the coordinates of identifed particles over the raw images they were obtained from."""

    plot, ax = plt.subplots(figsize=(16, 20))

    if plotraw:
        im = mpimg.imread(stub + str(start) + ".png")
        for image in range(start + 1, stop+1):
            im = im + mpimg.imread(stub + str(image) + ".png")
        im = im/(stop-start+1)
        plt.imshow(im)

    if plotallparticles:
        patches = []
        for particle in range(allparticles.shape[0]):
            patches.append(mpatches.CirclePolygon((allparticles[particle, 0], allparticles[particle, 1]), radius=5))
        p1 = PatchCollection(patches, alpha=0.2, color="blue")
        ax.add_collection(p1)

    if plotselected:
        patches = []
        for particle in range(finalparticles.shape[0]):
            patches.append(mpatches.CirclePolygon((finalparticles[particle, 0], finalparticles[particle, 1]), radius=1))
        p2 = PatchCollection(patches, alpha=1, color="red")
        ax.add_collection(p2)

    plt.xlim(0, 400)
    plt.ylim(0, 500)
    plt.savefig(filename + ".png")


def readinputfile(globstring):
    """
    Read multiple croker/grier style particle files, each representing a z-slice from a z-stack.
    Rerrange the data to my preference and return the data.
    """
    data = []

    # Add slice number to end of data.
    # Data is formatted: x, y, brightness, radius, eccentricity, frame number, particle number, slice number
    for slice_number, inputfile in enumerate(glob(globstring)):
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


def outputparticles(particles):
    """Take the data we have generated and output it as a text file"""
    pass


def main():

    distance_cutoff = 2
    depth_cutoff = 3
    num_timesteps = 2
    startslice = 1
    endslice = 7
    globstring = "slices/slice_[1-7]_raw_coords.txt"
    rawstub = "" #"images/FITC 19_t0001_z000"

    # Begin main
    print("Loading data.")
    data, n_frames, n_slices = readinputfile(globstring)
    print("Data loading complete.")

    for frametime in range(min(n_frames, num_timesteps)):
        outputparticles(particle_correlation(data[frametime], distance_cutoff, depth_cutoff, startslice, endslice, rawstub))

    print "Fin."


main()
