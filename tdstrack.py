import numpy as np
from glob import glob
from pandas import read_csv
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def particle_correlation(data, distance_cutoff, depth_cutoff, startslice, endslice, rawstub, timestep):
    """Take particle positions from multiple slices and return the coordinates of particles common to several slices."""
    particlelist = []

    # Build neighbour matrix
    correlations = cdist(data[:, 0:2], data[:, 0:2]) < distance_cutoff

    # Find particles present in more than x layers
    correlation_sum = (np.sum(correlations, axis=0) > depth_cutoff)
    filteredparticles = data[correlation_sum, 0:2]
    rejectedparticles = data[np.logical_not(correlation_sum), 0:2]
    linkagearray = linkage(filteredparticles)
    # plotclusterdistances(linkagearray)
    particleclusters = fcluster(linkagearray, t=distance_cutoff+1, criterion="distance")

    # Collect together particles in each cluster and average their coordinates to get final result
    for i in range(1, max(particleclusters)):
        particlelist.append(np.mean(filteredparticles[particleclusters == i], axis=0))

    if rawstub != "":
        plotresult(filteredparticles, rejectedparticles, np.array(particlelist), start=startslice, stop=endslice,
                   stub=rawstub, timestep=timestep, distancecutoff=distance_cutoff)

    return particlelist


def plotclusterdistances(data):
    """Plot the cophenetic distances from a linkage array as a histogram."""
    plt.figure(figsize=(25, 10))
    plt.hist(data[:, 2], bins=100)
    plt.show()


def plotresult(filteredparticles, rejectedparticles, finalparticles, start, stop, stub, timestep, distancecutoff):
    """Plot the coordinates of identifed particles over the raw images they were obtained from."""

    plotraw = True
    plotfilteredparticles = True
    plotrejectedparticles = True
    plotselected = True

    plot, ax = plt.subplots(figsize=(16, 20))

    if plotraw:
        im = mpimg.imread(stub + "t" + '{:04d}'.format(timestep) + "_z" + '{:04d}'.format(start) + ".png")
        for image in range(start + 1, stop+1):
            im = im + mpimg.imread(stub + "t" + '{:04d}'.format(timestep) + "_z" + '{:04d}'.format(image) + ".png")
        im = im/(stop-start+1)
        plt.imshow(im)

    if plotfilteredparticles:
        patches = []
        for particle in range(filteredparticles.shape[0]):
            patches.append(mpatches.CirclePolygon((filteredparticles[particle, 0], filteredparticles[particle, 1]), radius=5))
        p1 = PatchCollection(patches, alpha=0.2, color="blue")
        ax.add_collection(p1)

    if plotrejectedparticles:
        patches = []
        for particle in range(rejectedparticles.shape[0]):
            patches.append(mpatches.CirclePolygon((rejectedparticles[particle, 0], rejectedparticles[particle, 1]), radius=5))
        p1 = PatchCollection(patches, alpha=0.2, color="green")
        ax.add_collection(p1)


    if plotselected:
        patches = []
        for particle in range(finalparticles.shape[0]):
            patches.append(mpatches.CirclePolygon((finalparticles[particle, 0], finalparticles[particle, 1]), radius=1))
        p2 = PatchCollection(patches, alpha=1, color="red")
        ax.add_collection(p2)

    plt.xlim(0, 400)
    plt.ylim(0, 500)
    plt.savefig("t" + str(timestep) + "dist" + str(distancecutoff) + str(plotraw) + str(plotfilteredparticles) +
                str(plotrejectedparticles) + str(plotselected) + ".png")


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

    outputfile = open("processedcoords.txt", 'w')

    for framenumber, frame in enumerate(particles):
        for particle in frame:
            outputfile.write(str(particle[0]) + "\t" + str(particle[1]) + "\t" + str(framenumber) + "\n")


def main():

    distance_cutoff = 2
    depth_cutoff = 3
    num_timesteps = 2
    startslice = 1
    endslice = 7
    globstring = "slices/slice_[1-7]_raw_coords.txt"
    rawstub = "images/FITC 19_"

    # Begin main
    print("Loading data.")
    data, n_frames, n_slices = readinputfile(globstring)
    print("Data loading complete.")

    num_timesteps = min(n_frames, num_timesteps)
    selected_particles = []

    for frametime in range(num_timesteps):
        print("Processing frame " + str(frametime+1) + " of " + str(num_timesteps))
        selected_particles.append(particle_correlation(data[frametime], distance_cutoff, depth_cutoff, startslice, endslice, rawstub, frametime+1))

    outputparticles(selected_particles)
    print "Fin."


main()
