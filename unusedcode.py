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