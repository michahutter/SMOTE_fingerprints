from source.code import utils


def create_ranking_output(ranks, out_file):
    """ Creates a file containing the active compounds sorted by their summed Manhatten distances to the inactives.

    Parameters
    ----------
    ranks : list of list
        Ranking created by the ranking_manhatten() function
    out_file : str
        Name of the output file

    Returns
    -------
    list of list
        Each list contains a list, representing one entry row in the output data. One list contains
        the active's name and the corresponding Manhatten distance sum in that order.

    """
    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return None

    data = list()
    # First line of output file
    s = 'name distance_sum\n'
    s = s.encode('ascii')
    with open(out_file, 'wb') as f:
        f.write(s)
        # Generated output respecting first line
        for i in range(len(ranks)):
            line = ranks[i][0] + ' ' + str(ranks[i][1]) + '\n'
            line = line.encode('ascii')
            data.append([ranks[i][0], ranks[i][1]])
            f.write(line)

    f.close()
    return data


def ranking_manhatten(data, highest=True):
    """ Sorts the actives by their summed up Manhatten distances to the inactives.

    Parameters
    ----------
    data : list of list
        Each list contains a list with the data found in one row of a distance file. One list contains the name of an
        active, name of an inactive and the Manhatten distance of this pair, in this order. Data is provided by
        distances.read_distances() function.
    highest : bool
        True if ranking should be in descending order, false for ascending.

    Returns
    -------
    list of list
        Every list contains a list with the name of an active, and it's summed up Manhatten distances to all inactives.
        the outer lists are sorted by the sum of the Manhatten distance.
    """

    groups = group_manhatten(data)
    for i in range(len(groups)):
        groups[i][1] = sum(groups[i][1])

    # Highest value on top
    if highest:
        groups.sort(key=lambda x: x[1], reverse=True)

    # Lowest value on top
    else:
        groups.sort(key=lambda x: x[1], reverse=False)

    return groups


def group_manhatten(data):
    """ Groups all Manhatten distance entries in a distance file belonging to one active and sums the distances.

    Helper function for ranking_manhatten()

    Parameters
    ----------
    data : list of list
        Each list contains a list with the data found in one row of a distance file. One list contains the name of an
        active, name of an inactive and the Manhatten distance of this pair, in this order. Data is provided by
        distances.read_distances() function.

    Returns
    -------
    list of list
         Every list contains a list with the actives name and the summed manhattan distance, in this order.
    """

    active = data[0][0]
    man = list()
    sums = list()
    for i in range(len(data)):
        if data[i][0] != active:
            sums.append([active, man])
            active = data[i][0]
            man = list()
        man = man + [data[i][2]]

    sums.append([active, man])
    return sums


# ***************************************** Creation of combining files *****************************************


def read_ranking_output(rank_file):
    """ Reads generated rank files.

    Files generated with the ranks command/create_ranking_output function can be read with this function.

    Parameters
    ----------
    rank_file : str
        Path to the ranks file

    Returns
    -------
    list of list
        Every list contains a list with the actives name and the summed manhattan distance, in this order.
    """

    try:
        if type(rank_file) != str:
            raise TypeError

        data = list()
        with open(rank_file, 'r') as in_file:
            next(in_file)
            while True:

                line = in_file.readline()
                if not line:
                    break

                llist = line.split()
                data.append(llist)

    except FileNotFoundError:
        print('\nERROR: Input ranks_file: ' + str(rank_file) + ' could not be found.\n')
        data = None

    except TypeError:
        print('\nERROR: The function read_ranking_output got non str input.\n')
        data = None

    return data
