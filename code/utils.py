def check_extension_output(file, ext):
    """ Checks whether chosen output names have the correct extension.

    Parameters
    ----------
    file : str
        File name that should be checked for extension
    ext : str
        Extension the file should have
    Returns
    -------
    Str
        File name adjusted if necessary.
    """

    try:
        if type(file != str):
            file = str(file)

        if type(ext != str):
            ext = str(ext)

        if len(file) <= 0:
            raise NameError

        res = file
        end = file[-(len(ext)):]
        if end != ext:
            res = file + ext

    except NameError:
        print('\nERROR: Name of file needs to be at least one symbol long.\n')
        res = None

    return res


def check_extension_input(file, ext):
    try:
        if type(file != str):
            file = str(file)

        if type(ext != str):
            ext = str(ext)

        if len(file) <= 0 or len(ext) <= 0:
            raise NameError

        res = file
        end = file[-(len(ext)):]
        if end != ext:
            raise AssertionError

    except NameError:
        print('\nERROR: Name for output file needs to be at least one symbol long.\n')
        res = None

    except AssertionError:
        print('\nERROR: Input file: ' + file + ' must have extension: ' + ext)
        res = None

    return res


def cut_extension(name):
    """ Used to cut of extensions of output file names.

    Parameters
    ----------
    name : string
        Name to be shortened.

    Returns
    -------
    string
        The shortened string.
    """

    for i in range(len(name)):
        if name[len(name) - 1] != '.':
            name = name[:-1]
        else:
            name = name[:-1]
            break
    return name
