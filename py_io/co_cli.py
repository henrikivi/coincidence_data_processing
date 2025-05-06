"""Functions for running code with a CLI
"""

import argparse as arg


def parsing():
    """Argument Parser for coincidence analysis

    Returns
    -------
    parse
        The argument parser object
    """
    parse = arg.ArgumentParser(
        description="Purpose: Take coincidence data for multiple ionization and allow peaks to be analysed"
    )

    parse.add_argument(
        "-filepath",
        metavar="FILEPATH",
        type=str,
        help="Use when analysing a single experiment: Path to file (without extension) to load result from",
    )

    parse.add_argument(
        "-param_file",
        metavar="PARAM_FILE",
        type=str,
        help="Name of csv-file (without extension) containing peak bound coordinates"
    )

    parse.add_argument(
        "-data_folder",
        metavar="DATA_FOLDER",
        type=str,
        help="Use when analysing whole dataset: Folder containing entire EIEI dataset"
    )

    return parse
