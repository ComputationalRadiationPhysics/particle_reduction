import argparse


if __name__ == "__main__":
    """ Parse arguments from command line """

    parser = argparse.ArgumentParser(description="voronoi reduction")

    parser.add_argument("-hdf", metavar='hdf_file', type=str,
                        help="hdf file without patches")

    parser.add_argument("-hdf_re", metavar='hdf_file_reduction', type=str,
                        help="reducted hdf file")

    parser.add_argument("-algorithm", metavar='algorithm', type=str,
                        help="algorithm")

    parser.add_argument("-reduction_percent_start", metavar='reduction_percent_start', type=float,
                        help="starting reduction percent")

    parser.add_argument("-reduction_percent_end", metavar='reduction_percent_end', type=float,
                        help="ending reduction percent")

    parser.add_argument("-reduction_percent_step", metavar='reduction_percent_step', type=float,
                        help="step of reduction percent")


    args = parser.parse_args()
