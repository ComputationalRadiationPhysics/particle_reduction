import argparse
import sys
sys.path.append("../../")
import reduction_main


def multiple_runs_number_conservative(hdf_file_name, hdf_path_result, reduction_percent_start, reduction_percent_end, reduction_percent_step):

    number_of_runs = int((reduction_percent_end -reduction_percent_start)/reduction_percent_step)

    for i in range(0, number_of_runs):
        ratio_of_deleted_particles = reduction_percent_start + i * reduction_percent_step
        name_hdf_file_reduction = hdf_path_result + 'number_conservative_reduction_' + str(ratio_of_deleted_particles)[0:3] + '.h5'
        print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
        reduction_main.number_conservative_thinning_algorithm(hdf_file_name, name_hdf_file_reduction, ratio_of_deleted_particles)


def multiple_runs_energy_conservative(hdf_file_name, hdf_path_result, reduction_percent_start, reduction_percent_end, reduction_percent_step):

    number_of_runs = int((reduction_percent_end -reduction_percent_start)/reduction_percent_step)

    for i in range(0, number_of_runs):
        ratio_of_deleted_particles = reduction_percent_start + i * reduction_percent_step
        name_hdf_file_reduction = hdf_path_result + '/energy_conservative_reduction_' + str(ratio_of_deleted_particles)[0:3] + '.h5'
        print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
        reduction_main.energy_conservative_thinning_algorithm(hdf_file_name, name_hdf_file_reduction, ratio_of_deleted_particles)



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

    if args.algorithm == 'number_conservative':
        multiple_runs_number_conservative(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'energy_conservative':
        multiple_runs_energy_conservative(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)


