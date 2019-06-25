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


def multiple_runs_random(hdf_file_name, hdf_path_result, reduction_percent_start, reduction_percent_end, reduction_percent_step):

    number_of_runs = int((reduction_percent_end - reduction_percent_start)/reduction_percent_step)

    for i in range(0, number_of_runs):
        ratio_of_deleted_particles = reduction_percent_start + i * reduction_percent_step
        name_hdf_file_reduction = hdf_path_result + '/random_reduction_' + str(ratio_of_deleted_particles)[0:3] + '.h5'
        print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
        reduction_main.random_thinning_algorithm(hdf_file_name, name_hdf_file_reduction, ratio_of_deleted_particles)


def multiple_runs_k_means_clustering(hdf_file_name, hdf_path_result, reduction_percent_start, reduction_percent_end, reduction_percent_step):

    number_of_runs = int((reduction_percent_end -reduction_percent_start)/reduction_percent_step)

    for i in range(0, number_of_runs):
        ratio_of_deleted_particles = reduction_percent_start + i * reduction_percent_step
        name_hdf_file_reduction = hdf_path_result + '/k_means_clustering_' + str(ratio_of_deleted_particles)[0:3] + '.h5'
        print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
        reduction_main.k_means_cluster_algorithm(hdf_file_name, name_hdf_file_reduction, ratio_of_deleted_particles)


def multiple_runs_k_means_avg(hdf_file_name, hdf_path_result, reduction_percent_start, reduction_percent_end, reduction_percent_step):

    number_of_runs = int((reduction_percent_end - reduction_percent_start)/reduction_percent_step)

    for i in range(0, number_of_runs):
        ratio_of_deleted_particles = reduction_percent_start + i * reduction_percent_step
        name_hdf_file_reduction = hdf_path_result + '/k_means_avg_' + str(ratio_of_deleted_particles)[0:3] + '.h5'
        print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
        reduction_main.k_means_avg_algorithm(hdf_file_name, name_hdf_file_reduction, ratio_of_deleted_particles)


def multiple_runs_voronoi(hdf_file_name, hdf_path_result, tolerance_momentum_start,
                          tolerance_momentum_end, tolerance_momentum_step,
                          tolerance_position_start,
                          tolerance_position_end, tolerance_position_step):

    number_of_runs_momentum = int((tolerance_momentum_end - tolerance_momentum_start)/tolerance_momentum_step)
    number_of_runs_position = int((tolerance_position_end - tolerance_position_start) / tolerance_position_step)

    for i in range(0, number_of_runs_momentum):
        for j in range(0, number_of_runs_position):
            tolerance_momentum = tolerance_momentum_start + i * tolerance_momentum_step
            tolerance_position = tolerance_position_start + j * tolerance_position_step
            name_hdf_file_reduction = hdf_path_result + '/voronoi_'+ str(tolerance_momentum) + str(tolerance_position) + '.h5'
            print('name_hdf_file_reduction  ' + str(name_hdf_file_reduction))
            reduction_main.voronoi_algorithm(hdf_file_name, name_hdf_file_reduction, tolerance_momentum, tolerance_position)


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



    ## for voronoi

    parser.add_argument("-reduction_momentum_start", metavar='reduction_momentum_start', type=float,
                        help="starting reduction momentum")

    parser.add_argument("-reduction_momentum_end", metavar='reduction_momentum_end', type=float,
                        help="ending reduction momentum")

    parser.add_argument("-reduction_momentum_step", metavar='reduction_momentum_step', type=float,
                        help="step of reduction momentum")

    ## for position

    parser.add_argument("-reduction_position_start", metavar='reduction_position_start', type=float,
                        help="starting reduction position")

    parser.add_argument("-reduction_position_end", metavar='reduction_position_end', type=float,
                        help="ending reduction position")

    parser.add_argument("-reduction_position_step", metavar='reduction_position_step', type=float,
                        help="step of reduction position")


    args = parser.parse_args()

    if args.algorithm == 'number_conservative':
        multiple_runs_number_conservative(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'energy_conservative':
        multiple_runs_energy_conservative(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'random':
        multiple_runs_random(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'k_means_clustering':
        multiple_runs_k_means_clustering(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'k_means_avg':
        multiple_runs_k_means_avg(args.hdf, args.hdf_re, args.reduction_percent_start,
                                          args.reduction_percent_end, args.reduction_percent_step)

    elif args.algorithm == 'voronoi':
        print(args.reduction_momentum_start)
        print(args.reduction_momentum_end)
        print(args.reduction_momentum_step)

        print(args.reduction_position_start)
        print(args.reduction_position_end)
        print(args.reduction_position_step)
        multiple_runs_voronoi(args.hdf, args.hdf_re, args.reduction_momentum_start,
                                          args.reduction_momentum_end, args.reduction_momentum_step,
                              args.reduction_position_start, args.reduction_position_end, args.reduction_position_step)


