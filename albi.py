#!/usr/bin/python

import os
import sys
import argparse
from numpy import arange, loadtxt, savetxt, hstack, vstack, newaxis, concatenate, array
import albi_lib
import progress_bar

class AlbiArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

albi_USAGE  = \
"""
If you want to build_heritability_cis from a kinship_eigenvalues please specify        %(prog)s --kinship_eigenvalues  [kinship file]        --estimates_filename/--estimate_grid
If you want to build_heritability_cis from a distribution file please specify          %(prog)s --load_distributions   [distribuation file]  --estimates_filename/--estimate_grid
If you just want to calculate_probability_intervals please specify                     %(prog)s  --kinship_eigenvalues [kinship file]        --save_distributions [distribuation file]

There are optional arguments you can specify as you can see below:
"""

# opt2
def estimate_distributions(h2_values, H2_values, kinship_eigenvalues_data, samples = 1000, distributions_filename = None):
    print("Estimating distributions...")
    distributions = albi_lib.calculate_probability_intervals(h2_values = h2_values, 
                                                             H2_values = H2_values, 
                                                             kinship_eigenvalues = kinship_eigenvalues_data, 
                                                             n_random_samples = samples)
    
    print("Done estimating distributions.")
    if distributions_filename:
        header = ['ALBI'] + ['0'] + ["%f-%f" % (H2_values[i], H2_values[i+1]) for i in range(len(H2_values)-1)] + ['1']
        savetxt(fname = distributions_filename, 
                X = hstack([h2_values[:, newaxis], distributions]),
                fmt = "%1.5f",
                delimiter = '\t',
                header = '\t'.join(header),
                comments = '')
        print("Distributions are written to '%s'." % distributions_filename)
    return distributions

# opt3 
def build_heritability_cis_from_distributions(h2_values, H2_values, all_distributions_data, estimates, confidence = 0.95, output_filename = None):
    print("Building heritability CIs...")
    cis = albi_lib.build_heritability_cis(h2_values = H2_values, 
                                          H2_values = H2_values,
                                          all_distributions = all_distributions_data, 
                                          estimated_values = estimates, 
                                          confidence = confidence, 
                                          use_randomized_cis = False)
    print("Done building CIs.")
    if output_filename:
        header = ['Estimate', 'CI_lower_bound', 'CI_upper_bound']
        savetxt(fname = output_filename, 
                X = hstack([estimates[:, newaxis], cis]),
                fmt = "%1.5f",
                delimiter = '\t',
                header = '\t'.join(header),
                comments = '')
        print("Heritability CIs are written to '%s'." % output_filename)
    return cis

# opt1
def build_heritability_cis_from_kinship(h2_values, H2_values, kinship_eigenvalues_data, estimates, confidence = 0.95, samples = 1000, distributions_filename = None, output_filename = None):
    distribution = estimate_distributions(h2_values = h2_values,
                                                     H2_values = H2_values,
                                                     kinship_eigenvalues_data = kinship_eigenvalues_data,
                                                     samples = samples,
                                                     distributions_filename = distributions_filename
                                                   )

    return build_heritability_cis_from_distributions(h2_values = h2_values,
                                                      H2_values = H2_values, 
                                                      all_distributions_data = distribution, 
                                                      estimates = estimates, 
                                                      confidence = confidence, 
                                                      output_filename = output_filename
                                                     )
    print("Done")


def _get_estimates(estimate_grid, estimates_filename):
    if estimate_grid:
        estimates = arange(0, 1 + 1.0/estimate_grid, 1.0/estimate_grid)
    elif not os.path.exists(estimates_filename) :
        print("The file '%s' doesn't exist. Exiting" % estimates_filename)
        sys.exit(2)  
    else:
        estimates = file(args.estimates_filename, 'rb').read()
    return estimates


def run_albi(kinship_eigenvalues_filename = None,
              estimate_grid = None,
              estimates_filename = None,
              save_distributions_filename = None,
              load_distributions_filename = None,
              precision_h2 = 0.01,
              precision_H2 = 0.01,
              samples = 1000,
              confidence = 0.95,
              output_filename = None):
    """
    There are three options to run albi:
    opt1 - build_heritability_cis_from_kinship. The full run of albi. gets a kinship eigenvalues and returns the heritability CIs. (This is actually opt2 + opt3)
    opt2 - calculate_probability_intervals. gets a kinship eigenvalues and calculate probability intervals (distributions) and saves it for later.
    opt3 - build_heritability_cis_from_distributions. gets the distributions file from opt2 and returns the heritability CIs.

    This function parses the script arguments and runs the right option according to the specified arguments:
    for opt1 - you have to specify kinship_eigenvalues + estimates
    for opt2 - you have to specify kinship_eigenvalues + save_distributions
    for opt3 - you have to specify load_distributions + estimates

    """

    # if a kinship_eigenvalues wasn't speified - the user wants to run opt3 or he was wrong. If it was specified the user wants to run opt1 or opt2
    if kinship_eigenvalues_filename is None: # opt3 or error
        if load_distributions_filename is None:
            print("If you want to build_heritability_cis from a kinship_eigenvalues please specify        --kinship_eigenvalues  [kinship file]        --estimates_filename/--estimate_grid\n" \
                + "If you want to build_heritability_cis from a distribution file please specify          --load_distributions   [distribuation file]  --estimates_filename/--estimate_grid\n" \
                + "If you just want to calculate_probability_intervals please specify                     --kinship_eigenvalues  [kinship file]        --save_distributions [distribuation file]" 
               )
            return None
        elif not os.path.exists(load_distributions_filename):
            print("The file '%s' doesn't exist. Exiting" % load_distributions_filename) #use logging ifo?TODO
            return None
        elif estimates_filename is None and estimate_grid is None:
            print("Few arguments are missing.\n" \
                  + "If you want to build_heritability_cis from a distribution file please also specify --estimate_grid/--estimates_filename"
                 )
            return None
        else:
            # run opt3
            estimates = _get_estimates(estimate_grid, estimates_filename)

            # Read distributions file
            all_distributions_data = loadtxt(load_distributions_filename, delimiter='\t', skiprows=1)
            header = file(load_distributions_filename, 'rb').readline().strip().split()[1:]
            
            h2_values = all_distributions_data[:,0]
            H2_values = array([float(header[0])] + [float(field_name.split('-')[1]) for field_name in header[1:-2]] + [float(header[-1])])

            return build_heritability_cis_from_distributions(h2_values = h2_values,
                                                             H2_values = H2_values,
                                                             all_distributions_data = all_distributions_data[:,1:], 
                                                             estimates = estimates,
                                                             confidence = confidence,
                                                             output_filename = output_filename)

    elif not os.path.exists(kinship_eigenvalues_filename):
        print("The file '%s' doesn't exist. Exiting" % kinship_eigenvalues_filename)
        return None

    else: # opt1 or opt2
        if estimates_filename is None and estimate_grid is None: #opt2 or error
            if save_distributions_filename is None:
                print("Few arguments are missing.\n" \
                    + "If you want to build_heritability_cis please also specify --estimate_grid/--estimates_filename.\n" \
                    + "If you want to calculate_probability_intervals please also specify --save_distributions file"
                   )
                return None
            else:
                # run opt2
                return estimate_distributions(h2_values = arange(0, 1 + precision_h2, precision_h2),
                                              H2_values = arange(0, 1 + precision_H2, precision_H2),
                                              kinship_eigenvalues_data = loadtxt(kinship_eigenvalues_filename),
                                              samples = progress_bar.ProgressBarIter(samples),
                                              distributions_filename = save_distributions_filename
                                              )
        else:
            # run opt1
            estimates = _get_estimates(estimate_grid, estimates_filename)
            return build_heritability_cis_from_kinship(h2_values = arange(0, 1 + precision_h2, precision_h2),
                                                       H2_values = arange(0, 1 + precision_H2, precision_H2),
                                                       kinship_eigenvalues_data = loadtxt(kinship_eigenvalues_filename),
                                                       estimates = estimates, 
                                                       confidence = confidence,
                                                       samples = progress_bar.ProgressBarIter(samples),
                                                       distributions_filename = save_distributions_filename,
                                                       output_filename = output_filename


if __name__ == '__main__':
    # Parse arguments
    parser = AlbiArgumentParser(prog=os.path.basename(sys.argv[0]),  usage=albi_USAGE)
        
    group_load = parser.add_mutually_exclusive_group(required = True)
    group_load.add_argument('-l', '--load_dist_filename',                     type = str,                                   help = "default is None. This is a filename to which the program will load the output of calculate_probability_intervals, which is a matrix (loadtxt). This flag is mutually exclusive with --input")
    group_load.add_argument('-k', '--kinship_eigenvalues', help="Filename of a file containing the eigenvalues of the kinship matrix") 
    parser.add_argument('-p', '--precision',             nargs = '?',     type = float,   default = 0.01,    help = "default value is 0.1.  the intervals between the h^2 values albi checks")
    parser.add_argument('-d', '--distribution_precision',            nargs = '?',     type = float,   default = 0.01,    help = "default value is 0.1.  I dont know what this is for :p") 
    parser.add_argument('-n', '--samples',               nargs = '?',     type = int,     default = 1000,    help = "default value is 1000. Number of samples albi should take") #   monte_carlo_size
    parser.add_argument('-s', '--save_dist_filename',                     type = str,                                   help = "default is None. This is a filename to which the program will save the output of calculate_probability_intervals, which is a matrix (take a look at savetxt).")
    
    group_estimates = parser.add_mutually_exclusive_group(required = False)
    group_estimates.add_argument('-f', '--estimates_filename',                     type = str,                                   help = "A filename for the heritability estimates for which we want to build CIs. A text file with one estimate per row.")
    group_estimates.add_argument('-g', '--estimate_grid',                          type = int,                                 help = "How many 'jumps' to ahve between 0 and 1")
    parser.add_argument('-c', '--confidence',            nargs = '?',     type = float,   default = 0.95,    help = "default value is 0.95. How exact should the resualt CI be (percentage between 0-100")
    parser.add_argument('-o', '--output_filename',                        type = str,                                   help = "default is None. The filename to which we will output the results of albi. A text file of a matrix, where each row is (h^2 estimate, lower bound, upper bound")

    group_verbose = parser.add_mutually_exclusive_group(required = False)
    group_verbose.add_argument('-q', '--quiet', default = False, action="store_true", help = "Do not print progress information to screen.")
    
    args = parser.parse_args()

    # Validate arguments
    if args.kinship_eigenvalues is None and args.estimates_filename is None and args.estimate_grid is None:
        parser.print_help()
        sys.exit(2)
    
    # Run!
    run_albi(kinship_eigenvalues_filename = args.kinship_eigenvalues,
             estimate_grid = args.estimate_grid, 
             estimates_filename = args.estimates_filename,
             save_distributions_filename = args.save_dist_filename,
             load_distributions_filename = args.load_dist_filename,
             precision_h2 = args.precision,
             precision_H2 = args.distribution_precision,
             samples = args.samples,
             confidence = args.confidence,
             output_filename = args.output_filename)



