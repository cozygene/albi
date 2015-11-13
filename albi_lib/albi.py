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

albi_USAGE = ""
"""
If you want to build_heritability_cis from a kinship_eigenvalues please specify        %(prog)s --kinship_eigenvalues  [kinship file]        --estimates_filename/--estimate_grid
If you want to build_heritability_cis from a distribution file please specify          %(prog)s --load_distributions   [distribuation file]  --estimates_filename/--estimate_grid
If you just want to calculate_probability_intervals please specify                     %(prog)s  --kinship_eigenvalues [kinship file]        --save_distributions [distribuation file]

There are optional arguments you can specify as you can see below:
"""

# opt2
def estimate_distributions(h2_values, H2_values, kinship_eigenvalues_data, samples = 1000, distributions_filename = None):
    print("Estimating distributions...")
    distributions = albi_lib.estimate_distributions(h2_values = h2_values, 
                                                             H2_values = H2_values, 
                                                             kinship_eigenvalues = kinship_eigenvalues_data, 
                                                             n_random_samples = samples)
    
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
    cis = albi_lib.build_heritability_cis(h2_values = h2_values, 
                                          H2_values = H2_values,
                                          all_distributions = all_distributions_data, 
                                          estimated_values = estimates, 
                                          confidence = confidence, 
                                          use_randomized_cis = False)
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

    if kinship_eigenvalues_filename is None:
        if load_distributions_filename is None:
            # TODO: Add error
            return None
        elif not os.path.exists(load_distributions_filename):
            print("The file '%s' doesn't exist. Exiting" % load_distributions_filename) #use logging ifo?TODO
            return None
        elif estimates_filename is None and estimate_grid is None:
            # TODO: Add error
            return None
        else:            
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

    else:
        if estimates_filename is None and estimate_grid is None: #opt2 or error
            if save_distributions_filename is None:
                # TODO: Add error
                return None
            else:
                return estimate_distributions(h2_values = arange(0, 1 + precision_h2, precision_h2),
                                              H2_values = arange(0, 1 + precision_H2, precision_H2),
                                              kinship_eigenvalues_data = loadtxt(kinship_eigenvalues_filename),
                                              samples = progress_bar.ProgressBarIter(samples),
                                              distributions_filename = save_distributions_filename)
        else:
            estimates = _get_estimates(estimate_grid, estimates_filename)
            return build_heritability_cis_from_kinship(h2_values = arange(0, 1 + precision_h2, precision_h2),
                                                       H2_values = arange(0, 1 + precision_H2, precision_H2),
                                                       kinship_eigenvalues_data = loadtxt(kinship_eigenvalues_filename),
                                                       estimates = estimates, 
                                                       confidence = confidence,
                                                       samples = progress_bar.ProgressBarIter(samples),
                                                       distributions_filename = save_distributions_filename,
                                                       output_filename = output_filename)


if __name__ == '__main__':
    # Parse arguments
    parser = AlbiArgumentParser(prog=os.path.basename(sys.argv[0])) #,  usage=albi_USAGE)
        
    group_load = parser.add_mutually_exclusive_group(required=True)
    group_load.add_argument('-l', '--load_dist_filename',  type=str, help="Filename from which to load the estimated distributions.")
    group_load.add_argument('-k', '--kinship_eigenvalues', type=str, help="A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format.") 
    parser.add_argument('-p', '--precision', type=int, default=100, help="The number of grid points of the true heritability values, for which the estimator distributions are estimated. Effectively, this is the precision at which the CIs will be given (e.g., 100 grid points = 0.01 precision).")
    parser.add_argument('-d', '--distribution_precision', type=int, default=100, help="The number of grid points at which each estimator distribution is estimated. This controls the accuracy of estimation.") 
    parser.add_argument('-n', '--samples', type=int, default=1000, help="Number of random bootstrap samples to use for estimation.")
    parser.add_argument('-s', '--save_dist_filename', type=str, help="Filename at which to save the estimated distributions.")
    
    group_estimates = parser.add_mutually_exclusive_group(required=False)
    group_estimates.add_argument('-f', '--estimates_filename', type= str, help="A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one.")
    group_estimates.add_argument('-g', '--estimate_grid', type=int, help="A grid of heritability estimates. A CI will be calculated for each one (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1).")
    parser.add_argument('-c', '--confidence', type=float, default=0.95, help = "The required confidence level for the CIs.")
    parser.add_argument('-o', '--output_filename', type=str, help="File at which to write the calculated CIs.")
    
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
             precision_h2 = 1.0 / args.precision,
             precision_H2 = 1.0 / args.distribution_precision,
             samples = args.samples,
             confidence = args.confidence,
             output_filename = args.output_filename)



