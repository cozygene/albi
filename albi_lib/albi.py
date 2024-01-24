#!/usr/bin/python

import os
import sys
import argparse
from numpy import arange, loadtxt, savetxt, hstack, vstack, newaxis, concatenate, array

import progress_bar
import albi_lib

class AlbiArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

ALBI_USAGE = """
See https://github.com/cozygene/albi for full documentation about usage.

   python albi.py --kinship_eigenvalues filename 
                 [--use_eigenvectors_as_covariates <list of #>]
                 [--precision <# of grid points>] 
                 [--distribution_precision <# of grid points>] 
                 [--samples <# of random samples>] 
                 (--estimates_filename filename
                 (--estimates_filename filename] 
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                 [--output_filename filename]
or

   python albi.py --kinship_eigenvalues filename
                  --kinship_eigenvectors filename
                  --covariates filename
                 [--precision <# of grid points>] 
                 [--distribution_precision <# of grid points>] 
                 [--samples <# of random samples>] 
                 (--estimates_filename filename 
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                 [--output_filename filename]
"""

if __name__ == '__main__':
  #
  # Parse arguments
  #
  parser = AlbiArgumentParser(prog=os.path.basename(sys.argv[0]), usage=ALBI_USAGE)
      
  group_load = parser.add_mutually_exclusive_group(required=True)
  group_load.add_argument('-l', '--load_dist_filename',  type=str, help="Filename from which to load the estimated distributions.")
  group_load.add_argument('-k', '--kinship_eigenvalues', type=str, help="A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format.") 
  parser.add_argument('-u', '--use_eigenvectors_as_covariates', type=str, default='-1', help="A comma-separated list detailing which eigenvectors should be used as covariates.")
  parser.add_argument('-v', '--kinship_eigenvectors', type=str, help="A file containing the eigenvectors of the kinship matrix, one eigenvector per column, in text format.")
  parser.add_argument('-x', '--covariates', type=str, help="A file containing the covariates, one covariate per column, in text format.")
  parser.add_argument('-p', '--precision', type=int, default=100, help="The number of grid points of the true heritability values, for which the estimator distributions are estimated. Effectively, this is the precision at which the CIs will be given (e.g., 100 grid points = 0.01 precision).")
  parser.add_argument('-d', '--distribution_precision', type=int, default=100, help="The number of grid points at which each estimator distribution is estimated. This controls the accuracy of estimation.") 
  parser.add_argument('-n', '--samples', type=int, default=1000, help="Number of random bootstrap samples to use for estimation.")
  parser.add_argument('-s', '--save_dist_filename', type=str, help="Filename to which to save the estimated distributions.")
  
  group_estimates = parser.add_mutually_exclusive_group(required=False)
  group_estimates.add_argument('-f', '--estimates_filename', type= str, help="A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one.")
  group_estimates.add_argument('-g', '--estimate_grid', type=int, default=10, help="A grid of heritability estimates. A CI will be calculated for each one (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1).")
  parser.add_argument('-c', '--confidence', type=float, default=0.95, help = "The required confidence level for the CIs.")
  parser.add_argument('-o', '--output_filename', type=str, help="Filename to which to write the calculated CIs.")
  
  args = parser.parse_args()

  #
  # Validate arguments
  #
  for filename in [args.load_dist_filename, 
                   args.kinship_eigenvalues,
                   args.kinship_eigenvectors,
                   args.covariates,
                   args.estimates_filename]:
    if filename and not os.path.exists(filename):
      print("File %s does not exist." % filename); sys.exit(2)

  try:
    use_eigenvectors_as_covariates = list(map(int, args.use_eigenvectors_as_covariates.split(',')))
  except:
    print("Cannot parse --use_eigenvectors_as_covariates flag. It should be a comma-separated list of integers."); sys.exit(2)

  if args.precision <= 0 or args.distribution_precision <= 0:
    print("Precision should be a positive integer."); sys.exit(2)

  if args.samples <= 0:
    print("Number of samples should be a positive integer."); sys.exit(2)
  
  if args.estimate_grid <= 0:
    print("Grid size should be a positive integer."); sys.exit(2)
  
  if not (0 <= args.confidence <= 1):
    print("Confidence is a number between 0 and 1."); sys.exit(2)
  
  #
  # First, estimate the distribution or load it
  #
  if args.kinship_eigenvalues is None and args.load_dist_filename is None:
    print("A pre-estimated distribution file or kinship matrix file(s) are required."); sys.exit(2)

  if args.load_dist_filename:
    try:
      # Read distributions file
      all_distributions_data = loadtxt(args.load_dist_filename, delimiter='\t', skiprows=1)
      header = file(args.load_dist_filename, 'rb').readline().strip().split()[1:]
      
      h2_values = all_distributions_data[:,0]
      H2_values = array([float(header[0])] + [float(field_name.split('-')[1]) for field_name in header[1:-2]] + [float(header[-1])])
      distributions = all_distributions_data[:,1:]
    except:
      print("Failed reading distribution file."); raise

  else:
    precision_h2 = 1.0 / args.precision
    precision_H2 = 1.0 / args.distribution_precision 
    h2_values = arange(0, 1 + precision_h2, precision_h2)
    H2_values = arange(0, 1 + precision_H2, precision_H2)

    # We need to estimate the distribution
    try:
      kinship_eigenvalues = loadtxt(args.kinship_eigenvalues)
    except:
      print("Failed reading eigenvalues file."); raise

    # Decide if it's the simpler case or the general case
    if args.kinship_eigenvectors and not args.covariates:
        print("If using eigenvectors, covariates file must be supplied."); sys.exit(2)
    if args.covariates and not args.kinship_eigenvectors:
        print("If using covariates, eigenvectors file must be supplied."); sys.exit(2)
    if args.covariates and args.kinship_eigenvectors:
      # General case
      try:
        kinship_eigenvectors = loadtxt(args.kinship_eigenvectors)
      except:
        print("Failed reading eigenvectors file."); raise
      
      try:
        covariates = loadtxt(args.covariates)
      except:
        print("Failed reading covariates file."); raise
      
      print("Estimating distributions...")        
      distributions = albi_lib.estimate_distributions_general(h2_values = h2_values, 
                                                              H2_values = H2_values, 
                                                              kinship_eigenvalues = kinship_eigenvalues,
                                                              kinship_eigenvectors = kinship_eigenvectors,
                                                              covariates = covariates,
                                                              n_random_samples = progress_bar.ProgressBarIter(args.samples))

    if not args.covariates and not args.kinship_eigenvectors:
      # The simpler case
      print("Estimating distributions...")
      distributions = albi_lib.estimate_distributions_eigenvectors(h2_values = h2_values, 
                                                                   H2_values = H2_values, 
                                                                   kinship_eigenvalues = kinship_eigenvalues, 
                                                                   n_random_samples = progress_bar.ProgressBarIter(args.samples),
                                                                   eigenvectors_as_X = use_eigenvectors_as_covariates)

  #
  # If we want to do is save it, then just save it
  #
  if args.save_dist_filename:
      header = ['ALBI'] + ['0'] + ["%f-%f" % (H2_values[i], H2_values[i+1]) for i in range(len(H2_values)-1)] + ['1']
      savetxt(fname = args.save_dist_filename, 
              X = hstack([h2_values[:, newaxis], distributions]),
              fmt = "%1.5f",
              delimiter = '\t',
              header = '\t'.join(header),
              comments = '')
      print("Distributions are written to '%s'." % args.save_dist_filename)
      sys.exit(0)

  # 
  # Calculate CIs if needed
  #
  if args.estimates_filename:
    try:
      estimates = loadtxt(args.estimates_filename, ndmin=1)
    except:
      print("Failed reading estimates file."); raise
        
  else:
    estimates = arange(0, 1 + 1.0/args.estimate_grid, 1.0/args.estimate_grid)
  
  print("Building heritability CIs...")
  cis = albi_lib.build_heritability_cis(h2_values = h2_values, 
                                        H2_values = H2_values,
                                        all_distributions = distributions, 
                                        estimated_values = estimates, 
                                        confidence = args.confidence, 
                                        use_randomized_cis = False)
  header = ['Estimate', 'CI_lower_bound', 'CI_upper_bound']
  if args.output_filename:
    savetxt(fname = args.output_filename, 
            X = hstack([estimates[:, newaxis], cis]),
            fmt = "%1.5f",
            delimiter = '\t',
            header = '\t'.join(header),
            comments = '')
    print("Heritability CIs written to '%s'." % args.output_filename)
  else:
    # Output to screen
    print('\t'.join(header))
    for e, ci in zip(estimates, cis):
      print('\t'.join(["%1.5f" % x for x in [e] + list(ci)]))

  sys.exit(0)




