#!/usr/bin/python

import os
import sys
import argparse
from numpy import arange, loadtxt, savetxt, hstack, vstack, newaxis, concatenate, array, ones, mean
import fiesta_lib
import progress_bar

class FiestaArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

FIESTA_USAGE = """
See https://github.com/cozygene/albi for full documentation about usage.

   python fiesta.py --kinship_eigenvalues filename 
                 [--use_eigenvectors_as_covariates <list of #>]
                 [--kinship_eigenvectors filename]
                 [--covariates filename]
                 [--add_intercept True/False]                 
                 (--estimates_filename filename] 
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                 [--iterations <# of iterations>] 
                 [--output_filename filename]
"""

if __name__ == '__main__':
  #
  # Parse arguments
  #
  parser = FiestaArgumentParser(prog=os.path.basename(sys.argv[0]), usage=FIESTA_USAGE)
      
  parser.add_argument('-k', '--kinship_eigenvalues', type=str, help="A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format.") 
  parser.add_argument('-u', '--use_eigenvectors_as_covariates', type=str, default='-1', help="A comma-separated list detailing which eigenvectors should be used as covariates.")
  parser.add_argument('-v', '--kinship_eigenvectors', type=str, help="A file containing the eigenvectors of the kinship matrix, one eigenvector per column, in text format.")
  parser.add_argument('-x', '--covariates', type=str, help="A file containing the covariates, one covariate per column, in text format.")
  parser.add_argument('-i', '--add_intercept', type=eval, default=True, help="If using covariates, add an intercept covariate (or only use an intercept covariate, if not covariates file supplied.")
  parser.add_argument('-n', '--iterations', type=int, default=1000, help="Number of iterations to use for estimation.")
  
  group_estimates = parser.add_mutually_exclusive_group(required=False)
  group_estimates.add_argument('-f', '--estimates_filename', type= str, help="A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one.")
  group_estimates.add_argument('-g', '--estimate_grid', type=int, default=10, help="A grid of heritability estimates. A CI will be calculated for each one (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1).")
  parser.add_argument('-c', '--confidence', type=float, default=0.95, help = "The required confidence level for the CIs.")
  parser.add_argument('-o', '--output_filename', type=str, help="Filename to which to write the calculated CIs.")
  
  args = parser.parse_args()

  #
  # Validate arguments
  #
  for filename in [args.kinship_eigenvalues,
                   args.kinship_eigenvectors,
                   args.covariates,
                   args.estimates_filename]:
    if filename and not os.path.exists(filename):
      print("File %s does not exist." % filename); sys.exit(2)

  try:
    use_eigenvectors_as_covariates = map(int, args.use_eigenvectors_as_covariates.split(','))
  except:
    print("Cannot parse --use_eigenvectors_as_covariates flag. It should be a comma-separated list of integers."); sys.exit(2)

  if args.iterations <= 0:
    print("Number of iterations should be a positive integer."); sys.exit(2)
  
  if args.estimate_grid <= 0:
    print("Grid size should be a positive integer."); sys.exit(2)
  
  if not (0 <= args.confidence <= 1):
    print("Confidence is a number between 0 and 1."); sys.exit(2)
  
  if args.kinship_eigenvalues is None:
    print("Kinship matrix eigenvalues file is required."); sys.exit(2)

  try:
    kinship_eigenvalues = loadtxt(args.kinship_eigenvalues)
  except:
    print("Failed reading eigenvalues file."); raise

  # 
  # Calculate CIs if needed
  #
  if args.estimates_filename:
    try:
      estimates = loadtxt(args.estimates_filename)
    except:
      print("Failed reading estimates file."); raise
        
  else:
    estimates = arange(0, 1 + 1.0/args.estimate_grid, 1.0/args.estimate_grid)


  # Decide if it's the simpler case or the general case
  if args.kinship_eigenvectors and not args.covariates:
      print("If using eigenvectors, covariates file must be supplied."); sys.exit(2)
  if args.covariates and not args.kinship_eigenvectors:
      print("If using covariates, eigenvectors file must be supplied."); sys.exit(2)
  if (args.add_intercept or args.covariates) and args.kinship_eigenvectors:
    # General case
    try:
      kinship_eigenvectors = loadtxt(args.kinship_eigenvectors)
    except:
      print("Failed reading eigenvectors file."); raise
    
    if args.covariates is not None: 
        try:
          covariates = loadtxt(args.covariates)
        except:
          print("Failed reading covariates file."); raise
        if args.add_intercept:
            if not any(mean(covariates == 1, axis=0) == 1):
                covariates = hstack([ones((len(kinship_eigenvalues), 1)), covariates])
    else:
        covariates = ones((len(args.kinship_eigenvalues), 1))
  
    print("Building heritability CIs...")   
    cis = fiesta_lib.calculate_cis_general(h_hat_values = estimates, 
                                           kinship_eigenvalues = kinship_eigenvalues,                                          
                                           kinship_eigenvectors = kinship_eigenvectors,
                                           covariates = covariates,
                                           iterations = args.iterations,
                                           alpha=1-args.confidence, 
                                           tau=0.4, 
                                           use_convergence_criterion=False)

  if not args.covariates and not args.kinship_eigenvectors:
    # The simpler case
    
    print("Building heritability CIs...")
    cis = fiesta_lib.calculate_cis_eigenvectors(h_hat_values = estimates, 
                                                kinship_eigenvalues = kinship_eigenvalues, 
                                                eigenvectors_as_X = use_eigenvectors_as_covariates,                                                 
                                                iterations = args.iterations,
                                                alpha=1-args.confidence, 
                                                tau=0.4, 
                                                use_convergence_criterion=False)
    
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




