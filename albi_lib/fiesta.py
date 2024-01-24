#!/usr/bin/python

import os
import sys
import argparse
from numpy import *
import numpy.linalg

import progress_bar
import fiesta_lib

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
                 [--no_intercept]                 
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
  parser.add_argument('-u', '--use_eigenvectors_as_covariates', type=str, help="A comma-separated list detailing which eigenvectors should be used as covariates.")
  parser.add_argument('-v', '--kinship_eigenvectors', type=str, help="A file containing the eigenvectors of the kinship matrix, one eigenvector per column, in text format.")
  parser.add_argument('-x', '--covariates', type=str, help="A file containing the covariates, one covariate per column, in text format.")
  parser.add_argument('-i', '--no_intercept', action='store_true', help="If using covariates, don't add an intercept covariate.")
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

  if args.iterations <= 0:
    print("Number of iterations should be a positive integer."); sys.exit(2)
  
  if args.estimate_grid <= 0:
    print("Grid size should be a positive integer."); sys.exit(2)
  
  if not (0 <= args.confidence <= 1):
    print("Confidence is a number between 0 and 1."); sys.exit(2)

  if args.confidence < 0.5:
    print("Warning: Confidence is small; usually values as 0.95 are used. Make sure this is on purpose.")
  
  if args.kinship_eigenvalues is None:
    print("Kinship matrix eigenvalues file is required."); sys.exit(2)

  try:
    kinship_eigenvalues = loadtxt(args.kinship_eigenvalues)

    if any(kinship_eigenvalues < 0):
        print("Warning: Some eigenvalues are negative. Rounding up."); 
  except:
    print("Failed reading eigenvalues file."); raise

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


  if args.covariates and not args.kinship_eigenvectors:
      print("If using covariates, eigenvectors file must be supplied."); sys.exit(2)

  # Decide if it's the simpler case or the general case
  if args.kinship_eigenvectors and (args.covariates is not None or not args.no_intercept):
    # General case
    try:
      kinship_eigenvectors = loadtxt(args.kinship_eigenvectors)
    except:
      print("Failed reading eigenvectors file."); raise

    if args.use_eigenvectors_as_covariates is not None:
        print("If eigenvectors were given explicitly, cannot use --use_eigenvectors_as_covariates - use covariates directly."); sys.exit(2)

    if args.covariates is not None: 
        try:
          covariates = loadtxt(args.covariates)
        except:
          print("Failed reading covariates file."); raise

        if not any(mean(covariates == 1, axis=0) == 1):
            covariates = hstack([ones((len(kinship_eigenvalues), 1)), covariates])
    else:
        print("Note: No covariates supplied, using a constant intercept covariate.")
        covariates = ones((len(kinship_eigenvalues), 1))
  
    # Check pathological kinships will work
    # We need to check that each of the eigenvectors corresponding to zero eigenvalues are spanned by covariates
    Xdagger = numpy.linalg.pinv(covariates)

    is_dangerous = where(isclose(kinship_eigenvalues, 0) | (kinship_eigenvalues < 0))[0]
    dangerous_eigenvectors = kinship_eigenvectors[:, is_dangerous]

    if not allclose(dot(covariates, dot(Xdagger, dangerous_eigenvectors)), dangerous_eigenvectors):
        print("*** Warning ***: Some eigenvalues are zero (and do not exactly correspond to eigenvectors spanned by covariates) - This might create unstable results.")

    # For numerical stability
    kinship_eigenvalues = maximum(kinship_eigenvalues, 1e-10)

    print("Building heritability CIs...")   
    cis = fiesta_lib.calculate_cis_general(h_hat_values = estimates, 
                                           kinship_eigenvalues = kinship_eigenvalues,                                          
                                           kinship_eigenvectors = kinship_eigenvectors,
                                           covariates = covariates,
                                           iterations = args.iterations,
                                           alpha=1-args.confidence, 
                                           tau=0.4, 
                                           use_convergence_criterion=False)

  else:
    # The simpler case

    # A subcase if not covariates were given at all in the general case - fall back
    if args.kinship_eigenvectors and args.covariates is None and args.no_intercept:
        eigenvectors_as_X = []

    elif args.use_eigenvectors_as_covariates is None:
        eigenvectors_as_X = [-1]   # Default

    else:
        try:
            if args.use_eigenvectors_as_covariates == '':
                eigenvectors_as_X = []
            else:
                eigenvectors_as_X = list(map(int, args.use_eigenvectors_as_covariates.split(',')))
        except:
            print("Cannot parse --use_eigenvectors_as_covariates flag. It should be a comma-separated list of integers."); sys.exit(2)

    eigenvectors_as_X = array(eigenvectors_as_X) % len(kinship_eigenvalues)

    # Check pathological kinships will work
    is_dangerous = isclose(kinship_eigenvalues, 0) | (kinship_eigenvalues < 0)
    if not set(where(is_dangerous)[0]) <= set(eigenvectors_as_X):
        print("*** Warning ***: Some eigenvalues are zero (and do not correspond to eigenvectors used as covariates) - This might create unstable results.")

    # For numerical stability
    kinship_eigenvalues = maximum(kinship_eigenvalues, 1e-10)

    print("Building heritability CIs...")
    cis = fiesta_lib.calculate_cis_eigenvectors(h_hat_values = estimates, 
                                                kinship_eigenvalues = kinship_eigenvalues, 
                                                eigenvectors_as_X = eigenvectors_as_X,                                                 
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




