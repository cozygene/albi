import argparse
# import paper as ALBI


def _generate_parse_kinship_arguments( parser ):
    parser.add_argument( '-k', '--kinship_eigenvalues', dest = 'kinship_eigenvalues', required = True, help = "path to file containing the eigealues of the kinship matrix" ) 
    parser.add_argument( '-s', '--samples',             dest = 'samples',             nargs = '?',     type = int,     default = 1000, const = 1000, help = "Number of samples ALBI should take" ) #   monte_carlo_size
    parser.add_argument( '-p', '--precision',           dest = 'precision',           nargs = '?',     type = float,   default = 0.1,  const = 0.1,  help = "the intervals between the h^2 values ALBI checks" )
    
def _generate_build_CI_arguments( parser ):
    parser.add_argument( '-c', '--confidance',          dest = 'confidance',          required = True, help = "How exact should the resualt CI be (precentage between 0-100")
    group = parser.add_mutually_exclusive_group( required = True )
    group.add_argument( '-ef', '--estimates_filename', dest = 'estimates_filename',  help = "A filename for the heritability estimates for which we want to build CIs. A text file with one estimate per row." ) #estimate is h^2? to which func is this? isn't it the interval? do we need it if we have 'intervals'
    group.add_argument( '-eg', '--estimate_grid',      dest = 'estimate_grid',       nargs = '?',     type = int, const = 0.1, help = "Instead of giving a file with a list of estimates, we can ask ALBI to simply give us CIs over (0,a,2a,3a,...,1), for example with a=0.01. This parameter is the grid size (e.g., 0.01). This flag is mutually exclusive with estimates_filename." )
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='ALBI')
    subparsers = parser.add_subparsers(help='sub-command help')


    # OPTION 1
    parser_a = subparsers.add_parser( 'opt1', help = 'TODO: provide help' )
    _generate_parse_kinship_arguments(parser_a)
    _generate_build_CI_arguments( parser_a )

    # exclusive arguments
    parser_a.add_argument( '-d', '--save_distributions',  dest = 'save_distributions',  nargs = '?',     type = str,     default = None, const = None, help = "This is a filename to which the program will save the output of calculate_probability_intervals, which is a matrix (take a look at savetxt)." )
    parser_a.add_argument( '-o', '--output_filename',     dest = 'output_filename',     nargs = '?',     type = str,     default = None, const = None, help = "The filename to which we will output the results of ALBI. A text file of a matrix, where each row is (h^2 estimate, lower bound, upper bound" )

       
    # OPTION 2
    parser_b = subparsers.add_parser( 'opt2', help = 'TODO b help' )
    _generate_parse_kinship_arguments(parser_b)
    # exclusive arguments
    parser_b.add_argument( '-d', '--save_distributions',  dest = 'save_distributions',  required = True, help = "This is a filename to which the program will save the output of calculate_probability_intervals, which is a matrix (take a look at savetxt)." )


    # OPTION 3
    parser_c = subparsers.add_parser( 'opt3', help = 'TODO c help' )
    _generate_build_CI_arguments( parser_c )
    # exclusive arguments
    parser_c.add_argument( '-d', '--load_distributions',  dest = 'load_distributions',  required = True, help = "This is a filename to which the program will load the output of calculate_probability_intervals, which is a matrix (loadtxt). This flag is mutually exclusive with --input" )
    parser_c.add_argument( '-o', '--output_filename',     dest = 'output_filename',     nargs = '?',     type = str,     default = None, const = None, help = "The filename to which we will output the results of ALBI. A text file of a matrix, where each row is (h^2 estimate, lower bound, upper bound" )


    a = parser.parse_args()

    print 'kinship_eigenvalues', a.kinship_eigenvalues
    print 'samples', a.samples
    print 'precision', a.precision
    print 'confidance', a.confidance
    print 'estimate_grid', a.estimate_grid
    print 'estimates_filename', a.estimates_filename
    print 'save_distributions', a.save_distributions
    print 'load_distributions', a.load_distributions
    print 'output_filename', a.output_filename

    print dir(a)


# # OPTION 1 - full
# arguments.kinship
# arguments.samples
# arguments.precision
# arguments.confidance
# arguments.estimates_* # or

# arguments.save_distributions # optional
# arguments.save_output #??? yes?


# # OPTION 2 - create file
# arguments.kinship
# arguments.samples
# arguments.precision

# arguments.save_distributions

# # OPTION 3-  extract from filr
# arguments.confidance
# arguments.estimates_* # or

# arguments.load_distributions 
# arguments.save_output #??? yes?





"""
                                    precision
def calculate_probability_intervals(true_h2, interval_boundaries, kinship_eigenvalues, 
                                    eigenvectors_as_X=[-1], REML=True, monte_carlo_size=1000, n_chunks=1, seed=0):   (returns prob)


                                                    estimates_*
def build_heritability_cis(accept_regions, true_h2s, values):       (returns cis)


def build_acceptance_region_given_distribution(true_h2s, distribution, real_value, alpha, verbose=False):  (returns low, high, coverage)


"""



