#!/usr/bin/python

import os
import sys
import argparse
from numpy import arange, loadtxt, savetxt
import albi as ALBI 

class MyArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        print( "To see full help: %s -h/--help" % self.prog )
        sys.exit(2)

SHELL_WIDTH = 70

ALBI_USAGE  = \
"""
If you want to build_heritability_cis from a kinship_eigenvalues please specify        %(prog)s --kinship_eigenvalues  [kinship file]        --estimates_filename/--estimate_grid
If you want to build_heritability_cis from a distribution file please specify          %(prog)s --load_distributions   [distribuation file]  --estimates_filename/--estimate_grid
If you just want to calculate_probability_intervals please specify                     %(prog)s  --kinship_eigenvalues [kinship file]        --save_distributions [distribuation file]

There are optional arguments you can specify as you can see below:
"""


def ioctl_GWINSZ(fd):
    try:
        import fcntl, termios, struct
        _, width = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,'1234'))
    except:
        return
    return width

def get_shell_width():
    width = ioctl_GWINSZ(0)

    if not width:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            width = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass

    if not width:
        width = os.environ.get( 'COLUMNS' )

    return width if width else SHELL_WIDTH


class ProgressBarIter( object ):
    def __init__( self, length, stdout = sys.stdout, fill = '#', width = get_shell_width() ):
        self.length = float(length)
        if stdout.isatty():
          self.stdout = stdout
        else:
          self.stdout = sys.stdout
        self.current = 0
        self.fill = fill
        self.width = width
        self.prefix = '| '
        self.suffix = ' |'
        self.precentage = " {precentage}%"

    def __iter__( self ):
        return self

    def __next__( self ): # Python 3 uses __next__ for iterator
        self.next()

    def next( self ): # Python 2 uses next for iterator
        self.current += 1
        precentage = (self.current / self.length )
        precentage_str = self.precentage.format(precentage = precentage * 100)
        width_left = self.width - len(self.prefix) - len(self.suffix) - len(precentage_str)
        process_bar_fill_str = ( '#' * int(precentage * width_left) ).ljust( width_left )
        self._write( self.prefix + process_bar_fill_str + precentage_str + self.suffix )

        if precentage == 1:
            raise StopIteration

    def _write( self, line ) :
        self.stdout.write('\r')
        self.stdout.write( line )
        self.stdout.flush()

    def __len__( self ):
      return self.length

# opt2
def calculate_probability_intervals( precision_h2, precision_H2, kinship_eigenvalues_data, samples = 1000, distributions_filename = None):
    print( "Calculating probability intervals..." )
    distributions =  ALBI.calculate_probability_intervals( h2_values = arange(0, 1 + precision_h2, precision_h2), 
                                                           H2_values = arange(0, 1 + precision_H2, precision_H2), 
                                                           kinship_eigenvalues = kinship_eigenvalues_data, 
                                                           n_random_samples = samples
                                                         )
    print( "Done calculating" )
    if distributions_filename:
        savetxt( distributions_filename, distributions )
        print( "Distributions are written to '%s' " % distributions_filename )
    return distributions

# opt3 
def build_heritability_cis_from_distributions( precision_h2, precision_H2, all_distributions_data, estimates, confidence = 0.95, output_filename = None ):
    print( "Building heritability CIs..." )
    cis =  ALBI.build_heritability_cis( h2_values = arange(0, 1 + precision_h2, precision_h2), 
                                        H2_values = arange(0, 1 + precision_H2, precision_H2), 
                                        all_distributions = all_distributions_data, 
                                        estimated_values = estimates, 
                                        confidence = confidence, 
                                        use_randomized_cis = False
                                      )
    print( "Done building CIs" )
    if output_filename:
        savetxt( output_filename, cis )
        print( "Heritability CIs are written to '%s' " % output_filename )
    return cis

# opt1
def build_heritability_cis_from_kinship( precision_h2, precision_H2, kinship_eigenvalues_data, estimates, confidence = 0.95, samples = 1000, distributions_filename = None, output_filename = None ):
    distribution = calculate_probability_intervals(  precision_h2 = precision_h2,
                                                     precision_H2 = precision_H2,
                                                     kinship_eigenvalues_data = kinship_eigenvalues_data,
                                                     samples = samples,
                                                     distributions_filename = distributions_filename
                                                    )

    return build_heritability_cis_from_distributions( precision_h2 = precision_h2,
                                                      precision_H2 = precision_H2, 
                                                      all_distributions_data = distribution, 
                                                      estimates = estimates, 
                                                      confidence = confidence, 
                                                      output_filename = output_filename
                                                      )
    print( "Done" )


def _get_estimates( estimate_grid, estimates_filename ):
    if estimate_grid:
        estimates = arange(0, 1, 1.0/estimate_grid) #is that right? TODO
    elif not os.path.exists( estimates_filename ) :
        print( "The file '%s' doesn't exist. Exiting" % estimates_filename )
        sys.exit(2)  
    else:
        estimates = file( args.estimates_filename, 'rb' ).read()
    return estimates


def run_albi( kinship_eigenvalues_filename = None,
              estimate_grid = None,
              estimates_filename = None,
              save_distributions_filename = None,
              load_distributions_filename = None,
              precision_h2 = 0.1,
              precision_H2 = 0.1,
              samples = 1000,
              confidence = 0.95,
              output_filename = None,
              print_comments = False ):
    """
    There are three options to run ALBI:
    opt1 - build_heritability_cis_from_kinship. The full run of ALBI. gets a kinship eigenvalues and returns the heritability CIs. (This is actually opt2 + opt3)
    opt2 - calculate_probability_intervals. gets a kinship eigenvalues and calculate probability intervals (distributions) and saves it for later.
    opt3 - build_heritability_cis_from_distributions. gets the distributions file from opt2 and returns the heritability CIs.

    This function parses the script arguments and runs the right option according to the specified arguments:
    for opt1 - you have to specify kinship_eigenvalues + estimates
    for opt2 - you have to specify kinship_eigenvalues + save_distributions
    for opt3 - you have to specify load_distributions + estimates

    arguments:
    all except print_comments are ALBI arguments (you can see description in __main__ below)
    @ print_comments: if True prints to stdout the comments, otherwise doesn't print
    """
    
    if not print_comments:
        sys.stdout = open(os.devnull, "w") #set stdout to null

    # if a kinship_eigenvalues wasn't speified - the user wants to run opt3 or he was wrong. If it was specified the user wants to run opt1 or opt2
    if kinship_eigenvalues_filename is None: # opt3 or error
        if load_distributions_filename is None:
            print( "If you want to build_heritability_cis from a kinship_eigenvalues please specify        --kinship_eigenvalues  [kinship file]        --estimates_filename/--estimate_grid\n" \
                + "If you want to build_heritability_cis from a distribution file please specify          --load_distributions   [distribuation file]  --estimates_filename/--estimate_grid\n" \
                + "If you just want to calculate_probability_intervals please specify                     --kinship_eigenvalues  [kinship file]        --save_distributions [distribuation file]" 
                )
            return None
        elif not os.path.exists( load_distributions_filename ):
            print(  "The file '%s' doesn't exist. Exiting" % load_distributions_filename ) #use logging ifo?TODO
            return None
        elif estimates_filename is None and estimate_grid is None:
            print( "Few arguments are missing.\n" \
                  + "If you want to build_heritability_cis from a distribution file please also specify --estimate_grid/--estimates_filename"
                  )
            return None
        else:
            # run opt3
            estimates = _get_estimates( estimate_grid, estimates_filename )
            return build_heritability_cis_from_distributions( precision_h2 = precision_h2,
                                                              precision_H2 = precision_H2,
                                                              all_distributions_data = loadtxt( load_distributions_filename ), 
                                                              estimates = estimates,
                                                              confidence = confidence,
                                                              output_filename = output_filename
                                                            )

    elif not os.path.exists( kinship_eigenvalues_filename ) :
                print( "The file '%s' doesn't exist. Exiting" % kinship_eigenvalues_filename )
                return None

    else: # opt1 or opt2
        if estimates_filename is None and estimate_grid is None: #opt2 or error
            if save_distributions_filename is None:
                print( "Few arguments are missing.\n" \
                    + "If you want to build_heritability_cis please also specify --estimate_grid/--estimates_filename.\n" \
                    + "If you want to calculate_probability_intervals please also specify --save_distributions file"
                    )
                return None
            else:
                # run opt2
                return calculate_probability_intervals( precision_h2 = precision_h2,
                                                        precision_H2 = precision_H2,
                                                        kinship_eigenvalues_data = loadtxt( kinship_eigenvalues_filename ),
                                                        samples = ProgressBarIter( samples ),
                                                        distributions_filename = save_distributions_filename
                                                       )
        else:
            # run opt1
            estimates = _get_estimates( estimate_grid, estimates_filename )
            return build_heritability_cis_from_kinship(  precision_h2 = precision_h2,
                                                         precision_H2 = precision_H2,
                                                         kinship_eigenvalues_data = loadtxt( kinship_eigenvalues_filename ),
                                                         estimates = estimates, 
                                                         confidence = confidence,
                                                         samples = ProgressBarIter( samples ),
                                                         distributions_filename = save_distributions_filename,
                                                         output_filename = output_filename
                                                       )
    if not print_comments:
        sys.stdout = sys.__stdout__ # return stdout

if __name__ == '__main__':
    parser = MyArgumentParser(prog=os.path.basename(sys.argv[0]),  usage=ALBI_USAGE)

    parser.add_argument( '--kinship_eigenvalues',                                                                  help = "path to file containing the eigealues of the kinship matrix" ) 
    group = parser.add_mutually_exclusive_group( required = False )
    group.add_argument(  '--estimates_filename',                     type = str,                                   help = "A filename for the heritability estimates for which we want to build CIs. A text file with one estimate per row." )
    group.add_argument(  '--estimate_grid',                          type = int,                                 help = "How many 'jumps' to ahve between 0 and 1" )
    parser.add_argument( '--save_distributions',                     type = str,                                   help = "default is None. This is a filename to which the program will save the output of calculate_probability_intervals, which is a matrix (take a look at savetxt)." )
    parser.add_argument( '--load_distributions',                     type = str,                                   help = "default is None. This is a filename to which the program will load the output of calculate_probability_intervals, which is a matrix (loadtxt). This flag is mutually exclusive with --input" )
    
    parser.add_argument( '--precision',             nargs = '?',     type = float,   default = 0.1,  const = 0.1,  help = "default value is 0.1.  the intervals between the h^2 values ALBI checks" )
    parser.add_argument( '--precision2',            nargs = '?',     type = float,   default = 0.1,  const = 0.1,  help = "default value is 0.1.  I dont know what this is for :p" ) 
    parser.add_argument( '--samples',               nargs = '?',     type = int,     default = 1000, const = 1000, help = "default value is 1000. Number of samples ALBI should take" ) #   monte_carlo_size
    parser.add_argument( '--confidence',            nargs = '?',     type = float,   default = 0.95, const = 0.95, help = "default value is 0.95. How exact should the resualt CI be (precentage between 0-100")
    
    parser.add_argument( '--output_filename',                        type = str,                                   help = "default is None. The filename to which we will output the results of ALBI. A text file of a matrix, where each row is (h^2 estimate, lower bound, upper bound" )

    args = parser.parse_args()

    print( "Hello from ALBI" )
    print( "To see full help run: %s -h/--help\n" % parser.prog )
    
    if args.kinship_eigenvalues is None and args.estimates_filename is None and args.estimate_grid is None:
        parser.print_help()
    else:
        run_albi( kinship_eigenvalues_filename = args.kinship_eigenvalues,
                  estimate_grid = args.estimate_grid, 
                  estimates_filename = args.estimates_filename,
                  save_distributions_filename = args.save_distributions,
                  load_distributions_filename = args.load_distributions,
                  precision_h2 = args.precision,
                  precision_H2 = args.precision2,
                  samples = args.samples,
                  confidence = args.confidence,
                  output_filename = args.output_filename,
                  print_comments = True )



