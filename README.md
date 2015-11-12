## ALBI  (*A*ccurate *L*MM-*B*ased Confidence *I*ntervals)

ALBI is a method for the estimation of the distribution of the heritability estimator, and for the construction of accurate confidence intervals (CIs). ALBI can be used as an add-on to existing methods for heritability and variance components estimation. ALBI is described in the following [paper](http://).

ALBI's input is the eigenvalues of a kinship matrix, and it produces accurate confidence intervals for a set of heritability estimates.

ALBI can be used as a command-line program. It can be also used as a Python module. Code written by Reut Yedidim and Regev Schweiger.

### Simple Example

The following reads the eigenvalues from a file, calculates 95% CIs over a grid of heritability estimates, and outputs it to a file:

```
   python albi.py --kinship_eigenvalues ./data/eigenvalues.txt 
                  --estimate_grid 100 
                  --output_filename cis.txt
```
For more information and options, see below.

## Installation

There are several ways to install and use ALBI.

### 1. Download ZIP

If you only want to use ALBI as a command line, you can simply download the code and run it as a standard program:

1. Press the 'Download ZIP' button on the right side of the screen
2. Extract the ZIP file to folder
3. Run ALBI from the folder:
```
   python albi.py [flags]
```

**Note**: The [NumPy](http://www.numpy.org/) package is required for ALBI.

### 2. Install using `pip` or `setuptools`

In order to use ALBI and a Python package, it needs to be installed, using one of the following. If you have `pip` installed, run this on command line (from anywhere): 
```
   sudo pip install git+https://github.com/cozygene/albi
```

Alternatively, you can also clone the repository and do a manual install:
```
   git clone https://github.com/cozygene/albi
   sudo python setup.py install
```
### Uninstalling

To uninstall, run:
```
    sudo -H pip uninstall albi
```    
## Running ALBI

The following section describes all of ALBI's options and flags.

The method consists of two stages:

1. Estimating the distributions of the heritability estimator. This is the more computationally expensive stage, but it can be performed only once for each kinship matrix.
2. Using these distributions to build accurate CIs for a set of heritabilty estimates.
 
The two stages may be performed together, or they may be performed separately, with the distributions saved and loaded from a file. 

### 1. Estimating distributions and saving the results

The following command estimates the distributions, and saves the results to a file:
```
   python albi.py --kinship_eigenvalues filename 
                 [--precision <# of grid points>] 
                 [--distribution_precision <# of grid points>] 
                 [--samples <# of random samples>] 
                 --save_dist_filename filename
```

Details about each

prob0/1: It is sometimes useful... precision=1

### 2. Creating CIs from a file with pre-estimated distributions

The following command loads the distributions, builds CIs and saves them to a file:
```
   python albi.py --load_dist_filename filename 
                 (--estimates_filename filename
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                  --output_filename filename
```

Details

### 3. Performing both steps

Essentially the same, without saving/loading dists
```
   python albi.py --kinship_eigenvalues filename 
                 [--precision <# of grid points>] 
                 [--distribution_precision <# of grid points>] 
                 [--samples <# of random samples>] 
                 (--estimates_filename filename] 
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                  --output_filename filename
```

Details

### Full list of flags


Flag | Short | Value
------------ | -------------
`kinship_eigenvalues` | -k | Filename
`precision` | -p | fff

quiet


## Using ALBI as a Python library

```
    >>> import albi_lib
    >>> distributions = albi_lib.calculate_probability_intervals(...)
    >>> cis = albi_lib.build_heritability_cis(distributions, ...)

```
 

# ALBi // TODO !!
change aruments long name <br/>
change help and documentation
albi is everything you need


:hamster:
