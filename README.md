## ALBI  (*A*ccurate *L*MM-*B*ased Confidence *I*ntervals)

ALBI is a method for the estimation of the distribution of the heritability estimator, and for the construction of accurate confidence intervals (CIs). ALBI can be used as an add-on to existing methods for heritability and variance components estimation. ALBI is described in the following [paper](http://).

ALBI's input is the eigenvalues of a kinship matrix, and it produces accurate confidence intervals for a set of heritability estimates.

ALBI can be used as a command-line program. It can be also used as a Python module. Code written by Reut Yedidim and Regev Schweiger.

### Simple Example

The following reads the eigenvalues from a file, calculates 95% CIs over a grid of heritability estimates, and outputs it to a file:

```
   albi.py --kinship_eigenvalues ./data/eigenvalues.txt --estimate_grid 100 --output_filename cis.txt
```
For more information and options, see below.

## Installation

There are several ways to install and use ALBI.

### 1. Download ZIP

If you only want to use ALBI as a command line, you can simply download the code and run it as a standard program:

1. Press the 'Download ZIP' button on the right side of the screen
2. Extract the ZIP file to folder
3. Run ALBI

```
   albi.py [flags]
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
    sudo pip uninstall albi
```    

## Running
=======

### Command Line

Run from anywhere:
```
albi.py ---kinship_eigenvalues [file] --estimate_grid [number]
```
 
### Importing albi to a Python Script

```
    >>> import albi_lib
    >>> distributions = albi_lib.calculate_probability_intervals(...)
    >>> cis = albi_lib.build_heritability_cis( distributions, ...)

```

