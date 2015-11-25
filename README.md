## ALBI  (*A*ccurate *L*MM-*B*ased Confidence *I*ntervals)

ALBI is a method for the estimation of the distribution of the heritability estimator, and for the construction of accurate confidence intervals (CIs). ALBI can be used as an add-on to existing methods for heritability and variance components estimation. ALBI is described in the following [paper](http://biorxiv.org/content/early/2015/11/24/031492).

ALBI's input is the eigenvalues of a kinship matrix, and it produces accurate confidence intervals for a set of heritability estimates.

ALBI can be used as a command-line program. It can be also used as a Python module. Code written by Reut Yedidim and Regev Schweiger. For questions or comments, mail [schweiger@post.tau.ac.il](mailto:schweiger@post.tau.ac.il).

### Simple Example

The following reads the eigenvalues from a file, calculates 95% CIs over a grid of heritability estimates, and outputs it to a file:

```
   python albi.py --kinship_eigenvalues ./data/eigenvalues.txt  \
                  --estimate_grid 100                           \
                  --output_filename cis.txt
```

The output is the CIs for heritability estimates of 0, 0.01, ..., 0.99, 1:

```
Estimate   CI_lower_bound  CI_upper_bound
0.00000    0.00000         0.65000
0.01000    0.00000         0.65430
0.02000    0.00000         0.65860
0.03000    0.00000         0.73156

...

0.98000    0.39927         1.00000
0.99000    0.40494         1.00000
1.00000    0.41000         1.00000
```

For example, the CI for an estimate of 0.02 is [0, 0.6586].

For more information and options, see below.

## Installation

There are several ways to install and use ALBI.

### 1. Download ZIP

If you only want to use ALBI as a command line tool, you can simply download the code and run it as a standard program:

1. Press the 'Download ZIP' button on the right side of the screen
2. Extract the ZIP file to a folder
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

ALBI consists of two stages:

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

The flags are as follows:

* `kinship_eigenvalues` - A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format. This could be created, for example, with GCTA's `--pca` flag.
* `precision` - The number of grid points of the true heritability values, for which the estimator distributions are estimated. Effectively, this is the precision at which the CIs will be given (e.g., 100 grid points = 0.01 precision). Default is 100.
* `distribution_precision` - The number of grid points at which each estimator distribution is estimated. This controls the accuracy of estimation. Default is 100.
* `samples` - Number of random bootstrap samples to use for estimation. Default is 1000.
* `save_dist_filename` - Filename at which to save the estimated distributions.

One scenario at which the distributions may be useful in themselves, is using the boundary probabilities (0 and 1) for a preliminary assessment of CIs. This follows from the reasoning that narrow CIs translate to small boundary probabilities (see paper for more details). To calculate only the boundary probabilities, use `--distribution_precision 1`.

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

The flags are as follows:

* `load_dist_filename` - Filename from which to load the estimated distributions.
* `estimates_filename` - A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one. 
* `estimate_grid` - Alternatively, one can ask ALBI to calculate CIs for a grid of heritability estimates (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1).
* `confidence` - The required confidence level for the CIs. Default is 0.95 (95% CIs).
* `output_filename` - File at which to write the calculated CIs.

### 3. Performing both steps

The two steps may be performed consecutively, without saving or loading the estimated distributions from a file:
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


## Using ALBI as a Python library

ALBI may be used as a Python library. There are two main functions, corresponding to the two stages described above. The Python code is self-documenting:

```Python
   >>> import albi_lib
   >>> help(albi_lib.estimate_distributions)
```
```
   estimate_distributions(h2_values, H2_values, kinship_eigenvalues, n_random_samples=100, eigenvectors_as_X=[-1], REML=True, seed=0)
   
    Across a grid of possible estimated values H^2, approximately calculate the probability of either evaluating a boundary 
    estimate (for the boundaries of the grid) or the the estimate falling between each grid points. The probability is 
    defined for a true value of heritability h^2. The probability is estimated with a parametric bootstrap.

    Arguments:
        h2_values - a vector of size N, of all possible values of h^2
        H2_values - a vector of size M, of a grid of possible values of H^2
        kinship_eigenvalues - A vector of size K of the eigenvalues of the kinship matrix, in decreasing order.
        n_random_samples - The number of random samples to use for the parameteric bootstrap. Can be an int or an iterable with the __len__ func implemented
        eigenvectors_as_X - A list of indices, of which eigenvectors of the kinship matrix are fixed effects.
        REML - True is REML, False if ML.
        seed - A seed for the random generator used for the random samples.

    Returns:
        A matrix of size N x (M + 1) of the probabilities, where:
        - The cell at index (i, 0) is the probability of the estimate being smaller than the smallest grid point;
        - The cell at index (i, M) is the probability of the estimate being larger than the largest grid point;
        - The cell at index (i, j), for 0<j<M, is the probability of estimate being between the (j-1)-th and j-th grid points;

          all of the above are for the i-th h2 value.
```
```Python
   >>> help(albi_lib.build_heritability_cis)
```
```
   build_heritability_cis(h2_values, H2_values, all_distributions, estimated_values, confidence=0.95, use_randomized_cis=False, seed=0)

   Build confidence intervals for a set of estimated values, given estimator distributions.

    Arguments:
        h2_values - a vector of size N, of values of true h^2 values at which the estimator
                    distribution is given.
        H2_values - a vector of size M, of a grid of values of estimated values at which 
                    the estimator distribution is given, for each true value h^2.
        all_distributions - a matrix of size N x (M + 1), where the i-th row is the estimator
                            distribution of the i-th h^2 value (the output of estimate_distributions)        
        estimated_values - A vector of size P, of estimated heritability values, for which we wish to 
                           calculate confidence intervals.
        confidence - The required confidence level of the CIs we wish to construct (e.g., 95%).
        seed - A seed for the random generator used for the randomized CIs, if needed.
        use_randomized_cis - Should we use randomized CIs, if needed.

    Returns:
        A matrix of size (P X 2), of the confidence intervals for each required estimate value.
```

:hamster:
