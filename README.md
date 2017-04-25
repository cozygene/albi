## FIESTA (*F*ast confidence *I*nt*E*rvals using *ST*ochastic *A*pproximation) 

FIESTA is a method for the construction of accurate confidence intervals (CIs) for heritability. FIESTA can be used as an add-on to existing methods for heritability and variance components estimation. FIESTA is described in the following [paper]().

The Supplemental Information for the FIESTA RECOMB 2017 paper is available [here]().

If you want to actually estimate the distribution of the heritability estimator, or need CIs for a lot of estimates, check out FIESTA's predecessor, ALBI, available under the same package [here](ALBI.md).

FIESTA's input is the eigenvalues of a kinship matrix (possibly, also its eigenvectors, and general covariates), and it produces accurate confidence intervals for a set of heritability estimates.

Code written by Regev Schweiger, Eyal Fisher and Reut Yedidim. For questions or comments, mail [schweiger@post.tau.ac.il](mailto:schweiger@post.tau.ac.il).

### Simple Example

The following reads the eigenvalues from a file, calculates 95% CIs over a grid of heritability estimates, and outputs it:

```
   python fiesta.py --kinship_eigenvalues ./data/eigenvalues.txt
```
**Note**: Under Windows, replace ```fiesta.py``` with ```albi_lib/fiesta.py```.

The default output is the CIs for heritability estimates of 0, 0.1, ..., 0.9, 1:

```
Estimate   CI_lower_bound  CI_upper_bound
0.00000	0.00000	0.69230
0.10000	0.00000	0.81158
0.20000	0.00000	0.90069
0.30000	0.00000	1.00000
0.40000	0.00000	1.00000
0.50000	0.03170	1.00000
0.60000	0.07461	1.00000
0.70000	0.17200	1.00000
0.80000	0.21739	1.00000
0.90000	0.33654	1.00000
1.00000	0.42484	1.00000
```

For example, the CI for an estimate of 0.1 is [0, 0.81158].

If you use covariates or you generated your matrix in a non-standard way, you would need additional flags. For more information and options, see below.

## Installation

There are several ways to install and use FIESTA.

<details>
### 1. Download ZIP

You can simply download the code and run it as a standard program:

1. Download the ZIP file of the latest release here: [https://github.com/cozygene/albi/releases](https://github.com/cozygene/albi/releases).
2. Extract the ZIP file to a folder
3. Run FIESTA from the folder:
```
   python fiesta.py [flags]
```

**Note**: The [NumPy](http://www.numpy.org/) and SciPy packages are required for FIESTA.

### 2. Install using `pip` or `setuptools`

If you have `pip` installed, run this on command line (from anywhere): 
```
   sudo pip install git+https://github.com/cozygene/albi
```

Alternatively, you can also clone the repository and do a manual install:
```
   git clone https://github.com/cozygene/albi
   sudo python setup.py install
```

Note: If you are using anaconda, you could replace `sudo pip` and `sudo python` with the path to the binaries in the `bin` directory in your anaconda installation.

### Uninstalling

To uninstall, run:
```
    sudo -H pip uninstall albi
```    
</details>


## Running FIESTA

The following section describes all of FIESTA's options and flags.

There are two use cases, and it is important to distinguish between them: The case where all covariates are eigenvectors of the kinship matrix (this often includes the intercept!); and the case of general covariates. The first case allows for simpler and faster computation, while the second is more general.

#### 1. All covariates are eigenvectors 

Where all covariates are eigenvectors of the kinship matrix, there is a simple closed-form formula that produces faster computations. Another nice side effect is that, in this case, the eigenvectors themselves are not needed; only the eigenvalues, and the identity of the eigenvectors that are used as covariates, is needed. In particular, this includes the following special cases:

1. Only an **intercept term** (a vector of ones). This is the default in e.g., GCTA. Since the genotype matrix is mean-centered, it follows that the vector of ones is an eigenvector, corresponding to the last (smallest) eigenvector of the kinship matrix (see Supplemental Note S4.5 in the ALBI paper). **If you built your kinship matrix with default settings, this is probably the case relevant to you**.
2. **No fixed effects** at all (including an intercept). 
3. In addition to an intercept term, a common practice is adding as fixed effects the q **largest principal components**, corresponding to the first q eigenvectors of the kinship matrix.

The following command builds CIs and saves them to a file:
```
   python fiesta.py --kinship_eigenvalues filename
                 [--use_eigenvectors_as_covariates <list of #>]
                 (--estimates_filename filename
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                 [--iterations <# of iterations>] 
                 [--output_filename filename]
```
**Note**: Under Windows, replace ```fiesta.py``` with ```albi_lib/fiesta.py```.

The flags are as follows:

* `kinship_eigenvalues` (Shortcut: `-k`) - A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format. This could be created, for example, with GCTA's `--pca` flag.
* `use_eigenvectors_as_covariates` (`-u`) - A list detailing which eigenvectors should be used as covariates. For example: If you only use an intercept term, and your kinship matrix was constructed from a mean-centered SNP matrix, use `-1` (the last eigenvector is equal to the constant vector); if in addition, you add the first 3 PCs as covariates, use `0,1,2,-1`. Default is `-1`.
* `estimates_filename` (`-f`) - A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one. 
* `estimate_grid` (`-g`) - Alternatively, one can ask ALBI to calculate CIs for a grid of heritability estimates (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1). Default is 10.
* `confidence` (`-c`) - The required confidence level for the CIs. Default is 0.95 (95% CIs).
* `iterations` (`-n`) - The number of iterations used int the algorithm (see paper). Default is 1000.
* `output_filename` (`-o`) - File to which to write the calculated CIs. If not supplied, CIs will be printed.

#### 2. General covariates

When covariates are not eigenvectors of the kinship matrix, you must supply both the eigenvalues, the eigenvectors and the covariates explicitly. Computation will be slower (although still linear in the number of individuals):

The following command estimates the distributions, and saves the results to a file:
```
   python fiesta.py --kinship_eigenvalues filename
                  --kinship_eigenvectors filename
                 [--covariates filename]
                 [--add_intercept True/False]
                 (--estimates_filename filename
                    or
                  --estimate_grid <# of grid points>)
                 [--confidence <required confidence level>] 
                 [--iterations <# of iterations>] 
                 [--output_filename filename]
```
**Note**: Under Windows, replace ```fiesta.py``` with ```albi_lib/fiesta.py```.

The flags are as follows:

* `kinship_eigenvalues` (`-k`) - A file containing the eigenvalues of the kinship matrix, one eigenvalue per line, in text format. This could be created, for example, with GCTA's `--pca` flag.
* `kinship_eigenvectors`(`-v`) - A file containing the eigenvectors of the kinship matrix, one eigenvector per column, in text format. This could be created, for example, with GCTA's `--pca` flag.
* `covariates` (`-x`) - A file containing the covariates, one covariate per column, in text format. 
* `add_intercept` (`-i`) - Boolean flag (True/False) - whether to use a constant 1 (intercept) covariate. Can be used without an additional covariates file. Default is True.
* `estimates_filename` (`-f`) - A filename containing a list of heritability estimates (one per line) in text format. A CI will be calculated for each one. 
* `estimate_grid` (`-g`) - Alternatively, one can ask ALBI to calculate CIs for a grid of heritability estimates (e.g., a grid of 100, will calculate CIs for 0, 0.01, ..., 0.99, 1). Default is 10.
* `confidence` (`-c`) - The required confidence level for the CIs. Default is 0.95 (95% CIs).
* `iterations` (`-n`) - The number of iterations used int the algorithm (see paper). Default is 1000.
* `output_filename` (`-o`) - File to which to write the calculated CIs. If not supplied, CIs will be printed.

:hamster:
