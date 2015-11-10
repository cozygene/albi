
# ALBi // TODO !!
change aruments long name <br/>
change help and documentation
albi is everything you need

## Installation

Three installation options:

### 1. Download ZIP

a. Press the 'Download ZIP' button on the right side of the screen
b. Extract the ZIP file to a albi's folder
c. Run on command line (from the folder):
```
   sudo python setup.py install
```

### 2. or Install using pip

If you have pip installed, run this on command line (from anywhere): 
```
   sudo pip install git+https://github.com/cozygene/albi
```

### 3. or Clone folder

Run this on command line (from anywhere): 
You can also clone the repository and do a manual install:
```
   git clone https://github.com/cozygene/albi
   sudo python setup.py install
```


 * To uninstal run 'sudo pip uninstall albi' <br/>
** After installing you can delete the cloned folder.
TODO maybe remove this


## Running

### Command Line

Run from anywhere:
```
albi.py ---kinship_eigenvalues [file] --estimate_grid [number]
```
 
### Importing albi

```
    >>> from albi import calcCI
    >>> from albi import buildHeritability
    >>> from albi import run # combines both calcCI and buildHeritability
    >>> calcCI(...)

```
 
