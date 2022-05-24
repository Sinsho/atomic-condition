# Detecting Floating-Point Errors via Atomic Conditions MultiVar Extension
## Disclaimer
This repository is forked from [here](https://github.com/FP-Analysis/atomic-condition) and extended to work with
multivariate functions. This work changes fundamentally how functions are read and processed and as such is not
easily merged to the original repository.


## Setup
Install the [docker](https://www.docker.com/). Here are some guides for install [docker for Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/) and [docker for MacOS](https://docs.docker.com/docker-for-mac/install/)

Clone this repo to your local workspace:
```
$ git clone https://github.com/Sinsho/atomic-condition
```

The `docker` folder contains a Dockerfile that can be used to build our implementation. If an error about permission denied occurs, try to run docker in root `sudo docker ...`.
```
$ cd atomic-condition/docker
$ docker build -t atomic .
```
It takes around 40 minutes for installing and building all necessary packages and compiling the whole GSL Library.

## Usage
Run this docker's container with interactive mode, the working directory is at `/atom`.
```
$ docker run -it atomic /bin/bash
```

### Run on GSL Functions
```
$ make
$ bin/gslSolver.out gsl <function_index>
```

#### Compute Relative Error (Only support GSL functions for now)
Using the oracle from `mpmath` to compute the relative error:
```
$ make
$ bin/gslSolver.out gsl <function_index> && python3 script/oracleMpmathMultVar.py
```

Example: GSL function `gsl_sf_gegenpoly_3` has the index 7,
```
$ bin/gslSolver.out gsl 7 && python3 script/oracleMpmath.py
...
Function Index: 7
Max Relative Error:
  Input: -6.1512035535082410e-01, 1.0407336901334387e+00, 
  Output: -1.4214773298447788e-16
  Oracle: -1.8897179271988252e-16
        Relative Error: 2.47783e-01
```

### Reproduce the Results
#### Getting the relative errors of GSL functions
You can check GSL functions one by one as mentioned in above section.
You can also check all of the results in single run with:
```
$ make
$ bin/gslSolver.out gsl all
$ python3 script/oracleMpmathMultVar.py
```