# ODE\_MCMC\_tools

ODE\_MCMC\_tool provides a C++ framework for Markov chain Monte Carlo (MCMC) and parallel tempering MCMC with ordinary differential equation models.

## Requirements

The software has currently been tested with Linux (RHEL 7) and Mac OS X with MacPorts (without the static compile flag). The software requires the Ceres Solver, Eigen3 and Boost C++ libraries (including the development header files). The build system currently uses SCons.

## Compiling the C++ code

To compile the software go into the main directory and run the command:
```
scons -j1 
```

Where the -j option specifies the number of simultaneous tasks you wish to use to speed up compilation.

## Running MCMC on the example data 

Some example ODE models and real experimental data are provided from some Bacillus megaterium cell-free transcription and translation experiments.

Go into the example\_data/bmeg/xylose\_xylr or example\_data/bmeg/competition\_experiment directories and follow the example command lines provided in the README.md files.

## License

This software is distributed under the GNU GPL license, version 3.

(C) James T. MacDonald, 2017. 
Imperial College London.





