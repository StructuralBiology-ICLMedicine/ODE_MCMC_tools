# Examples


## Running standard MCMC on the promoter competition data

```
../../../bin/mcmc_bmeg_competition --params:bounds parameter_bounds.dat --exp_data:list experimental_inputs.dat --params:initial start_params.dat --prior:use_uninf --mcmc:chain_output_freq 10 --mcmc:inverse_temperature 1  --mcmc:burn_in 1  --proposal:use_mvn --proposal:mvn_mapping mapping.dat --proposal:mvn_sigma ProposalSigma.dat   --random_seed:rand_init --likelihood:student_t --likelihood:degf 5 --output_root OUTPUT
```

Where ProposalSigma.dat is the multivariate normal (MVN) proposal distribution covariance matrix, mapping.dat describes the mapping of columns/rows the ProposalSigma.dat to parameter names or initial values, experimental\_inputs.dat is a file listing the experimental conditions and corresponding experimental data files (in the form of a vector of means and a covariance matrix), start\_params.dat is a vector describing the initial start parameters. Descriptions of the options can be found by using the --help command line option.

## Running parallel tempering MCMC on the promoter competition data

```
../../../bin/ptmcmc_bmeg_competition --params:bounds parameter_bounds.dat --exp_data:list experimental_inputs.dat --params:initial_list start_params_list.dat --prior:use_uninf --mcmc:chain_output_freq 10 --mcmc:inverse_temperature 1 --threads 16 --pt:ladder_size 16 --mcmc:burn_in 1000  --proposal:use_mvn --proposal:mvn_mapping mapping.dat --proposal:mvn_sigma ProposalSigma.dat --random_seed:rand_init --sim_mgr:exchange_freq 500 --sim_mgr:job_chunks 500 --likelihood:student_t --likelihood:degf 5  --output_root PT_OUTPUT
```

Where the input files are the same as for the mcmc\_bmeg\_competition command, but with extra options to control the parallel tempering and multithreading behaviour. An arbitrary number of threads can be made available to the program but they will not all be used if the temperature ladder is smaller. The --sim\_mgr:job\_chunks argument should be set to the same as the --sim\_mgr:exchange\_freq argument, which controls how often exchanges are attempted between different temperatures. The output file PT\_OUTPUT\_beta\_0\_model\_0\_std\_mcmc\_chain.dat will contain the MCMC chain at thermodynamic beta (inverse temperature) 1, will the mcmc\_chain.dat output files containing progressively higher temperature chains (i.e. lower inverse temperatures).


## License

This software is distributed under the GNU GPL license, version 3.

(C) James T. MacDonald, 2017. 
Imperial College London.


