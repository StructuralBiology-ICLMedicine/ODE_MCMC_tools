#!/bin/sh
#   qsub  -q james  -lnodes=1:ppn=1:james -v PAR=start_params_v9_single_deg_run.dat   submit_bmeg_mcmc.sh
cd $PBS_O_WORKDIR
#         echo `pwd`

echo using par $PAR



/csb/home/jmacdona/development/ode_mcmc_code/bin/mcmc_bmeg_v9_GFP_only --params:bounds parameter_bounds_v9_GFP_only_xyl.dat --exp_data:list XylRXylose_ts25_Mat_GFP_list_of_inputs.dat --params:initial $PAR --prior:use_mvn --prior:mvn_sigma PriorSigma.dat --prior:mvn_means PriorMeans.dat --prior:mvn_mapping prior_mapping_v9 --mcmc:chain_output_freq 10 --mcmc:inverse_temperature 1  --mcmc:burn_in 10000  --proposal:use_mvn --proposal:mvn_mapping mapping_v9_xyl --proposal:mvn_sigma Sigma_v9_xyl.dat   --random_seed:rand_init --likelihood:student_t --likelihood:degf 5 --output_root prod_v9_GFP_only_xyl_job_$PBS_JOBID

#/csb/home/jmacdona/development/ode_mcmc_code/bin/mcmc_bmeg_v9_single_deg --params:bounds parameter_bounds_v9_single_deg_comp.dat --exp_data:list all_msd05_inputs.dat --params:initial $PAR --prior:use_uninf --mcmc:chain_output_freq 10 --mcmc:inverse_temperature 1  --mcmc:burn_in 1  --proposal:use_mvn --proposal:mvn_mapping mapping_v9_single_deg --proposal:mvn_sigma Sigma_v9_new.dat   --random_seed:rand_init --likelihood:student_t --likelihood:degf 5 --output_root prod_comp_test_v9_single_deg_job_$PBS_JOBID

