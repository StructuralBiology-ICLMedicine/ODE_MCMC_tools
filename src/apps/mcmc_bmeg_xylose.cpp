//#include "mcmc_multispecies_utils.h"

#define MODEL_DEF_HEADER_FILE "models/bmegModel_resource_v9_GFP_only_edit.h"
#define MODEL_CLASS_NAMESPACE MODEL1
#define MODEL_CLASS_NAME ode_model_functor
#include "model_class.h"
#undef MODEL_CLASS_NAME
#undef MODEL_CLASS_NAMESPACE
#undef MODEL_DEF_HEADER_FILE

#include "mvn_RAM_move_functor.h"
#include "utils/options_store.h"
#include "ode_likelihood_function.h"
#include "mcmc_ptr_algorithms.h"

#include "common_variables.h"

// #include "minimiser.h"


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <memory>




using std::cout;
using std::endl;
using std::cerr;
using std::ios;
using std::ifstream;
using std::ofstream;
using std::istream;


using std::string;



using namespace MODEL1;

int main(int argc, char **argv){


    typedef MCMC::ode_model_system_interface model_sys_intf_type;
    typedef MCMC::prior_functor_interface<MCMC::model_system_interface::vector_type> prior_intf_type;
    typedef MCMC::proposal_functor_interface<model_sys_intf_type::vector_type, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > proposal_intf_type;
    typedef MCMC::loglikelihood_functor_interface<model_sys_intf_type> ll_intf_type;

    typedef std::shared_ptr<model_sys_intf_type>  model_sys_intf_ptr_type;
    typedef std::shared_ptr<prior_intf_type>  prior_intf_ptr_type;
    typedef std::shared_ptr<proposal_intf_type>  proposal_intf_ptr_type;
    typedef std::shared_ptr<ll_intf_type>  ll_intf_ptr_type;

    std::shared_ptr<MCMC::prior_functor_interface<ode_model_functor::vector_type> > uninf_prior(new ODE::uninf_prior_functor<ode_model_functor::vector_type>());


    cout << "#"<<  endl <<  "# COMMAND LINE:" << endl;
    cout << "# ";
    for (int ii = 0; ii < argc; ii++){
        cout << argv[ii] << " ";
    }
    cout << endl << "#" << endl;

    UTILS::options_manager opt_mgr;

    opt_mgr.register_option("output_root",string(""),"root for output files");

    opt_mgr.register_option("random_seed", (unsigned long)315875, "random seed" );
    opt_mgr.register_option("random_seed:rand_init", (bool)false, "from random seed /dev/urandom" );

    opt_mgr.register_option("mcmc:steps", (unsigned long)(std::numeric_limits< unsigned long >::max()), "MCMC steps to run" );
    opt_mgr.register_option("mcmc:burn_in", (unsigned long)1, "MCMC burn in before outputting chain" );
    opt_mgr.register_option("mcmc:chain_output_freq", (unsigned long)1, "MCMC output frequency for chain" );
    opt_mgr.register_option("mcmc:inverse_temperature", (double)1, "MCMC inverse temperature" );
    opt_mgr.register_option("mcmc:min_adapt_step", (unsigned long)1, "minimum adaptation step for adaptive algorithms" );



    opt_mgr.register_option("help",bool(false),"output help");
    opt_mgr.register_option("output_param_bounds",bool(false),"output help");

    opt_mgr.register_option("proposal:mvn_sigma",string(""),"MVN sigma matrix for proposal distribution");
    opt_mgr.register_option("proposal:mvn_mapping",string(""),"parameter mapping for MVN sigma matrix for proposal distribution");
    opt_mgr.register_option("proposal:use_mvn",bool(false),"use input MVN for initial proposal distribution");
    opt_mgr.register_option("proposal:use_bounds",bool(false),"use information from the parameter bounds file to set up initial MVN proposal distribution");
    opt_mgr.register_option("proposal:indep_sampler",bool(false),"do uniform independence sampling");

    opt_mgr.register_option("params:bounds",string(""),"parameter bounds");
    opt_mgr.register_option("exp_data:list",string(""),"file containing location of experimental data for likelihood function");
    opt_mgr.register_option("likelihood:degf",double(5),"degrees of freedom setting if using multivariate Students's t distribution for likelihood function");
    opt_mgr.register_option("likelihood:student_t",bool(false),"use multivariate Students's t distribution for likelihood function");
    opt_mgr.register_option("params:initial",string(""),"initial parameters file");
    opt_mgr.register_option("params:random",bool(false),"initialise to random start");


    opt_mgr.register_option("prior:use_uninf",bool(false),"use uniform (uninformative) prior");
    opt_mgr.register_option("prior:use_mvn",bool(false),"use multivariate normal prior");
    opt_mgr.register_option("prior:mvn_means",string(""),"vector of means over parameters");
    opt_mgr.register_option("prior:mvn_sigma",string(""),"MVN sigma matrix for prior distribution");
    opt_mgr.register_option("prior:mvn_mapping",string(""),"MVN sigma matrix for prior distribution");


    std::vector<string> compulsory_flags;
    compulsory_flags.push_back("params:bounds");
    compulsory_flags.push_back("exp_data:list");
    compulsory_flags.push_back("params:initial");
    compulsory_flags.push_back("output_root");





    cout << "# parsing cmd line" << endl;
    if (opt_mgr.parse_cmd_line(argc, argv))
    {
        cout << "# parsed OK." << endl;
    }
    else {
        cerr << "# parsing errors" << endl;
        return -1;
    }

    //opt_mgr.list_options_full_format(cout);

    if (opt_mgr.get_option_value<bool>("help")){
        opt_mgr.list_options_full_format(cout);
        return 0;
    }




    if (opt_mgr.get_option_value<bool>("output_param_bounds")){
    	ode_model_functor mod;
    	std::vector<std::string> nvec = mod.get_joint_para_and_state_name_vec();
    	ode_model_functor::joint_para_x_y_state_type vals = mod.get_joint_para_and_state();
    	for (unsigned int ii = 0; ii < nvec.size(); ii++){
    		cout << nvec[ii] << "\t"
    				<< vals[ii] << "\t"
					<< 0 << "\t"
					<< 0 << "\t"
					<< endl;
    	}
    	return 0;
    }



    if (opt_mgr.is_set(""))
    {
        //! print flag-less options if they are set
        cerr << "# ERROR: You also provided disallowed free options - all options must have flags: "
            << opt_mgr.get_option_value<string>("")
            << endl;
        return -1;
    }

    bool compulsory_ok = true;
    for (std::vector<string>::const_iterator it = compulsory_flags.begin(); it != compulsory_flags.end(); it++){
        if (!opt_mgr.is_set(*it)){
            cerr << "ERROR: flag: --" << *it << " must be set" << endl;
            compulsory_ok = false;
        }
    }

    if (compulsory_ok == false){
        return -1;
    }


    uint64_t seed = opt_mgr.get_option_value<unsigned long>("random_seed");

    if (opt_mgr.get_option_value<bool>("random_seed:rand_init")){
    	cout << "# setting random seed from /dev/urandom" << std::endl;
    	ifstream f("/dev/urandom");
    	f.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    }

    const string output_root = opt_mgr.get_option_value<string>("output_root");

    const unsigned long steps = opt_mgr.get_option_value<unsigned long>("mcmc:steps");//100000000;
	const unsigned long burn_in = opt_mgr.get_option_value<unsigned long>("mcmc:burn_in");//1; //5000
	const unsigned long output_freq = opt_mgr.get_option_value<unsigned long>("mcmc:chain_output_freq");//1;  //10
    const double inv_temp = opt_mgr.get_option_value<double>("mcmc:inverse_temperature");//1;

    //const bool use_RAM = opt_mgr.get_option_value<bool>("mcmc:use_RAM");

    const double degf = opt_mgr.get_option_value<double>("likelihood:degf");
    const bool use_student_t = opt_mgr.get_option_value<bool>("likelihood:student_t");


    typedef std::shared_ptr<ode_model_functor> system_ptr_;
    system_ptr_ sys(new ode_model_functor());
	ode_model_functor::vector_type joint_state = sys->new_joint_para_and_state_type();

    const string prior_sigma_file = opt_mgr.get_option_value<string>("prior:mvn_sigma");//"uninf";
    const string prior_mean_file = opt_mgr.get_option_value<string>("prior:mvn_means");//"uninf";
    const string prior_mapping_file = opt_mgr.get_option_value<string>("prior:mvn_mapping");//"uninf";
    const bool using_mvn_prior = opt_mgr.get_option_value<bool>("prior:use_mvn");
    const bool using_uninf_prior = opt_mgr.get_option_value<bool>("prior:use_uninf");

    if ((!using_mvn_prior && !using_uninf_prior) || (using_mvn_prior && using_uninf_prior)){
        cerr << "ERROR: please set one of the --prior:use_mvn or --prior:use_uninf flags but not both" << endl;
        return -1;
    }

    //seed = lexical_cast<uint64_t>(argv[1]);
    cout << "# seed:\t" << seed << endl;




    const bool using_mvn_prop = opt_mgr.get_option_value<bool>("proposal:use_mvn");
    const bool using_bounds_prop = opt_mgr.get_option_value<bool>("proposal:use_bounds");
    const bool do_indep_sample = opt_mgr.get_option_value<bool>("proposal:indep_sampler");


    if ((!using_mvn_prop && !using_bounds_prop) || (using_mvn_prop && using_bounds_prop)){
        cerr << "ERROR: please set one of the --proposal:use_mvn or --proposal:use_bounds flags but not both" << endl;
        return -1;
    }

    const string sigma_file = opt_mgr.get_option_value<string>("proposal:mvn_sigma");//string(argv[2]);
    const string move_mapping = opt_mgr.get_option_value<string>("proposal:mvn_mapping");//string(argv[3]);
    const string mover_bounds_file = opt_mgr.get_option_value<string>("params:bounds");//string(argv[4]);


	//ODE::ode_loglikelihood_mvn_functor<ode_model_functor> ll_functor;
	std::shared_ptr<ODE::ode_loglikelihood_mvn_ptr_functor<MCMC::ode_model_system_interface> > ll_functor(new ODE::ode_loglikelihood_mvn_ptr_functor<MCMC::ode_model_system_interface>());

	//sim_condition_map conds_map;
	//double_vec times;

	const string data_list_file = opt_mgr.get_option_value<string>("exp_data:list");//string(argv[5]);
	const string param_file = opt_mgr.get_option_value<string>("params:initial");//string(argv[6]);
	const bool random_start_params = opt_mgr.get_option_value<bool>("params:random");//string(argv[6]);







    //output_freq = lexical_cast<unsigned long>(argv[7]);
    //inv_temp = lexical_cast<double>(argv[8]);

    cout << "# output_freq: " << output_freq << " inv_temp: " << inv_temp << endl;


	if (using_mvn_prior){
	    //cout << "# reading MVN prior from files: mean: " << prior_mean_file << " sigma: " << prior_sigma_file << " mapping: " << prior_mapping_file << endl;
        if (opt_mgr.is_set("prior:mvn_means") && opt_mgr.is_set("prior:mvn_sigma") && opt_mgr.is_set("prior:mvn_mapping")){
            cout << "# reading MVN prior from files: mean: " << prior_mean_file << " sigma: " << prior_sigma_file << " mapping: " << prior_mapping_file << endl;
        }
        else {
            cerr << "ERROR: the --prior:mvn_mean --prior:mvn_sigma and --prior:mvn_mapping flags need to be set with the --prior:use_mvn option." << endl;
            return -1;
        }
	}








	cout << "# parsing input data list file" << endl;
	bool parseOK = ODE::parse_MVN_list_file<ode_model_functor>(data_list_file, ll_functor->conds, ll_functor->times, degf);
	if (parseOK == false){
        cerr << "ERROR: problem parsing file: " << data_list_file << "\n Depending on the error, you probably need to check your file paths" << endl;
        exit(EXIT_FAILURE);
	}

	if (use_student_t){
		ll_functor->use_t_dist = true;
		cout << "# Using multivariate Student's t distribution for the likelihood function with degrees of freedom = " << degf << endl;
	}
	else {
		ll_functor->use_t_dist = false;
		cout << "# Using multivariate normal distribution for the likelihood function" << endl;
	}

    cout << "# number of time steps: " << ll_functor->times.size() << endl;

	COMMON::rand_gen.seed(seed);
	// burn in random num generator
	for (int ii = 0; ii < 1000; ii++){
		COMMON::rand_gen.randInt();
	}

	auto rand_ptr = std::shared_ptr<MTRand>(new MTRand());
    rand_ptr->seed(seed);
	COMMON::rand_gen.seed(seed);
	// burn in random num generator
	for (int ii = 0; ii < 1000; ii++){
		COMMON::rand_gen.randInt();
		rand_ptr->randInt();
	}

	cout << "# parsing parameters" << endl;
	ODE::parse_initial_params(param_file, joint_state);
	sys->set_joint_para_and_state(joint_state);



	std::shared_ptr<RAM_move_parameters_functor<ode_model_functor> > ram_move_functor(new RAM_move_parameters_functor<ode_model_functor>(sys, rand_ptr));


	ram_move_functor->max_adapt_step = std::numeric_limits< unsigned long>::max();//1000;
	cout << "# INFO: MAX_ADAPT_STEP = " <<  ram_move_functor->max_adapt_step;

    cout << "# reading parameter bounds for mover_functor" << endl;
    ram_move_functor->loadBounds(*sys, mover_bounds_file);


    if (do_indep_sample){
    	cout << "# doing uniform independence sampling" << endl;
    	ram_move_functor->set_up_default_covar();
    	ram_move_functor->do_uniform_random = do_indep_sample;
    }
    else if (using_mvn_prop){
        cout << "# loading MVN proposal sigma from command line" << endl;
        //use_RAM ? ram_move_functor->use_mvn_mover = true : move_functor->use_mvn_mover = true;
        cout << "# reading MVN sigma for mover_functor: " <<  sigma_file << endl;
        ram_move_functor->loadSigma(sigma_file);
        cout << "# loading move_mapping for mover_functor: "<<  move_mapping << endl;
        ram_move_functor->load_mapping(*sys, move_mapping);
	}
	else {
        cout << "# using parameter bounds file to set up initial MVN mover" << endl;
        ram_move_functor->set_up_default_covar();
        //use_RAM ? ram_move_functor->use_mvn_mover = true : move_functor->use_mvn_mover = true;//false;
	}

    if (random_start_params){
    	cout << "# setting start parameters to random" << endl;
    	joint_state = ram_move_functor->uniform_random(joint_state);
    	sys->set_joint_para_and_state(joint_state);
    }





    if (ram_move_functor->check_bounds(joint_state, true) == false){
        cerr << "ERROR: check your starting parameter values - they are not within the defined bounds" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "# starting MCMC" << endl;
    if (using_mvn_prior == false){


    	double ll = ll_functor->operator()(sys);
    	cout << "# Initial LL: " << ll << endl;

    	cout << "# using RAM adaptive mover" << endl;
    	ram_move_functor->target_acc_rate = 0.234;//0.3;
    	//model_sys_intf_ptr_type sys_base = sys;
    	ll_intf_ptr_type ll_base = ll_functor;

    	//minimise(sys, ll_functor, uninf_prior, ram_move_functor, 10000);

    	auto mcmc_functor =  MCMC::PTR::adaptive_mcmc_interface<model_sys_intf_ptr_type, ll_intf_ptr_type,
    			prior_intf_ptr_type, proposal_intf_ptr_type >(
    					sys,
						ll_functor, uninf_prior, ram_move_functor,
						output_root, steps, burn_in, output_freq, inv_temp);
    	mcmc_functor.greedy_start = false;
    	mcmc_functor.min_adaptation_step = opt_mgr.get_option_value<unsigned long>("mcmc:min_adapt_step");
    	mcmc_functor(sys, rand_ptr, cout, cerr);




    }
    else {
    	cout << "# parsing MVN prior files" << endl;
    	Eigen::Matrix<double, Eigen::Dynamic, 1> prior_mean = readDynamicMatrix(prior_mean_file);
    	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> prior_sigma = readDynamicMatrix(prior_sigma_file);
    	std::shared_ptr<ODE::multigauss_prior_functor<MCMC::ode_model_system_interface> > prior_functor(new ODE::multigauss_prior_functor<MCMC::ode_model_system_interface>() );
    	prior_functor->load_mapping(*sys, prior_mapping_file);

    	prior_functor->set_covar(prior_sigma);
    	prior_functor->set_mean(prior_mean);
    	bool priorOK = prior_functor->verify();
    	if (priorOK == false){
    		throw std::exception();
    	}

    	double ll = ll_functor->operator()(sys);
    	cout << "# Initial LL: " << ll << endl;
    	cout << "# Initial prior: " << prior_functor->operator()(joint_state) << endl;

    	cout << "# using RAM adaptive mover" << endl;
    	ram_move_functor->target_acc_rate = 0.234;//0.3;
    	prior_intf_ptr_type prior_base = prior_functor;
    	auto mcmc_functor =  MCMC::PTR::adaptive_mcmc_interface<model_sys_intf_ptr_type, ll_intf_ptr_type,
    			prior_intf_ptr_type, proposal_intf_ptr_type >(
    					sys,
						ll_functor, prior_functor, ram_move_functor,
						output_root, steps, burn_in, output_freq, inv_temp);
    	mcmc_functor.greedy_start = false;
    	mcmc_functor.min_adaptation_step = opt_mgr.get_option_value<unsigned long>("mcmc:min_adapt_step");
    	mcmc_functor(sys, rand_ptr, cout, cerr);

    }
}




