//#include "mcmc_multispecies_utils.h"
#include "ode_likelihood_function.h"

#define MODEL_DEF_HEADER_FILE "models/bmegModel_resource_v9_GFP_only_edit.h"
#define MODEL_CLASS_NAMESPACE MODEL1
#define MODEL_CLASS_NAME ode_model_functor
#include "model_class.h"
#undef MODEL_CLASS_NAME
#undef MODEL_CLASS_NAMESPACE
#undef MODEL_DEF_HEADER_FILE

//#include "mvn_AM_move_functor.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "utils/options_store.h"

using namespace MODEL1;
using namespace std;
using namespace ODE;

int main(int argc, char **argv){


    UTILS::options_manager opt_mgr;
    opt_mgr.register_option("times",string(""),"times");
    opt_mgr.register_option("params",string(""),"single parameters");
    opt_mgr.register_option("cond",string(""),"condition definition");

    opt_mgr.register_option("cond_species_only",bool(false),"output only the species defined in the condition file");

    opt_mgr.register_option("help",bool(false),"output help");

    std::vector<string> compulsory_flags;
    compulsory_flags.push_back("times");
    compulsory_flags.push_back("params");
    //compulsory_flags.push_back("cond");

    //cout << "# parsing cmd line" << endl;
    if (opt_mgr.parse_cmd_line(argc, argv))
    {
    	//cout << "# parsed OK." << endl;
    }
    else {
    	cerr << "# parsing errors" << endl;
    	return -1;
    }

    if (opt_mgr.get_option_value<bool>("help")){
    	opt_mgr.list_options_full_format(cout);
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

    if (!opt_mgr.is_set("cond") && opt_mgr.is_set("cond_species_only")){
    	cerr << "ERROR: can't use --cond_species_only option without also defining --cond";
    	return -1;
    }

    const bool output_cond_species = opt_mgr.get_option_value<bool>("cond_species_only");


    ode_model_functor::joint_para_x_y_state_type joint_state = ode_model_functor::new_joint_para_x_y_state_type();
	vector<double> times = {0, 300};

	string times_file = opt_mgr.get_option_value<string>("times");//string(argv[1]);
	string param_file = opt_mgr.get_option_value<string>("params");//string(argv[2]);
	string cond_file = opt_mgr.get_option_value<string>("cond");//string(argv[3]);

    times = ODE::readTimes(times_file);

	ODE::parse_initial_params(param_file, joint_state);

	ode_model_functor sys(joint_state);


    unsigned int species = NAME_NOT_FOUND;
    double scale = 1.0;
    double min_variance = 0;
    double weight = 1;

    if (opt_mgr.is_set("cond")){
    	sim_condition this_cond = parse_condition_spec_file<ode_model_functor>(cond_file, species, scale, min_variance, weight);
    	apply_condition<ode_model_functor>(sys, this_cond);
    }


	vector<typename ode_model_functor::joint_x_y_state_type> traj = ODE::simple_run_trajectory(sys,  times);

    for (vector<typename ode_model_functor::joint_x_y_state_type>::const_iterator it = traj.begin(); it != traj.end(); it++){
            const unsigned int index = it - traj.begin();
            std::cout << times[index];
            if (!output_cond_species){
            	for (ode_model_functor::joint_x_y_state_type::const_iterator sit = it->begin(); sit != it->end(); sit++){
            		std::cout << "\t" << *sit;
            	}
            }
            else {
            	std::cout << "\t" << (((*it)[species])/scale);
            }
            std::cout << std::endl;
    }


}
