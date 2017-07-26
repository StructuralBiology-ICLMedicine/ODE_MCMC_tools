//#include "mcmc_multispecies_utils.h"
#include "ode_likelihood_function.h"

#define MODEL_DEF_HEADER_FILE "models/bmegModel_resource_v9_single_deg_edit.h"
#define MODEL_CLASS_NAMESPACE MODEL_ORIG
#define MODEL_CLASS_NAME ORIG_ode_model_functor
#include "model_class.h"
#undef MODEL_CLASS_NAME
#undef MODEL_CLASS_NAMESPACE
#undef MODEL_DEF_HEADER_FILE

#define MODEL_DEF_HEADER_FILE "models/bmegModel_resource_v9_single_deg_track_edit.h"
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
using namespace MODEL_ORIG;
using namespace std;
using namespace ODE;

ode_model_functor::joint_para_x_y_state_type convert_params(ORIG_ode_model_functor::joint_para_x_y_state_type orig_joint_state){
    ode_model_functor::joint_para_x_y_state_type joint_state = ode_model_functor::new_joint_para_x_y_state_type();
    ORIG_ode_model_functor old_mod;
    std::vector<std::string> name_vec = old_mod.get_joint_para_and_state_name_vec();
    for (unsigned int ii = 0 ; ii < orig_joint_state.size(); ii++){
    	string this_name = name_vec[ii];
    	unsigned int new_index = ode_model_functor::get_joint_para_x_y_state_index(this_name);
    	joint_state[new_index] = orig_joint_state[ii];
    }

    return joint_state;
}


typedef std::vector< ode_model_functor::joint_para_x_y_state_type > para_type_vector;
typedef std::vector< ORIG_ode_model_functor::joint_para_x_y_state_type > orig_para_type_vector;

int main(int argc, char **argv){


    UTILS::options_manager opt_mgr;
    opt_mgr.register_option("times",string(""),"times");
    opt_mgr.register_option("param_list",string(""),"list of parameters (MCMC line output format)");
    opt_mgr.register_option("cond",string(""),"condition definition");
    opt_mgr.register_option("species",string(""),"use this species instead of one defined in condition file");

    //opt_mgr.register_option("cond_species_only",bool(false),"output only the species defined in the condition file");

    opt_mgr.register_option("help",bool(false),"output help");

    std::vector<string> compulsory_flags;
    compulsory_flags.push_back("times");
    compulsory_flags.push_back("param_list");
    compulsory_flags.push_back("cond");

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


	vector<double> times = {0, 300};

	string times_file = opt_mgr.get_option_value<string>("times");//string(argv[1]);
	string param_list_file = opt_mgr.get_option_value<string>("param_list");//string(argv[2]);
	string cond_file = opt_mgr.get_option_value<string>("cond");//string(argv[3]);

    times = ODE::readTimes(times_file);

	//ODE::parse_initial_params(param_file, joint_state);
    orig_para_type_vector orig_para_vec = ODE::parse_parameters_list<ORIG_ode_model_functor>(param_list_file);
    para_type_vector para_vec;
    for (auto  it = orig_para_vec.begin(); it != orig_para_vec.end(); it++ ){
    	ode_model_functor::joint_para_x_y_state_type new_pars = convert_params(*it);
    	para_vec.push_back(new_pars);
    }

    unsigned int species = NAME_NOT_FOUND;
    double scale = 1.0;
    double min_variance = 0;
    double weight = 1;

    sim_condition this_cond = parse_condition_spec_file<ode_model_functor>(cond_file, species, scale, min_variance, weight);

    for (auto it = times.begin(); it != times.end() ; it++){
    	if (times.end() - it > 1){
    		std::cout << *it << "\t";
    	}
    	else{
    		std::cout << *it;
    	}
    }
    std::cout << std::endl;

    if (opt_mgr.is_set("species")){
    	string spec_name = opt_mgr.get_option_value<string>("species");
    	species = ode_model_functor::get_joint_x_y_state_index(spec_name);
    	scale = 1;
    	if (species == NAME_NOT_FOUND){
    		cerr << "ERROR: species label " << spec_name << " in --species not found" << endl;
    		return -1;
    	}
    }


    for (auto it = para_vec.begin(); it != para_vec.end(); it++ ){
    	ode_model_functor::joint_para_x_y_state_type joint_state = *it;//ode_model_functor::new_joint_para_x_y_state_type();
    	ode_model_functor sys(joint_state);

    	apply_condition<ode_model_functor>(sys, this_cond);

    	vector<typename ode_model_functor::joint_x_y_state_type> traj = ODE::simple_run_trajectory(sys,  times);

    	for (vector<typename ode_model_functor::joint_x_y_state_type>::const_iterator it = traj.begin(); it != traj.end(); it++){
    		//const unsigned int index = it - traj.begin();
    		//std::cout << times[index];
    		if ((traj.end() - it) > 1 ){
    			std::cout << (((*it)[species])/scale) << "\t";
    		}
    		else {
    			std::cout << (((*it)[species])/scale) ;
    		}
    	}
    	std::cout << std::endl;
    }


}
