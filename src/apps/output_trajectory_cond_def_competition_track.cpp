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
#include <boost/algorithm/string/split.hpp>


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

int main(int argc, char **argv){


    UTILS::options_manager opt_mgr;
    opt_mgr.register_option("times",string(""),"times");
    opt_mgr.register_option("params",string(""),"single parameters");
    opt_mgr.register_option("cond",string(""),"condition definition");

    opt_mgr.register_option("cond_species_only",bool(false),"output only the species defined in the condition file");
    opt_mgr.register_option("output_species_list",string(""),"comma separated list of species names");

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
    const bool using_species_list = opt_mgr.is_set("output_species_list");

    if (output_cond_species && using_species_list){
    	cerr << "ERROR: can't use --cond_species_only option simultaneously with the --output_species_list option" << endl;
    	return -1;
    }

    //ode_model_functor::joint_x_y_state_type output_mapping = ode_model_functor::new_joint_x_y_state_type();
    vector<unsigned int> species_order;
    /*
    for (unsigned int ii = 0 ; ii < output_mapping.size(); ii++){
    	output_mapping[ii] = using_species_list ? 0 : 1;
    }
    */
    if (using_species_list){
    	string csv_spec_list = opt_mgr.get_option_value<string>("output_species_list");
    	string_vector species_list;
    	boost::split( species_list, csv_spec_list, boost::is_any_of(",") );
    	for (unsigned int ii = 0 ; ii < species_list.size(); ii++){
    		const unsigned int sp_index = ode_model_functor::get_joint_x_y_state_index(species_list[ii]);
    		if (sp_index != NAME_NOT_FOUND){
    			//output_mapping[sp_index] = 1;
    			species_order.push_back(sp_index);
    		}
    		else {
    			cerr << "ERROR: species label " << species_list[ii] << " in --output_species_list not found" << endl;
    			return -1;
    		}
    	}
    }






   //ode_model_functor::joint_para_x_y_state_type joint_state = ode_model_functor::new_joint_para_x_y_state_type();
	vector<double> times = {0, 300};

	string times_file = opt_mgr.get_option_value<string>("times");//string(argv[1]);
	string param_file = opt_mgr.get_option_value<string>("params");//string(argv[2]);
	string cond_file = opt_mgr.get_option_value<string>("cond");//string(argv[3]);

    times = ODE::readTimes(times_file);

    ORIG_ode_model_functor::joint_para_x_y_state_type orig_joint_state = ORIG_ode_model_functor::new_joint_para_x_y_state_type();
	ODE::parse_initial_params(param_file, orig_joint_state);
	ode_model_functor::joint_para_x_y_state_type joint_state = convert_params(orig_joint_state);

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

            	for (vector<unsigned int>::const_iterator ui_it = species_order.begin(); ui_it != species_order.end(); ui_it++){
            		const unsigned int sp_index = *ui_it;
            		std::cout << "\t" << (*it)[sp_index];
            	}
            	/*
            	for (ode_model_functor::joint_x_y_state_type::const_iterator sit = it->begin(); sit != it->end(); sit++){
            		const unsigned int sp_index = sit - it->begin();
            		if (output_mapping[sp_index] > 0){
            			std::cout << "\t" << *sit;
            		}
            	}
            	*/
            }
            else {
            	std::cout << "\t" << (((*it)[species])/scale);
            }
            std::cout << std::endl;
    }


}
