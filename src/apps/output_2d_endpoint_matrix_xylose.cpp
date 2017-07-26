//#include "mcmc_multispecies_utils.h"
#include "ode_likelihood_function.h"



#define MODEL_DEF_HEADER_FILE "models/bmegModel_resource_v9_GFP_only_edit.h"
#define MODEL_CLASS_NAMESPACE MODEL_ORIG
#define MODEL_CLASS_NAME ORIG_ode_model_functor
#include "model_class.h"
#undef MODEL_CLASS_NAME
#undef MODEL_CLASS_NAMESPACE
#undef MODEL_DEF_HEADER_FILE

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


//! less than comparison
template <typename T>
class less_than_comp {
public:
	unsigned int N;
	bool operator() (T i,T j) { return i[N]<j[N]; }
};

int main(int argc, char **argv){


    UTILS::options_manager opt_mgr;
    opt_mgr.register_option("times",string(""),"times file");
    opt_mgr.register_option("params",string(""),"single parameters");
    opt_mgr.register_option("cond",string(""),"condition definition");

    opt_mgr.register_option("cond_species_only",bool(false),"output only the species defined in the condition file");
    opt_mgr.register_option("species",string(""),"species name to output");

    opt_mgr.register_option("para1",string(""),"parameter name to vary");
    opt_mgr.register_option("para2",string(""),"parameter name to vary");

    opt_mgr.register_option("para1:start",double(0),"parameter start range");
    opt_mgr.register_option("para2:start",double(0),"parameter start range");

    opt_mgr.register_option("para1:end",double(1),"parameter end range");
    opt_mgr.register_option("para2:end",double(1),"parameter end range");

    opt_mgr.register_option("para1:steps",(unsigned int)(100),"parameter range steps");
    opt_mgr.register_option("para2:steps",(unsigned int)(100),"parameter range steps");


    opt_mgr.register_option("output_min_max",bool(false),"output min/max values of trajectory instead of time points");

    opt_mgr.register_option("help",bool(false),"output help");

    std::vector<string> compulsory_flags;
    compulsory_flags.push_back("times");
    compulsory_flags.push_back("params");
    compulsory_flags.push_back("para1");
    compulsory_flags.push_back("para2");
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
    	cerr << "ERROR: can't use --cond_species_only option without also defining --cond\n";
    	return -1;
    }

    if ( !(opt_mgr.is_set("cond") || opt_mgr.is_set("species")) ){
    	cerr << "ERROR: at least one of --cond or --species needs to be set\n";
    	return -1;
    }

    const bool output_cond_species = opt_mgr.get_option_value<bool>("cond_species_only");
    const bool using_species = opt_mgr.is_set("species");

    if (output_cond_species && using_species){
    	cerr << "ERROR: can't use --cond_species_only option simultaneously with the --species option" << endl;
    	return -1;
    }

    const bool output_min_max = opt_mgr.get_option_value<bool>("output_min_max");

    //const double endtime = opt_mgr.get_option_value<double>("end_time");
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

    joint_state = sys.get_joint_para_and_state();
    const ode_model_functor::joint_para_x_y_state_type start_joint_state = joint_state;

    if (using_species){
    	string species_name = opt_mgr.get_option_value<string>("species");
    	const unsigned int sp_index = ode_model_functor::get_joint_x_y_state_index(species_name);
    	if (sp_index != NAME_NOT_FOUND){
    		species = sp_index;
    		scale = 1;
    	}
    	else {
    		cerr << "ERROR: species label " << species_name << " in --output_species_list not found" << endl;
    		return -1;
    	}
    }

    cout << "# outputting species: " <<  ode_model_functor::get_joint_state_name_vec()[species] << std::endl;


    //get_joint_para_x_y_state_index()
    const string para1_name = opt_mgr.get_option_value<string>("para1");
    const string para2_name = opt_mgr.get_option_value<string>("para2");

    const unsigned int p1 = ode_model_functor::get_joint_para_x_y_state_index(para1_name);
    if (p1 == NAME_NOT_FOUND){
    	cerr << "ERROR: species label " << para1_name << " in --para1 not found" << endl;
    	return -1;
    }
    const unsigned int p2 = ode_model_functor::get_joint_para_x_y_state_index(para2_name);
    if (p2 == NAME_NOT_FOUND){
    	cerr << "ERROR: species label " << para2_name << " in --para2 not found" << endl;
    	return -1;
    }

    const double para1_start = opt_mgr.get_option_value<double>("para1:start");
    const double para2_start = opt_mgr.get_option_value<double>("para2:start");
    const double para1_end = opt_mgr.get_option_value<double>("para1:end");
    const double para2_end = opt_mgr.get_option_value<double>("para2:end");
    const unsigned int para1_steps = opt_mgr.get_option_value<unsigned int>("para1:steps");
    const unsigned int para2_steps = opt_mgr.get_option_value<unsigned int>("para2:steps");
    const double para1_stepsize = (para1_end - para1_start) / (double(para1_steps - 1));
    const double para2_stepsize = (para2_end - para2_start) / (double(para2_steps - 1));
    if (para1_stepsize <= 0 || para1_stepsize <= 0){
    	cerr << "ERROR: check the ranges for para1 and/or para2" << endl;
    	return -1;
    }


    less_than_comp<typename ode_model_functor::joint_x_y_state_type> comp_obj;
    comp_obj.N = species;

    if (output_min_max){
    	cout << "#" << para1_name << "\t" << para2_name << "\tmin\tmax" << std::endl;
    }
    else {
    	cout << "#" << para1_name << "\t" << para2_name << "\t0\t" << times.back() << std::endl;
    }

    for (unsigned int ii = 0; ii < para1_steps; ii++){
    	const double para1_val = para1_start + (double(ii)*para1_stepsize);
    	for (unsigned int jj = 0; jj < para1_steps; jj++){
    		joint_state = start_joint_state;
    		const double para2_val = para2_start + (double(jj)*para2_stepsize);
    		//std::cout	<<  para1_val << "\t" << para2_val << "\t" << p1 << "\t" << p2 << endl;
    		//std::cout	<< joint_state.size() << endl;
    		joint_state[p1] = para1_val;
    		joint_state[p2] = para2_val;
    		sys.set_joint_para_and_state(joint_state);
    		//cout << "debug1" << endl;
    		vector<typename ode_model_functor::joint_x_y_state_type> traj = ODE::simple_run_trajectory(sys,  times);
    		if (output_min_max){
    			const double max_val = (*std::max_element(traj.begin(), traj.end(), comp_obj))[species];
    			const double min_val = (*std::min_element(traj.begin(), traj.end(), comp_obj))[species];
    			std::cout	<<  para1_val << "\t" << para2_val << "\t"
    			    						<< min_val << "\t" <<  max_val << std::endl;

    		}
    		else {
    			std::cout	<<  para1_val << "\t" << para2_val << "\t"
    						<< traj.front()[species] << "\t" <<  traj.back()[species] << std::endl;
    		}
    	}
    }
}
