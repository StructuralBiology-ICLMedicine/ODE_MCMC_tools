#ifndef ODE_LIKELIHOOD_FUNCTION_H_INCLUDED
#define ODE_LIKELIHOOD_FUNCTION_H_INCLUDED


#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdio>
#include <sstream>
#include <exception>
#include <utility>
#include <deque>
#include <functional>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include <boost/math/distributions/normal.hpp>
#include  <cmath>

#include "MersenneTwister.h"

//#include "hist_utils.h"
//#include "hist2d_interp.h"

#include "multivariate_normal.h"
#include "multivariate_t_dist_pdf.h"

#include "eigen_matrix_utils.h"

#include "common_model_defs.h"

#include "model_numerical_jacobian.h"

#include "mcmc_interfaces.h"

namespace ODE {


typedef std::map< unsigned int, double > uint_double_map;
typedef std::map< unsigned int, MultivariateNormalPDF > MultivariateNormalPDF_map;
typedef std::map< unsigned int, MultivariateStudentsTPDF > MultivariateStudentsTPDF_map;
typedef std::vector<double> double_vec;
typedef std::vector<unsigned int> uint_vec;

typedef std::map< unsigned int, double_vec > uint_double_vec_map;



	template<class state_type_>
	class observer_whole_traj_functor {
	public:
		std::vector<double>& store_times;
		std::vector<state_type_>& store_traj;

	public:
		observer_whole_traj_functor(std::vector<double>& store_times_, std::vector<state_type_>& store_traj_) : store_times(store_times_), store_traj(store_traj_){
		}

		void operator()(const state_type_ &x , const double t){
			store_traj.push_back(x);
			store_times.push_back(t);
		}

	};

// TODO change this to store multiple species trajectories
class sim_condition{

public:
	// note all conc units in uM

	// uint double map mapping from joint_para_ct_state_type index to parameter value
	uint_double_map cond_vals;



    // This maps species index to a MVN PDF object
	MultivariateNormalPDF_map mvn_pdf_map;
	MultivariateStudentsTPDF_map mvn_t_pdf_map;

	// uint double map mapping from species index to weight (default to 1)
	uint_double_map weights;

	uint_double_vec_map times_map;
	double_vec merged_times;


	sim_condition() {
		merged_times.push_back(0);
		//exp_traj = double_vec_map();
		//exp_traj_stddev = double_vec_map();
	}

	void add_to_merged_times(double_vec &vec){
		double_vec temp;
		std::set_union(merged_times.begin(), merged_times.end(),
				vec.begin(), vec.end(),
				std::back_inserter(temp));
		merged_times = temp;
		/*
		std::cout << "# merged_times:";
		for (auto it = merged_times.begin(); it != merged_times.end(); it++){
			std::cout << "\t" << *it;
		}
		std::cout << std::endl;
		*/
	}

	void merge_times(){
		for (uint_double_vec_map::iterator it = times_map.begin(); it != times_map.end(); it++ ){
			/*
			double_vec temp;
			std::set_union(merged_times.begin(), merged_times.end(),
					it->second.begin(), it->second.end(),
					std::back_inserter(temp));
			merged_times = temp;
			*/
			add_to_merged_times(it->second);
		}
		std::cout << "# merged_times:";
		for (auto it = merged_times.begin(); it != merged_times.end(); it++){
			std::cout << "\t" << *it;
		}
		std::cout << std::endl;
	}



	void print_summary() const{
	    /*
		cout << "xylr_conc:" << xylr_conc << "\t"
		     << "xylose_conc:" << xylose_conc << "\t"
		     << "dna_conc:" << dna_conc << "\t"
		     << "exp_traj.size():" << exp_traj.size() << "\t"
		     << "exp_traj_stddev.size():" << exp_traj_stddev.size() << "\t";
        */
	}

	static std::string make_hash(const uint_double_map condition_vals) {
        uint_double_map::const_iterator it;
        std::stringstream shash(std::stringstream::in | std::stringstream::out);
        shash << "cond";
        for (it = condition_vals.begin(); it != condition_vals.end(); it++){
            shash << "|" << it->first << "->" << it->second;
        }
        return shash.str();
    }

    std::string get_hash() const {
        return make_hash(cond_vals);
    }

};

typedef std::vector<sim_condition> sim_condition_vec;
typedef std::map<std::string, sim_condition> sim_condition_map;

// mvn log likelihood functor
    template<typename system_>
    class ode_loglikelihood_mvn_functor {
    private:


    public:

        double_vec times;
        sim_condition_map conds;




        std::vector<typename system_::joint_x_y_state_type> run_trajectory(system_ &sys, double_vec times)
        {



        	sys->initial_assignments();
            std::pair< system_,  model_jacobian<system_>  > sys_pair = std::make_pair(sys, model_jacobian<system_>(sys));
            //model_jacobian jaco(sys);
            //std::pair< std::reference_wrapper<system_>,  std::reference_wrapper<model_jacobian<system_> > > sys_pair(std::reference_wrapper<system_>(sys),
             //                                                                                                                         std::reference_wrapper<model_jacobian<system_> > (model_jacobian<system_>(sys)));

            //typename system_::state_type x = { 0.0 , 0.0 , sys.p[9], 0.0 }; // initial conditions: set mRNA capacity of  to initial mRNA capacity, otherwise all zero

            //typedef runge_kutta_cash_karp54< typename system_::state_type > error_stepper_type;
            //typedef runge_kutta_fehlberg78< typename system_::state_type > error_stepper_type;     //typedef runge_kutta_cash_karp54< typename system_::state_type > error_stepper_type;
            //typedef runge_kutta_dopri5< typename system_::state_type > error_stepper_type;
			typedef boost::numeric::odeint::rosenbrock4< double > error_stepper_type;
            //typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
            //controlled_stepper_type controlled_stepper;
            //auto controlled_stepper = make_controlled( 1.0e-6 , 1.0e-6 , error_stepper_type() );;
            auto controlled_stepper = boost::numeric::odeint::rosenbrock4_controller<error_stepper_type>();

            double_vec store_times;
            std::vector< typename system_::state_type >  store_traj;

            observer_whole_traj_functor< typename system_::state_type > obs = observer_whole_traj_functor< typename system_::state_type >(store_times, store_traj);
            //system_ sys = system_(p);
            integrate_times( controlled_stepper , sys_pair, sys_pair.first.x , times.begin() , times.end() , 0.1 , obs );
            //integrate_times( controlled_stepper , sys_pair, sys.x , times.begin() , times.end() , 0.1 , obs );

            std::vector< typename system_::joint_x_y_state_type >  joint_store_traj;
            joint_store_traj.reserve(store_traj.size());
            for (typename std::vector< typename system_::state_type >::const_iterator it = store_traj.begin(); it != store_traj.end(); it++)
            {
                typename system_::state_y_type y_vec = sys.get_y_vector(*it);
                typename system_::joint_x_y_state_type joint_state = system_::new_joint_x_y_state_type();
                std::copy(it->begin(), it->end(), joint_state.begin());
                std::copy(y_vec.begin(), y_vec.end(), joint_state.begin() + it->size());
                joint_store_traj.push_back(joint_state);
            }



            if (store_times[0] != 0)
            {
                std::cerr << "ERROR !!!!!!!!!!!!!!!!!!!! I think there's an odeint bug when the time vector does not start at t=0. " << store_times[0] << "\t" << times[0] << std::endl ;
            }

            return joint_store_traj;
        }





        double get_LL(std::vector<typename system_::joint_x_y_state_type>& sim_traj, MultivariateNormalPDF_map& mvnpdf_map)
        {

            //const double min_abs_stddev = 0.2;

            double sum_ll = 0;

            //uint_vec::const_iterator iter;
            for (MultivariateNormalPDF_map::const_iterator iter = mvnpdf_map.begin(); iter != mvnpdf_map.end(); iter++)
            {
                const unsigned int this_species = iter->first;

                /*
                if (sim_traj.size() != exp_traj[this_species].size()){
                        cerr << "ERROR" << endl;
                        return -9999;
                }
                // TODO this seems to be same thing as previous if so what is the point of it?
                else if (exp_traj[this_species].size() != sim_traj.size()){
                        cerr << "ERROR" << endl;
                        return -9999;
                }
                */

                //cout << sim_traj.size() << endl;



                Eigen::Matrix<double,Eigen::Dynamic,1> sim_traj_mat(sim_traj.size());
                for (unsigned long i = 0; i < sim_traj.size(); i++)
                {
                    sim_traj_mat(i) = sim_traj[i][this_species];
                }
                const double this_ll = mvnpdf_map[this_species].getLogPdf(sim_traj_mat); //sim_traj_mat
                //cout << "DETERMINANT: " << mvnpdf_map[this_species].covar.determinant() << endl;
                //cout << "sim_traj " << endl;
                //cout << sim_traj_mat << endl;
                //cout << "means " << endl;
                //cout << mvnpdf_map[this_species].mean << endl;
                //cout << "covar " << endl;
                //cout << mvnpdf_map[this_species].covar_inv << endl;
                //cout << "this_ll: "<< this_ll << endl;
                sum_ll += this_ll;


                /*
                for (unsigned long i = 0; i < sim_traj.size(); i++){
                        normal norm = normal(exp_traj[this_species][i],
                                             std::max(min_abs_stddev[iter],std::max(exp_traj_stddev[this_species][i],
                                            (min_stddev_percent[iter] * exp_traj[this_species][i]))));
                        double val = pdf(norm, sim_traj[i][this_species]);
                        sum_ll += log(val);
                }
                */
            }

            return sum_ll;
        }




        double get_loglikelihood(system_ sys, double_vec times, MultivariateNormalPDF_map& mvnpdf_map )
        {

            std::vector<typename system_::joint_x_y_state_type> store_traj = this->run_trajectory(sys,  times);

            return this->get_LL(store_traj,   mvnpdf_map);

        }





        double operator()(system_ sys)
        {
            sim_condition_map::iterator iter;
            double ll_sum = 0;
            //typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };
            const typename system_::joint_para_x_y_state_type  common_paras = sys.get_joint_para_and_state();
            for (iter = conds.begin(); iter != conds.end(); iter++)
            {

                sys.set_joint_para_and_state( common_paras);
                apply_condition(sys, iter->second);

                ll_sum += this->get_loglikelihood(sys, times, iter->second.mvn_pdf_map);
            }
            return ll_sum;
        }



    };


bool readMeans(std::string filename, double_vec& times, double_vec& means){
    times.clear();
    means.clear();
    std::ifstream input(filename.c_str(), std::ios::in);
    string_vector SplitVec;
    std::string lineStr;
    unsigned long length, lineNum = 0;

    if (!input.is_open()){
        std::cerr << "ERROR: readMeans(): couldn't read file: " << filename << std::endl;
        throw std::exception();
    }

    while ( !input.eof() )
    {
        getline(input, lineStr);
        std::string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 1)
        {
            boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
            const double time = boost::lexical_cast<double>(SplitVec[0]);
            const double val = boost::lexical_cast<double>(SplitVec[1]);
            times.push_back(time);
            means.push_back(val);
        }
    }

    return true;
}

double_vec readTimes(std::string filename){
    double_vec times;
    std::ifstream input(filename.c_str(), std::ios::in);
    string_vector SplitVec;
    std::string lineStr;
    unsigned long length, lineNum = 0;

    while ( !input.eof() )
    {
        getline(input, lineStr);
        std::string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 0)
        {
            boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
            const double time = boost::lexical_cast<double>(SplitVec[0]);
            times.push_back(time);
        }
    }

    return times;
}

template <typename system_>
sim_condition parse_condition_spec_file(const std::string filename, unsigned int& species, double &scale,
		double &min_variance, double &weight){
    sim_condition cond;


    std::ifstream input(filename.c_str(), std::ios::in);
    string_vector SplitVec;
    std::string lineStr;
    unsigned long length, lineNum = 0;

    species = NAME_NOT_FOUND;
    weight = 1;

    if (!input.is_open()){
        std::cerr << "ERROR: parse_condition_spec_file(): couldn't open file: " << filename  << std::endl;
        throw std::exception();
    }

    while ( !input.eof() )
    {
        getline(input, lineStr);
        std::string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 0)
        {
            //std::cout << lineStr << endl;
            //std::cout << lineNum << endl;
            boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
            if (SplitVec.size() >= 2 && lineStr[0] != '#'){
                const std::string type = SplitVec[0];
                if (type.compare("SPECIES") == 0){
                    //std::cout << "# SPECIES " << endl;
                    species = system_::get_joint_x_y_state_index(SplitVec[1]);
                    if (species == NAME_NOT_FOUND){
                        std::cerr << "ERROR: parse_condition_spec_file(): couldn't find species name: " << SplitVec[1] << std::endl;
                        throw std::exception();
                    }
                }
                else if (type.compare("UNIT_SCALE") == 0){
                    //std::cout << "# UNIT_SCALE " << endl;
                    scale = boost::lexical_cast<double>(SplitVec[1]);
                }
                else if (type.compare("CONDITION") == 0){
                    //std::cout << "# CONDITION " << endl;
                    const unsigned int this_id = system_::get_joint_para_x_y_state_index(SplitVec[1]);
                    const double value = boost::lexical_cast<double>(SplitVec[2]);
                    cond.cond_vals[this_id] = value;
                }
                else if (type.compare("MIN_STDDEV") == 0){
                    //std::cout << "# MIN_STDDEV " << endl;
                    min_variance = scale*scale*boost::lexical_cast<double>(SplitVec[1])*boost::lexical_cast<double>(SplitVec[1]);
                }
                else if (type.compare("WEIGHT") == 0){
                	weight = boost::lexical_cast<double>(SplitVec[1]);
                }
            }
        }
    }
    return cond;
}



template <typename system_>
bool parse_MVN_list_file(std::string filename, sim_condition_map& conds, double_vec& times, double t_deg_free = 3)
{
    //const bool use_new_key = true;
    //times.clear();
    std::ifstream input(filename.c_str(), std::ios::in);
    string_vector SplitVec;
    std::string lineStr;
    unsigned long length, lineNum = 0;

    double_vec combined_times;

    unsigned int species = NAME_NOT_FOUND;
    double scale = 1.0;
    double min_variance = 0;
    double weight = 1;

    if (!input.is_open()){

    }

    while ( !input.eof() )
    {
        getline(input, lineStr);
        std::string resStr;
        lineNum++;

        length = lineStr.length();


        if (length > 0)
        {
            boost::split( SplitVec, lineStr, boost::is_any_of("\t") );


            if (SplitVec.size() >= 2 && lineStr[0] != '#'){
                // parsing trajectories
                //string old_key = "d" + SplitVec[0] + "_x" + SplitVec[1] + "_r" + SplitVec[2];
                //double dna_conc = lexical_cast<double>(SplitVec[0])  /1000.0;
                //double xylose_conc = lexical_cast<double>(SplitVec[1]) * 1000.0;
                //double xylr_conc = lexical_cast<double>(SplitVec[2]);
                std::string cond_spec_filename = SplitVec[0];
                std::string sigma_filename = SplitVec[1];
                std::string means_filename = SplitVec[2];

                std::cout << "# reading condition specs file: " << cond_spec_filename << std::endl;
                sim_condition this_cond = parse_condition_spec_file<system_>(cond_spec_filename, species, scale, min_variance, weight);
                const std::string key = this_cond.get_hash();//make_new_key(dna_conc, xylose_conc, xylr_conc);
                //std:stringstream skey(std::stringstream::in | std::stringstream::out);
                //skey << "d"  << setfill('0') << setw(5) << (dna_conc*1000000) << "_x" << setfill('0') << setw(5) << (xylose_conc) << "_r" << setfill('0') << setw(5) << (xylr_conc*1000);
                //string key = use_new_key ? new_key : old_key;// skey.str();//

                if (conds.find(key) == conds.end())
                {
                    // create new element
                    std::cout << "# creating new condition: " << key  << std::endl;
                    conds[key] = this_cond;//sim_condition(xylr_conc, xylose_conc, dna_conc);
                }

                std::cout << "# reading LL Sigma: " << sigma_filename << std::endl;
                Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> cov_mat = (scale*scale) * readDynamicMatrix(sigma_filename);
                for (int i = 0; i < cov_mat.rows(); i++){
                    if (cov_mat(i,i) < min_variance){
                        cov_mat(i,i) = min_variance;
                    }
                }
                double_vec means_vec;
                std::cout << "# reading LL Means: " << means_filename << std::endl;
                readMeans(means_filename, times, means_vec);
                Eigen::Matrix<double,Eigen::Dynamic,1> means_mat(means_vec.size(), 1);
                for (unsigned int i = 0; i< means_vec.size(); i++){
                    means_mat(i) = scale * means_vec[i];
                }

                if (means_mat.rows() != cov_mat.rows()){
                    std::cerr << "ERROR: Dimensions of means and covs do not match: " << means_mat.rows() << "\t" << cov_mat.rows() << std::endl;
                    return false;
                }

                conds[key].mvn_pdf_map[species].setMean(means_mat);
                conds[key].mvn_pdf_map[species].setCovar(cov_mat);

                conds[key].mvn_t_pdf_map[species].setMean(means_mat);
                conds[key].mvn_t_pdf_map[species].setCovar(cov_mat, t_deg_free);

                conds[key].weights[species] = weight;

                std::cout << "# weight: " << weight << std::endl;

                conds[key].times_map[species] = times;
                conds[key].add_to_merged_times(times);

                {
                	double_vec temp;
                	std::set_union(conds[key].merged_times.begin(), conds[key].merged_times.end(),
                			combined_times.begin(), combined_times.end(),
							std::back_inserter(temp));
                	combined_times = temp;
                }

            }
        }
    }
    times = combined_times;

    std::cout << "# combined_times:";
    for (auto it = combined_times.begin(); it != combined_times.end(); it++){
    	std::cout << "\t" << *it;
    }
    std::cout << std::endl;

    std::cout << "# Number of different conditions:\t" << conds.size() << std::endl;

    return true;
}






template <typename system_>
void apply_condition(system_ &sys, sim_condition &cond)
{
    typename system_::joint_para_x_y_state_type  paras = sys.get_joint_para_and_state();
    uint_double_map::const_iterator cond_val_it;
    for (cond_val_it = cond.cond_vals.begin(); cond_val_it != cond.cond_vals.end(); cond_val_it++)
    {
        paras[cond_val_it->first] = cond_val_it->second;
    }
    sys.set_joint_para(paras);
}

template <typename system_ptr_>
void apply_condition_ptr(system_ptr_ &sys_ptr, sim_condition &cond)
{
    typename system_ptr_::element_type::vector_type  paras = sys_ptr->get_joint_para_and_state();
    uint_double_map::const_iterator cond_val_it;
    for (cond_val_it = cond.cond_vals.begin(); cond_val_it != cond.cond_vals.end(); cond_val_it++)
    {
        paras[cond_val_it->first] = cond_val_it->second;
    }
    sys_ptr->set_joint_para_and_state(paras);
}



	// uninformative prior function
	template <typename para_type_>
	class uninf_prior_functor : public MCMC::prior_functor_interface<para_type_> {
	public:
		double operator()(const para_type_&  p) override{
			return 0.0;
		}
		
		virtual std::shared_ptr< MCMC::prior_functor_interface<para_type_> > clone() override{
			std::shared_ptr< MCMC::prior_functor_interface<para_type_> > new_ptr(new uninf_prior_functor<para_type_>(*this));
			return new_ptr;
		}
	};



	// multivariate gaussian prior functor
	//template <typename para_type_, unsigned int DIM>
	template<typename system_>
	class multigauss_prior_functor : public MCMC::prior_functor_interface<typename system_::vector_type> {
	private:
		uint_vec mapping;
		Eigen::Matrix<double,Eigen::Dynamic,1> mapped_vec;
		MultivariateNormalPDF pdf;

	public:
		
		virtual std::shared_ptr< MCMC::prior_functor_interface<typename system_::vector_type> > clone() override{
			std::shared_ptr< MCMC::prior_functor_interface<typename system_::vector_type> > new_ptr(new multigauss_prior_functor<system_>(*this));
			return new_ptr;
		}


		void set_mapping(uint_vec mapping_){
			mapping = mapping_;
		}

		bool load_mapping(const system_& sys,const std::string filename)
		{
			mapping.clear();
			std::ifstream input(filename.c_str(), std::ios::in);

			string_vector SplitVec;
			std::string lineStr;
			unsigned long length, lineNum = 0;

			while ( !input.eof() )
			{
				getline(input, lineStr);

				length = lineStr.length();

				if (length > 0)
				{
					boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
					if (SplitVec.size() >= 1 && lineStr[0] != '#')
					{
						const std::string name = SplitVec[0];

						const unsigned int joint_index = sys.get_joint_para_and_state_index(name);//system_::get_joint_para_x_y_state_index(name);

						std::cout << "# Reading prior mapping:\t" << name << "\t" << lineNum << "\t" << joint_index << std::endl;
						if (joint_index == NAME_NOT_FOUND){
							throw std::exception();
						}
						mapping.push_back(joint_index);
						lineNum++;
					}
				}
			}
			mapped_vec.resize(mapping.size());
			return true;
		}


		void set_covar(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_){
			pdf.setCovar(covar_);
		}

		void set_mean(Eigen::Matrix<double,Eigen::Dynamic,1> mean_){
			pdf.setMean(mean_);
		}

		bool verify() const{
			if (pdf.get_mean().rows() != pdf.get_covar().rows()){
				std::cerr << "ERROR: multigauss_prior_functor: prior mean and covar dimensions don't match" << std::endl;
				return false;
			}
			if (pdf.get_covar().rows() != mapping.size()){
				std::cerr << "ERROR: multigauss_prior_functor: prior covar and mapping dimensions don't match" << std::endl;
				return false;
			}
			return true;
		}

		double operator()(const typename system_::vector_type&  p) override{
			for (unsigned int ii = 0; ii < mapping.size(); ii++){
				mapped_vec(ii) = p[mapping[ii]];
			}
			const double logpdf = pdf.getLogPdf(mapped_vec);
			//cout << "# logPrior:\t" << logpdf << endl;
			//cout << "# params:\t" << mapped_vec.transpose() << endl;
			//cout << "# means:\t" << pdf.get_mean().transpose() << endl;
			return logpdf;
		}
	};


	template <typename para_type_>
	bool parse_initial_params(std::string filename, para_type_& p){
		const  int len = p.size();
		std::ifstream input(filename.c_str(), std::ios::in);
		string_vector SplitVec;
		std::string lineStr;
		unsigned long length, lineNum = -1;

		if (!input.is_open()){
			std::cerr << "ERROR: parse_initial_params(): couldn't open file: " << filename << std::endl;
			throw std::exception();
		}

		while ( !input.eof() ) {
			getline(input, lineStr);
			boost::trim(lineStr);
			std::string resStr;
			lineNum++;

			length = lineStr.length();


			if (length > 0) {
				boost::split( SplitVec, lineStr, boost::is_any_of("\t ") );
				if (lineNum < len){
					p[lineNum]= boost::lexical_cast<double>(SplitVec[0]);
				}
				else {
					std::cerr << "ERROR: parse_initial_params(): line number out of bounds: expected len=" << len
							<< " lineNum=" << lineNum
							<< std::endl;
					throw std::exception();
				}
				//cout << lineNum << "\t" << p[lineNum] << endl;
			}
		}
		return true;
	}








	template<typename system_>
    class ode_loglikelihood_mvn_ptr_functor : public MCMC::loglikelihood_functor_interface<system_> {
    private:


    public:

        double_vec times;
        sim_condition_map conds;
    	bool use_t_dist;

    	ode_loglikelihood_mvn_ptr_functor() : use_t_dist(false){

    	}

		virtual std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > clone() override{
			std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > new_ptr(new ode_loglikelihood_mvn_ptr_functor<system_>(*this));
			return new_ptr;
		}




        std::vector<typename system_::vector_type> run_trajectory(std::shared_ptr<  system_ > sys_ptr, double_vec times)
        {


        	sys_ptr->initial_assignments();
            MCMC::ode_model_system_wrapper<std::shared_ptr<  system_ > > sys(sys_ptr);
            model_jacobian_wrapper<std::shared_ptr<  system_ > > jacob(sys_ptr);

            auto sys_pair = std::make_pair(sys, jacob );
            //std::pair< system_,  model_jacobian<system_>  > sys_pair = std::make_pair(*sys_ptr, model_jacobian<system_>(*sys_ptr));
            //model_jacobian jaco(sys);
            //std::pair< std::reference_wrapper<system_>,  std::reference_wrapper<model_jacobian<system_> > > sys_pair(std::reference_wrapper<system_>(sys),
             //                                                                                                                         std::reference_wrapper<model_jacobian<system_> > (model_jacobian<system_>(sys)));

            //typename system_::state_type x = { 0.0 , 0.0 , sys.p[9], 0.0 }; // initial conditions: set mRNA capacity of  to initial mRNA capacity, otherwise all zero

            //typedef runge_kutta_cash_karp54< typename system_::state_type > error_stepper_type;
            //typedef runge_kutta_fehlberg78< typename system_::state_type > error_stepper_type;     //typedef runge_kutta_cash_karp54< typename system_::state_type > error_stepper_type;
            //typedef runge_kutta_dopri5< typename system_::state_type > error_stepper_type;
			typedef boost::numeric::odeint::rosenbrock4< double > error_stepper_type;
            //typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
            //controlled_stepper_type controlled_stepper;
            //auto controlled_stepper = make_controlled( 1.0e-6 , 1.0e-6 , error_stepper_type() );;
            auto controlled_stepper = boost::numeric::odeint::rosenbrock4_controller<error_stepper_type>();

            double_vec store_times;
            std::vector< typename system_::vector_type >  store_traj;

            observer_whole_traj_functor< typename system_::vector_type > obs = observer_whole_traj_functor< typename system_::vector_type >(store_times, store_traj);
            //system_ sys = system_(p);
            typename system_::vector_type x0 = sys_ptr->get_state();//sys_pair.first.x ;
            integrate_times( controlled_stepper , sys_pair, x0 , times.begin() , times.end() , 0.1 , obs );
            //integrate_times( controlled_stepper , sys_pair, sys.x , times.begin() , times.end() , 0.1 , obs );

            std::vector< typename system_::vector_type >  joint_store_traj;
            joint_store_traj.reserve(store_traj.size());
            for (typename std::vector< typename system_::vector_type >::const_iterator it = store_traj.begin(); it != store_traj.end(); it++)
            {
                typename system_::vector_type y_vec = sys_ptr->get_dependent_vars(*it);
                typename system_::vector_type joint_state = sys_ptr->new_joint_para_and_state_type();
                std::copy(it->begin(), it->end(), joint_state.begin());
                std::copy(y_vec.begin(), y_vec.end(), joint_state.begin() + it->size());
                joint_store_traj.push_back(joint_state);
            }



            if (store_times[0] != 0)
            {
                std::cerr << "ERROR !!!!!!!!!!!!!!!!!!!! I think there's an odeint bug when the time vector does not start at t=0. " << store_times[0] << "\t" << times[0] << std::endl ;
            }

            return joint_store_traj;
        }



        double get_LL(std::vector<typename system_::vector_type>& sim_traj, MultivariateNormalPDF_map& mvnpdf_map)
        {
            double sum_ll = 0;
            for (MultivariateNormalPDF_map::const_iterator iter = mvnpdf_map.begin(); iter != mvnpdf_map.end(); iter++)
            {
                const unsigned int this_species = iter->first;
                Eigen::Matrix<double,Eigen::Dynamic,1> sim_traj_mat(sim_traj.size());
                for (unsigned long i = 0; i < sim_traj.size(); i++)
                {
                    sim_traj_mat(i) = sim_traj[i][this_species];
                }
                const double this_ll = mvnpdf_map[this_species].getLogPdf(sim_traj_mat); //sim_traj_mat
                sum_ll += this_ll;
            }
            return sum_ll;
        }

        double get_LL(std::vector<typename system_::vector_type>& sim_traj, MultivariateStudentsTPDF_map& mvnpdf_map)
        {
        	double sum_ll = 0;
        	for (MultivariateStudentsTPDF_map::const_iterator iter = mvnpdf_map.begin(); iter != mvnpdf_map.end(); iter++)
        	{
        		const unsigned int this_species = iter->first;
        		Eigen::Matrix<double,Eigen::Dynamic,1> sim_traj_mat(sim_traj.size());
        		for (unsigned long i = 0; i < sim_traj.size(); i++)
        		{
        			sim_traj_mat(i) = sim_traj[i][this_species];
        		}
        		const double this_ll = mvnpdf_map[this_species].getLogPdf(sim_traj_mat); //sim_traj_mat
        		sum_ll += this_ll;
        	}
        	return sum_ll;
        }

        double get_loglikelihood(std::shared_ptr<  system_ > sys_ptr, double_vec times, MultivariateNormalPDF_map& mvnpdf_map )
        {
            std::vector<typename system_::vector_type> store_traj = this->run_trajectory(sys_ptr,  times);
            return this->get_LL(store_traj,   mvnpdf_map);
        }

        double get_loglikelihood(std::shared_ptr<  system_ > sys_ptr, double_vec times, MultivariateStudentsTPDF_map& mvnpdf_map )
        {
        	std::vector<typename system_::vector_type> store_traj = this->run_trajectory(sys_ptr,  times);
        	return this->get_LL(store_traj,   mvnpdf_map);
        }


        double get_LL(std::vector<typename system_::vector_type>& sim_traj, sim_condition& cond)
        {
        	double sum_ll = 0;
        	for (uint_double_vec_map::const_iterator iter = cond.times_map.begin(); iter != cond.times_map.end(); iter++)
        	{
        		const unsigned int this_species = iter->first;
        		const double_vec& times = iter->second;
        		Eigen::Matrix<double,Eigen::Dynamic,1> sim_traj_mat(times.size());//sim_traj.size());
        		unsigned int pos = 0;
        		for (unsigned long i = 0; i < times.size(); i++)
        		{
        			for (unsigned long j = pos; j < cond.merged_times.size(); j++){
        				if (times[i] == cond.merged_times[j]){
        					sim_traj_mat(i) = sim_traj[j][this_species];
        					pos++;
        					break;
        				}
        			}
        		}
        		const double this_weight = cond.weights[this_species];
        		const double this_ll = use_t_dist ? this_weight*cond.mvn_t_pdf_map[this_species].getLogPdf(sim_traj_mat) : this_weight*cond.mvn_pdf_map[this_species].getLogPdf(sim_traj_mat); //sim_traj_mat
        		sum_ll += this_ll;
        	}
        	return sum_ll;
        }


        double get_loglikelihood(std::shared_ptr<  system_ > sys_ptr, sim_condition& cond){
        	std::vector<typename system_::vector_type> store_traj = this->run_trajectory(sys_ptr,  cond.merged_times);
        	return this->get_LL(store_traj, cond);
        }




        virtual double operator()( std::shared_ptr<  system_ > sys_ptr ) override{
            sim_condition_map::iterator iter;
            double ll_sum = 0;
            const typename system_::vector_type  common_paras = sys_ptr->get_joint_para_and_state();

            for (iter = conds.begin(); iter != conds.end(); iter++)
            {

                sys_ptr->set_joint_para_and_state( common_paras);
                apply_condition_ptr(sys_ptr, iter->second);

                ll_sum += get_loglikelihood(sys_ptr, iter->second);

                /*
                if (!use_t_dist){
                	ll_sum += this->get_loglikelihood(sys_ptr, times, iter->second.mvn_pdf_map);
                }
                else {
                	ll_sum += this->get_loglikelihood(sys_ptr, times, iter->second.mvn_t_pdf_map);
                }
                */
            }

            sys_ptr->set_joint_para_and_state( common_paras);
            return ll_sum;
        }

    };



	template <typename system_>
	std::vector<typename system_::joint_x_y_state_type> simple_run_trajectory(system_ sys, double_vec times){



		sys.initial_assignments();
		auto sys_pair = std::make_pair(sys, model_jacobian<system_>(sys));


		typedef boost::numeric::odeint::rosenbrock4< double > error_stepper_type;

		auto controlled_stepper = boost::numeric::odeint::rosenbrock4_controller<error_stepper_type>();

		double_vec store_times;
		std::vector< typename system_::state_type >  store_traj;

		observer_whole_traj_functor< typename system_::state_type > obs = observer_whole_traj_functor< typename system_::state_type >(store_times, store_traj);
		//system_ sys = system_(p);
		integrate_times( controlled_stepper , sys_pair, sys_pair.first.x , times.begin() , times.end() , 0.1 , obs );

		std::vector< typename system_::joint_x_y_state_type >  joint_store_traj;
		joint_store_traj.reserve(store_traj.size());
		for (typename std::vector< typename system_::state_type >::const_iterator it = store_traj.begin(); it != store_traj.end(); it++){
			typename system_::state_y_type y_vec = sys.get_y_vector(*it);
			typename system_::joint_x_y_state_type joint_state = system_::new_joint_x_y_state_type();
			std::copy(it->begin(), it->end(), joint_state.begin());
			std::copy(y_vec.begin(), y_vec.end(), joint_state.begin() + it->size());
			joint_store_traj.push_back(joint_state);
		}



		if (store_times[0] != 0){
			std::cerr << "ERROR !!!!!!!!!!!!!!!!!!!! I think there's an odeint bug when the time vector does not start at t=0. " << store_times[0] << "\t" << times[0] << std::endl ;
		}

		return joint_store_traj;
	}

	template <typename system_>
	std::vector<typename system_::vector_type> parse_parameters_list(std::string filename, const unsigned int start_pos = 7)
	{

		//const unsigned int start_pos = 7;
		std::ifstream input(filename.c_str(), std::ios::in);
		string_vector SplitVec;
		std::string lineStr;
		unsigned long length, lineNum = -1;

		std::vector<typename system_::joint_x_y_state_type> rtn_vec;

		while ( !input.eof() )
		{
			getline(input, lineStr);
			boost::trim(lineStr);
			std::string resStr;
			lineNum++;

			length = lineStr.length();

			//para_type_ p;
			if (length > 0)
			{
				boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
				//cout << lineStr << endl;
				if ( SplitVec[0].substr(0,1).compare("#") != 0 )
				{
					typename system_::vector_type p = system_().new_joint_para_and_state_type();//system_::new_joint_x_y_state_type();
					for (int i = start_pos; i < (start_pos + p.size()); i++)
					{
						if (i < SplitVec.size()){
							//cout << i << "\t" << SplitVec[i] << endl;
							try {
								double val = boost::lexical_cast<double>(SplitVec[i]);
								p[i-start_pos]= val;
							}
							catch (std::exception e) {
								std::cerr << "ERROR: problem parsing line: " << lineNum
										<< " position: " << i
										<< " value '" <<  SplitVec[i] << "'"
										<< std::endl;
								throw e;
							}
						}
						else {
							std::cerr << "ERROR: not enough parameters parsing line: " << lineNum
									<< " position: " << i
									<< std::endl;
							throw std::exception();
						}
					}
					rtn_vec.push_back(p);
				}
				//cout << lineNum << "\t" << p[lineNum] << endl;
			}
		}
		return rtn_vec;
	}



}// namespace ODE
#endif // ODE_LIKELIHOOD_FUNCTION_H_INCLUDED
