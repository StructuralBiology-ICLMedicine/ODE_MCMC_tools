//
//  mcmc_algorithms.h
//  ode_mcmc_code
//
//  Created by jmacdon on 11/09/2016.
//  Copyright © 2016 James T. MacDonald. All rights reserved.
//

#ifndef mcmc_ptr_mod_sel_algorithms_h
#define mcmc_ptr_mod_sel_algorithms_h




#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdio>
#include <sstream>
#include <exception>
#include <utility>
#include <deque>
#include <memory>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

#include <boost/math/distributions/normal.hpp>
#include  <cmath>

#include "MersenneTwister.h"


#include "multivariate_normal.h"

#include "eigen_matrix_utils.h"

#include "common_model_defs.h"

#include "model_numerical_jacobian.h"

#include "common_variables.h"

#include "mcmc_interfaces.h"

#include <functional>






namespace MCMC{
namespace PTR{
namespace MODEL_SEL{


    typedef std::shared_ptr<std::ofstream> ofstream_shared_ptr;
    typedef std::shared_ptr<MTRand> MTRand_shared_ptr;
	
	typedef std::vector<double> double_vec;
	typedef std::vector<unsigned long> ulong_vec;
	typedef std::vector<unsigned int> uint_vec;


    template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
    class mcmc_interface_v2{

    public:
		
		typedef std::vector<system_ptr_> system_ptr_vector;
		typedef std::vector<ll_functor_ptr_> ll_functor_ptr_vector;
		typedef std::vector<prior_ptr_> prior_ptr_vector;
		typedef std::vector<proposal_ptr_> proposal_ptr_vector;
		

    public:

        MTRand_shared_ptr rand_gen;
		
		std::ostream* stdout_;
		std::ostream* stderr_;
		
		uint_vec model_id_mapping;


		


		system_ptr_vector sys_vec;
        ll_functor_ptr_vector get_LL_vec;
        prior_ptr_vector log_prior_vec;
        proposal_ptr_vector mover_vec;

        unsigned int num_models;
        unsigned int current_model, proposed_model, last_model;


        std::string output_root;
        std::vector<typename system_ptr_::element_type::vector_type> acc_params;
        std::vector<typename system_ptr_::element_type::vector_type> prop_params;
		std::vector<typename system_ptr_::element_type::vector_type> best_params;



        unsigned long steps;
        unsigned long burn_in;
        unsigned long output_freq;
        double inv_temp;

        double_vec accepted_ll, accepted_prior, accepted_sum;
        double_vec best_accepted_ll_prior;

        double_vec this_ll, this_prior, this_sum;

        unsigned long accepted_count;
		ulong_vec current_step;
		ulong_vec covar_calc_step;

		unsigned long total_steps;

		unsigned long output_acc_rate_freq;
		unsigned long output_acc_rate_window;

		double_vec acc_rate_total, acc_rate_rolling;
		std::vector<std::deque<unsigned int> > accepted_history;
		ulong_vec accepted_history_count;

		double acceptance_probability;

		std::vector<ofstream_shared_ptr> output_chain;

		int move_type;

		double log_proposal_fwd, log_proposal_rev;

		double jump_move_probability; // probability of making an jump to different model
		
		bool setup_ok;
		
		bool is_running;

		boost::mutex is_running_mutex;
		
		bool supress_foutput;
		bool supress_mcmc_foutput;
		bool supress_prob_foutput;
		bool supress_chkpt_foutput;

		unsigned long steps_since_jump, min_steps_before_jump, min_jump_burnin;

		Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> jump_prop_count;
		Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> jump_acc_count;

		bool do_indep_jumps;

    public:
        mcmc_interface_v2(MTRand_shared_ptr rand_ptr, std::ostream& stdout__, std::ostream& stderr__, uint_vec model_id_mapping_,
						  system_ptr_vector sys_vec_, ll_functor_ptr_vector get_LL_vec_, prior_ptr_vector log_prior_vec_, proposal_ptr_vector mover_vec_, std::string output_root_,
                       unsigned long steps_, unsigned long burn_in_, unsigned long output_freq_, const double inv_temp_) :
                           rand_gen(rand_ptr),stdout_(&stdout__), stderr_(&stderr__), model_id_mapping(model_id_mapping_), sys_vec(sys_vec_), get_LL_vec(get_LL_vec_), log_prior_vec(log_prior_vec_), mover_vec(mover_vec_),
                           num_models(sys_vec_.size()), current_model(0), proposed_model(0), last_model(0),
                           output_root(output_root_), acc_params(num_models), prop_params(num_models), best_params(num_models),
                           steps(steps_), burn_in(burn_in_), output_freq(output_freq_), inv_temp(inv_temp_),
                           accepted_ll(num_models,0), accepted_prior(num_models,0), accepted_sum(num_models,0), best_accepted_ll_prior(num_models,0),
                           this_ll(num_models,0), this_prior(num_models,0), this_sum(num_models,0),
                           accepted_count(0), current_step(num_models,1), covar_calc_step(num_models,1), total_steps(1),
                           output_acc_rate_freq(1000), output_acc_rate_window(1000), acc_rate_total(num_models,0), acc_rate_rolling(num_models,0), accepted_history(num_models, std::deque<unsigned int>(output_acc_rate_window, 0)),
                           accepted_history_count(num_models, 0),
                           acceptance_probability(0), move_type(0),
                           log_proposal_fwd(0), log_proposal_rev(0),
                           jump_move_probability(0.5), setup_ok(false), is_running(false), supress_foutput(false),
						   supress_mcmc_foutput(false), supress_prob_foutput(false), supress_chkpt_foutput(false),
						   steps_since_jump(0), min_steps_before_jump(5), min_jump_burnin(10), do_indep_jumps(true) {

            // check vector sizes
            if (sys_vec_.size() != get_LL_vec.size() || get_LL_vec.size() != log_prior_vec.size() || log_prior_vec.size() != mover_vec.size() || model_id_mapping.size() != sys_vec_.size()){
                std::cerr << "ERROR: MCMC::PTR::MODEL_SEL::mcmc_interface_v2(): functor pointer vector sizes do not match!" << std::endl;
                throw("ERROR: MCMC::PTR::MODEL_SEL::mcmc_interface_v2(): functor pointer vector sizes do not match!");
            }

        }
		
		void set_rand_gen(MTRand::MTRand_shared_ptr rand_ptr) {
			rand_gen = rand_ptr;
		}
		
		std::vector<unsigned int> get_mapping() const  {
			return model_id_mapping;
		}
		
		unsigned int get_current_model_mapped_id() const  {
			return model_id_mapping[current_model];
		}
		
		unsigned int get_current_model_unmapped_id() const  {
			return current_model;
		}
		
		typename system_ptr_::element_type::vector_type get_current_model_acc_params() const {
			return acc_params[this->current_model];
		}
		
		void set_current_model_acc_params(typename system_ptr_::element_type::vector_type params){
			proposed_model = current_model;
			prop_params[this->current_model] = params;
			acc_params[this->current_model] = params;
			sys_vec[this->current_model]->set_joint_para_and_state(params);
			evaluate_ll(this->current_model);
			this->accepted();
		}
		
		double get_current_acc_sum_ll() {
			// need to run evaluate in order to update
			this->evaluate_ll(current_model);
			return accepted_sum[current_model];
		}
		
		double get_inv_temp() const {
			return inv_temp;
		}
		
		double get_current_model_summed_ll(typename system_ptr_::element_type::vector_type params){
			return get_summed_ll(params, this->current_model);
		}



        // debugging function
        virtual double get_summed_ll( typename system_ptr_::element_type::vector_type params, const unsigned int model_number){
			if (sys_vec[model_number]->get_joint_para_and_state().size() != params.size()){
				*stderr_ << "ERROR!!!!!!!!!!!!!!!: get_summed_ll(): parameter sizes don't match!!!!!" << std::endl;
				return 99999999999;
			}
			else {
				sys_vec[model_number]->set_joint_para_and_state(params);
				return evaluate_ll(model_number);
			}
        }


        virtual double evaluate_ll(const unsigned int model_number){
            /*
            std::cout << "DEBUG_LL: 0" << std::endl;
            for (int ii = 0 ; ii < sys->get_joint_para_and_state().size(); ii++)
            {
                    std::cout << "\t" << sys->get_joint_para_and_state()[ii];
            }
            */


        	const bool bounds_ok = this->mover_vec[model_number]->check_bounds(this->sys_vec[model_number]->get_joint_para_and_state());
        	if (bounds_ok){

        		//std::cout << std::endl;
        		this_prior[model_number] = log_prior_vec[model_number]->operator()(sys_vec[model_number]->get_joint_para_and_state());
        		//std::cout << "DEBUG_LL: 1" << std::endl;
        		this_ll[model_number] = get_LL_vec[model_number]->operator()(sys_vec[model_number]);
        		//std::cout << "DEBUG_LL: 2" << std::endl;

        		//std::cout << "DEBUG_LL: 3" << std::endl;
        		this_sum[model_number] = this_ll[model_number] + this_prior[model_number];
        		//std::cout << "DEBUG_LL: 4" << std::endl;
        		return this_sum[model_number];
        		//return 1;
        	}
        	else {
        		this->this_prior[model_number] = std::numeric_limits<double>::lowest();
        		this->this_sum[model_number] = std::numeric_limits<double>::lowest();
        	}


        	return this->this_sum[model_number];

        }

        virtual int initialise(){
            current_step = ulong_vec(num_models, 0);
            covar_calc_step = ulong_vec(num_models, 0);
            this->jump_prop_count = Eigen::MatrixXd::Identity(num_models, num_models);
            this->jump_acc_count = Eigen::MatrixXd::Identity(num_models, num_models);
            for (current_model = 0; current_model < num_models; current_model++){
                        proposed_model = current_model;
                        acc_params[current_model] = sys_vec[current_model]->get_joint_para_and_state();
                        evaluate_ll(current_model);
                        accepted_ll[current_model] = this_ll[current_model];
                        accepted_prior[current_model] = this_prior[current_model];
                        accepted_sum[current_model] = this_sum[current_model];

                        best_accepted_ll_prior[current_model] = (accepted_ll[current_model] + accepted_prior[current_model]);
                        //acc_params = sys->get_joint_para_and_state(); // get initial conditions and parameters
                        best_params[current_model] = acc_params[current_model];
                        //this->mover_vec[current_model]->set_para_mean_vec(this->acc_params[current_model]);
                        this->mover_vec[current_model]->check_bounds(this->sys_vec[current_model]->get_joint_para_and_state(), true);

                        *stdout_ << "# " << output_root
                        		<< ":\tSTART\tstep_num:\t" << current_step[current_model]
								<< "\tcurrent_model:\t"<< this->model_id_mapping[current_model]
								<< "\tLL:\t" << best_accepted_ll_prior[current_model];

                        for (int ii = 0 ; ii < best_params[current_model].size(); ii++)
                        {
                        	*stdout_ << "\t" << best_params[current_model][ii];
                        }
                        *stdout_ << std::endl;


            }
			// start at a random model
            current_model = this->rand_gen->randInt(this->num_models - 1);//0;
			proposed_model = current_model;//0;
			steps_since_jump = 0;
            return 1;
        }

        virtual int initialise_output_files(){

            bool all_ok = true;
            output_chain.clear();
			if (!supress_mcmc_foutput){
            for (current_model = 0; current_model < num_models; current_model++){
                        std::string output_chain_filename = output_root + "_model_" + std::to_string(model_id_mapping[current_model]) + "_std_mcmc_chain.dat";
                        boost::filesystem::path p (output_chain_filename);

                        if (!boost::filesystem::exists(p))
                        {
                            output_chain.push_back(ofstream_shared_ptr(new std::ofstream(output_chain_filename.c_str(), std::ios::out)));
                        }
                        else {
                            std::cerr << "ERROR: overwriting file: " << output_chain_filename << std::endl;
                            all_ok = false;
                        }


                        if (output_chain.back()->is_open())
                        {
                            std::cout << "# writing MCMC chain to: " << output_chain_filename << std::endl;
                            //return 1;
                        }
                        else {
                            std::cerr << "ERROR: couldn't open file: " << output_chain_filename << std::endl;
                            all_ok = false;
                        }
            }
			}
            //current_model = 0;
            //proposed_model = 0;
			// start at a random model
			current_model = this->rand_gen->randInt(this->num_models - 1);//0;
			proposed_model = current_model;//0;
            return all_ok ? 1 : 0;
        }

        virtual bool indep_reversible_jump_move_move(){
        	 this->log_proposal_rev = 0;
        	            this->log_proposal_fwd = 0;



        	            if (this->num_models > 1){
        	                this->proposed_model = this->current_model;
        	                while (this->proposed_model == this->current_model){
        	                    this->proposed_model = this->rand_gen->randInt(this->num_models - 1);
        	                }
        	            }
        	            else {
        	                std::cerr << "ERROR: reversible_jump_move_move: you only have one model" << std::endl;
        	                return false;
        	            }

        	            //std::cout << "\nRJ: moving " << current_model << "->" << proposed_model << std::endl;

        	            //std::cout << "current size: " <<acc_params[current_model].size() << std::endl;

        	            Eigen::MatrixXd params_current = this->mover_vec[this->current_model]->get_mapped_parameters(this->acc_params[this->current_model]);//Eigen::MatrixXd::Zero(acc_params[current_model].size(), 1);
        	            /*
        	            for (unsigned int ii = 0; ii < acc_params[current_model].size(); ii++){
        	                params_current(ii) = acc_params[current_model][ii];
        	            }
        	            */

        	            Eigen::Matrix<double,Eigen::Dynamic, 1> mean_current = this->mover_vec[this->current_model]->get_mean();
        	            //cheat here:
        	            Eigen::Matrix<double,Eigen::Dynamic, 1> mean_proposed = this->mover_vec[this->proposed_model]->get_mean();
        	            //Eigen::Matrix<double,Eigen::Dynamic, 1> mean_proposed = this->mover_vec[proposed_model]->get_mapped_parameters(acc_params[proposed_model]);//this->mover_vec[proposed_model]->get_mean();

        	            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_current = this->mover_vec[this->current_model]->get_B_matrix(); // get_S_matrix
        	            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_proposed = this->mover_vec[this->proposed_model]->get_B_matrix(); //get_S_matrix

        	            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_current_inv = B_current.inverse();

        	            //std::cout << "B_curr:\n" << B_current << std::endl;
        	            //std::cout << "B_curr_inv:\n" << B_current_inv << std::endl;
        	            //std::cout << "params_curr:\n" << params_current << std::endl;
        	            //std::cout << "mean_curr:\n" << mean_current << std::endl;

        	            //std::cout << "B_proposed:\n" << B_proposed << std::endl;
        	            //std::cout << "mean_proposed:\n" << mean_proposed << std::endl;






        	            Eigen::Matrix<double,Eigen::Dynamic, 1> z_vec = B_current_inv * (params_current - mean_current);
        	            //std::cout << "z_vec:\n" << z_vec << std::endl;







        				// lets correct this to remove diagonal terms from B matrices corresponding to P1 and P2
        				const double log_jacobian = std::log(B_proposed.determinant() / B_current.determinant());
        				//const double new_log_jacobian = std::log(B_det_corrected_proposed/ B_det_corrected_current);
        				//const double new_new_correction = new_log_jacobian - log_jacobian;




        	            const int n_current = mean_current.rows();
        	            const int n_proposed = mean_proposed.rows();

        	            //const int n_diff = n_proposed - n_current;



        	            // random permutation matrix
        	            Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(std::max(n_proposed, n_current));
        	            perm.setIdentity();
        	            for (unsigned int ii = perm.cols() - 1; ii > 0; ii-- ){
        	                unsigned int jj = this->rand_gen->randInt(ii);
        	                std::swap(perm.indices().data()[ii], perm.indices().data()[jj]);
        	            }





        	            Eigen::Matrix<double,Eigen::Dynamic, 1> u_prop = Eigen::MatrixXd::Zero(n_proposed, 1);
        	            for (int ii = 0 ; ii < n_proposed; ii++){
        	            	u_prop(ii) = this->rand_gen->randNorm();
        	            }
        	            MultivariateNormalPDF pdf_prop;
        	            pdf_prop.setMean(Eigen::MatrixXd::Zero(n_proposed, 1));
        	            pdf_prop.setCovar(Eigen::MatrixXd::Identity(n_proposed, n_proposed)); //WTF??
        	            const double g_prop = pdf_prop.getLogPdf(u_prop);

        	            /*
        	            Eigen::Matrix<double,Eigen::Dynamic, 1> u_curr = Eigen::MatrixXd::Zero(n_current, 1);
        	            for (int ii = 0 ; ii < n_current; ii++){
        	            	u_curr(ii) = this->rand_gen->randNorm();
        	            }
        	            */
        	            MultivariateNormalPDF pdf_curr;
        	            pdf_curr.setMean(Eigen::MatrixXd::Zero(n_current, 1));
        	            pdf_curr.setCovar(Eigen::MatrixXd::Identity(n_current, n_current)); //WTF??

        	            const double g_curr = pdf_curr.getLogPdf(z_vec);


        	            const double log_G = g_curr - g_prop;
        	            Eigen::Matrix<double,Eigen::Dynamic, 1> proposed_params = mean_proposed + (B_proposed * u_prop);//Eigen::MatrixXd::Zero(n_proposed, 1);



        	            //std::cout << "proposed_params:\n" << proposed_params << std::endl;

        	            this->prop_params[this->proposed_model] = this->mover_vec[this->proposed_model]->map_changes_back(this->acc_params[this->proposed_model], proposed_params);
        	            this->sys_vec[this->proposed_model]->set_joint_para_and_state(this->prop_params[this->proposed_model]);



        				/*(this->mover_vec[this->proposed_model]->get_upper_bound(P1[this->proposed_model]) - this->mover_vec[this->proposed_model]->get_lower_bound(P1[this->proposed_model]))
        	                                            * (this->mover_vec[this->proposed_model]->get_upper_bound(P2[this->proposed_model]) - this->mover_vec[this->proposed_model]->get_lower_bound(P2[this->proposed_model]));*/

        				//this->stdout_ << "current_model: " << this->current_model << " vol: " << current_hist_volume << " proposed_model: " << this->proposed_model <<" vol: " << proposed_hist_volume << std::endl;







        				//const double new_log_final = new_log_jacobian + log_G;
        				const double log_final_no_corr = log_jacobian + log_G ;
        				const double log_final = log_final_no_corr;// + dim_correction;//new_new_correction;




        				this->log_proposal_rev = log_final;//old_log_final_no_corr;//old_log_final;//new_log_final;//old_log_final;//log_jacobian + log_G + correction;//+ new_correction;// - std::log(proposed_hist_volume);
        				this->log_proposal_fwd = 0;//-std::log(current_hist_volume);


        				/*
        				if (this->total_steps % 1000 == 0) {
        					std::cout << this->output_root << " MODS:" << this->current_model << "->" << this->proposed_model
        					<< " OLDFINAL_correction: " << old_log_final
        					<< " NEWFINAL_correction: " << new_log_final
        					<< std::endl;
        				}
        				 */


        	            //std::cout << "log_jacobian: " << log_jacobian << " log_G: " << log_G << std::endl;

        	            return true;

        }

        virtual bool reversible_jump_move_move(){
			
			
			
            this->log_proposal_rev = 0;
            this->log_proposal_fwd = 0;



            if (num_models > 1){
                proposed_model = current_model;
                while (proposed_model == current_model){
                    proposed_model = this->rand_gen->randInt(this->num_models - 1);
                }
            }
            else {
                std::cerr << "ERROR: reversible_jump_move_move: you only have one model" << std::endl;
                return false;
            }

            //std::cout << "\nRJ: moving " << current_model << "->" << proposed_model << std::endl;

            //std::cout << "current size: " <<acc_params[current_model].size() << std::endl;

            Eigen::MatrixXd params_current = this->mover_vec[current_model]->get_mapped_parameters(acc_params[current_model]);//Eigen::MatrixXd::Zero(acc_params[current_model].size(), 1);
            /*
            for (unsigned int ii = 0; ii < acc_params[current_model].size(); ii++){
                params_current(ii) = acc_params[current_model][ii];
            }
            */

            Eigen::Matrix<double,Eigen::Dynamic, 1> mean_current = this->mover_vec[current_model]->get_mean();
            //cheat here:
            Eigen::Matrix<double,Eigen::Dynamic, 1> mean_proposed = this->mover_vec[proposed_model]->get_mean();
            //Eigen::Matrix<double,Eigen::Dynamic, 1> mean_proposed = this->mover_vec[proposed_model]->get_mapped_parameters(acc_params[proposed_model]);//this->mover_vec[proposed_model]->get_mean();

            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_current = this->mover_vec[current_model]->get_B_matrix(); //  ->get_S_matrix();
            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_proposed = this->mover_vec[proposed_model]->get_B_matrix(); //  ->get_S_matrix();

            Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B_current_inv = B_current.inverse();

            //std::cout << "B_curr:\n" << B_current << std::endl;
            //std::cout << "B_curr_inv:\n" << B_current_inv << std::endl;
            //std::cout << "params_curr:\n" << params_current << std::endl;
            //std::cout << "mean_curr:\n" << mean_current << std::endl;
            //std::cout << "mean_proposed:\n" << mean_proposed << std::endl;






            Eigen::Matrix<double,Eigen::Dynamic, 1> z_vec = B_current_inv * (params_current - mean_current);
            //std::cout << "z_vec:\n" << z_vec << std::endl;

            const double log_jacobian = std::log(B_proposed.determinant() / B_current.determinant());



            const int n_current = mean_current.rows();
            const int n_proposed = mean_proposed.rows();

            const int n_diff = n_proposed - n_current;


            // random permutation matrix
            Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> perm(std::max(n_proposed, n_current));
            perm.setIdentity();
            for (unsigned int ii = perm.cols() - 1; ii > 0; ii-- ){
                unsigned int jj = rand_gen->randInt(ii);
                std::swap(perm.indices().data()[ii], perm.indices().data()[jj]);
            }
            //std::cout << "\nperm:\n" << perm.toDenseMatrix() << std::endl;


            double log_G = 0;
            Eigen::Matrix<double,Eigen::Dynamic, 1> proposed_params = Eigen::MatrixXd::Zero(n_proposed, 1);
            if (n_proposed < n_current){
                    proposed_params = mean_proposed + (B_proposed * (perm * z_vec).head(n_proposed));
                    Eigen::Matrix<double,Eigen::Dynamic, 1> u = (perm * z_vec).tail(n_diff);//Eigen::MatrixXd::Zero(std::abs(n_diff), 1);
                    /*
                    for (int ii = 0 ; ii< std::abs(n_diff); ii++){
                        u(ii) = rand_gen->randNorm();
                    }
                    */
                    MultivariateNormalPDF pdf;
                    pdf.setMean(Eigen::MatrixXd::Zero(std::abs(n_diff), 1));
                    pdf.setCovar(Eigen::MatrixXd::Identity(std::abs(n_diff), 1));
                    log_G = pdf.getLogPdf(u);

            }
            else if (n_proposed == n_current){
                proposed_params = mean_proposed + (B_proposed * (perm*z_vec));
                log_G = 0;
            }
            else {
                // n_proposed > n_current
                //proposed_params = Eigen::MatrixXd::Zero(n_proposed, 1);
                Eigen::Matrix<double,Eigen::Dynamic, 1> u = Eigen::MatrixXd::Zero(std::abs(n_diff), 1);
                for (int ii = 0 ; ii< std::abs(n_diff); ii++){
                    u(ii) = rand_gen->randNorm();
                }
                for (int ii = 0 ; ii< std::abs(n_diff); ii++){
                        u(ii) = rand_gen->randNorm();
                    }
                    MultivariateNormalPDF pdf;
                    pdf.setMean(Eigen::MatrixXd::Zero(std::abs(n_diff), 1));
                    pdf.setCovar(Eigen::MatrixXd::Identity(std::abs(n_diff), std::abs(n_diff)));
                    log_G = -pdf.getLogPdf(u);
                proposed_params << z_vec, u;
                //std::cout << "proposed_params before:\n" << proposed_params << std::endl;
                proposed_params = mean_proposed + (B_proposed * (perm*proposed_params));

            }

            //std::cout << "proposed_params:\n" << proposed_params << std::endl;

            prop_params[proposed_model] = mover_vec[proposed_model]->map_changes_back(acc_params[proposed_model], proposed_params);
            sys_vec[proposed_model]->set_joint_para_and_state(prop_params[proposed_model]);

            this->log_proposal_rev = log_jacobian + log_G;


            //std::cout << "log_jacobian: " << log_jacobian << " log_G: " << log_G << std::endl;

            return true;
        }

        virtual bool switch_model_move(){

            this->log_proposal_rev = 0;
            this->log_proposal_fwd = 0;

            if (num_models > 1){
                proposed_model = current_model;
                while (proposed_model == current_model){
                    proposed_model = this->rand_gen->randInt(this->num_models - 1);
                }
            }

            prop_params[proposed_model] = acc_params[proposed_model];


            sys_vec[proposed_model]->set_joint_para_and_state(prop_params[proposed_model]);

            //std::cout << "moving " << current_model << "->" << proposed_model << std::endl;

            return true;

            /*
            proposed_model = this->rand_gen->randInt(this->num_models - 1);
            if (proposed_model == current_model){
                return false;
            }
            bool accept = false;
            if (accepted_sum[proposed_model] > accepted_sum[current_model]){
                accept = true;
            }
            else {
                double p_acc = std::exp(inv_temp*(accepted_sum[proposed_model] - accepted_sum[current_model]));
                double rand = this->rand_gen->rand();
                    if (rand < p_acc)
                    {
                        accept = true;
                        //std::cout << "accepted:\t" << p_acc << "\taccepted_sum:\t" << accepted_sum << "\tnew_ll:\t" << this_sum << std::endl;
                    }
                    else
                    {
                        //std::cout << "rejected:\t" << p_acc << "\taccepted_sum:\t" << accepted_sum << "\tnew_ll:\t" << this_sum << std::endl;
                    }
            }

            if (accept == true){
                //std::cout << "# INFO: switching_models" << "\told_model_num:\t" << current_model << "\tnew_model_num:\t" << proposed_model
                //        << "\tdiff:\t" << ((accepted_sum[proposed_model] - accepted_sum[current_model])) << std::endl;
                current_model = proposed_model;
                return true;
            }
            //proposed_model = current_model;
            return false;
            */

        }

        virtual int proposal(){

            this->log_proposal_rev = 0;
            this->log_proposal_fwd = 0;

            const unsigned int n_dims = this->mover_vec[this->current_model]->get_n_dims();
            const unsigned long jump_burnin = std::max(n_dims*this->min_steps_before_jump, this->min_jump_burnin);


            move_type = 0;
            if (num_models > 1 && (this->steps_since_jump > jump_burnin)){
                const double prob = rand_gen->rand();
                if (prob < this->jump_move_probability){
                    move_type = this->rand_gen->randInt(1);
                }
            }

            if (move_type == 0){
                // normal move
                //std::cout << "current_model: " << current_model << std::endl;
                //std::cout << "acc_params.size(): " << acc_params[current_model].size() << std::endl;
                prop_params[current_model] = mover_vec[current_model]->operator()(acc_params[current_model]);
                sys_vec[current_model]->set_joint_para_and_state(prop_params[current_model]);
                proposed_model = current_model;
            }
            else if (move_type == 1) {
            	if (this->do_indep_jumps){
            		indep_reversible_jump_move_move();
            	}
            	else {
            		reversible_jump_move_move();
            	}
                //switch_model_move(sys_vec);
            }

            /*
            std::cout << "step:" << total_steps << "\tDEBUG:proposal:";
            for (int ii = 0 ; ii < prop_params[proposed_model].size(); ii++){
            	std::cout << "\t" << prop_params[proposed_model][ii];
            }
            std::cout << std::endl;
            */

            return 1;

            /*  OLD Method
            proposed_model = current_model;//this->rand_gen->randInt(this->num_models - 1);
            //std::cout << "DEBUG: prop " << current_model << "->" << proposed_model << std::endl;
            //std::cout << "in: " << acc_params[proposed_model][0] << " " << acc_params[proposed_model][1] << std::endl;
            prop_params[proposed_model] = mover_vec[proposed_model]->operator()(acc_params[proposed_model]);
            sys_vec[proposed_model]->set_joint_para_and_state(prop_params[proposed_model]);
            return 1;
            */
        }

        virtual bool to_accept(){
            bool accept = false;
			acceptance_probability = 0;

            const bool bounds_ok = this->mover_vec[this->proposed_model]->check_bounds(this->prop_params[this->proposed_model]);
            if (bounds_ok == false){
                //acceptance_probability = 0;
                return false;
            }

            const double proposal_diff = this->log_proposal_rev - this->log_proposal_fwd;

            const double prop_sum = this_sum[proposed_model] + proposal_diff;//move_type == 0 ? this_sum[proposed_model] : accepted_sum[proposed_model];



            //std::cout << "\nto_accept:start" << std::endl;
			if (!std::isinf(prop_sum) && !std::isnan(prop_sum))
            {
                //const double debug_curr_sum = get_summed_ll(sys_vec, acc_params[current_model], current_model);
                //const double debug_new_sum = get_summed_ll(sys_vec, prop_params[proposed_model], proposed_model);
                //std::cout << "DEBUG: debug_curr_sum: " << debug_curr_sum << " -> " << debug_new_sum << std::endl;


                //if (prop_sum != debug_new_sum){
                //    std::cout << "ERROR!!!! sum do not match up!!"   << std::endl;
                //}


                if (prop_sum > accepted_sum[current_model])
                {
                    accept = true;
                    acceptance_probability = 1;
                    /*
                    std::cout << "moving " << current_model << "->" << proposed_model << " accepted_better:\t" << "\told_ll:\t" << accepted_sum[current_model]
                        << "\prop_sum:\t" << prop_sum
                        << "\tnew_ll:\t" << this_sum[proposed_model]
                        << std::endl;
                        */
                }
                else
                {
                    double p_acc = std::exp(inv_temp*(prop_sum - accepted_sum[current_model]));
                    acceptance_probability = p_acc;
                    double rand = this->rand_gen->rand();
                    if (rand < p_acc)
                    {
                        accept = true;
                        /*
                        std::cout << "moving " << current_model << "->" << proposed_model << " accepted:\t" << p_acc << "\taccepted_sum:\t" << accepted_sum[current_model]
                            << "\tprop_sum:\t" << prop_sum
                            << "\tnew_ll:\t" << this_sum[proposed_model]
                            << std::endl;
                            */
                    }
                    else
                    {
                        /*
                        std::cout << "moving " << current_model << "->" << proposed_model << " rejected:\t" << p_acc << "\taccepted_sum:\t" << accepted_sum[current_model]
                            << "\tprop_sum:\t" << prop_sum
                            << "\tnew_ll:\t" << this_sum[proposed_model]
                            << std::endl;
                            */
                    }
                }
            }
			else{
				std::cerr << "# WARNING: this_sum[this->proposed_model] is nan or inf" << std::endl;
			}

            return accept;
        }

        virtual int post_metropolis_eval(const bool accept){
            // this needs thinking about
            if (proposed_model == current_model){
                        accepted_history[current_model].push_back(accept ? 1 : 0);
                        accepted_history_count[current_model] += (accept ? 1 : 0);
                        accepted_history_count[current_model] -= accepted_history[current_model].front();
                        accepted_history[current_model].pop_front();
                        acc_rate_rolling[current_model] = static_cast<double>(accepted_history_count[current_model])/ static_cast<double>(output_acc_rate_window);
            }
            accepted_count += accept ? 1 : 0;
            acc_rate_total[current_model] = static_cast<double>(accepted_count) / static_cast<double>(total_steps);

            if (this->move_type == 1){


            	// jump move
            	this->jump_prop_count(current_model, proposed_model)++;
            	if (accept){
            		this->jump_acc_count(current_model, proposed_model)++;
            	}

            }

            return 1;
        }

        virtual int post_metropolis_eval_output(const bool accept){
            // this needs thinking about
            if (current_step[current_model] % output_acc_rate_freq == 0){
                *stdout_ << "# " << output_root << ":\tstep_num:\t" << current_step[current_model] << "\tcurrent_model:\t"<< this->model_id_mapping[current_model] << "\taccepted_num:\t" << accepted_count << "\tAcceptance_rate:\t" << acc_rate_total[current_model]
                        << "\tAcceptance_rate_running:\t" << acc_rate_rolling[current_model]
                        << std::endl;
                Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> rates_mat =  this->jump_acc_count.cwiseProduct(this->jump_prop_count.cwiseInverse());

                //*stdout_ << "jump_prop_count:\n" << this->jump_prop_count << std::endl;
                //*stdout_ << "jump_acc_count:\n" << this->jump_acc_count << std::endl;
                *stdout_ << "jump_transition_rates:\n" << rates_mat << std::endl;

                /*
                for (unsigned int ii = 0 ; ii < num_models; ii++){
                	for (unsigned int jj = 0 ; jj < num_models; jj++){
                		const  int n_dims_curr = this->mover_vec[ii]->get_n_dims();
                		const  int n_dims_prop = this->mover_vec[jj]->get_n_dims();
                		const  int n_diff = n_dims_prop - n_dims_curr;
                		std::cout << "DEBUG_FACTOR:\t"
                				<<   n_dims_curr << "\t"
								<< n_dims_prop << "\t"
                				<< rates_mat(ii,jj) << "\t"
								<< rates_mat(jj,ii) <<  "\t"
								<< n_diff  << "\t"
								<< rates_mat(ii,jj) / rates_mat(jj,ii) << "\t"
								<< std::log(rates_mat(ii,jj) / rates_mat(jj,ii)) << "\t"
								<< std::endl;
                	}
                }
                */

            }



            /*
            stdout_ << "#accepted_energies";
            for (unsigned int ii = 0; ii < num_models; ii++){
                stdout_ << "\t" << accepted_sum[ii];
            }
            stdout_ << std::endl;
            */

            return 1;
        }

        virtual int accepted(){
            last_model = current_model;
            current_model = proposed_model;
            //if (move_type == 0){
            accepted_ll[proposed_model] = this_ll[proposed_model];
            accepted_prior[proposed_model] = this_prior[proposed_model];
            accepted_sum[proposed_model] = this_sum[proposed_model];
            acc_params[proposed_model] = prop_params[proposed_model];//sys->get_joint_para_and_state();
            //}
            // else {
            // do nothing
            //}



            return 1;
        }

        virtual int accepted_output(){
            if ( (accepted_ll[current_model] + accepted_prior[current_model]) > best_accepted_ll_prior[current_model]){
                    best_accepted_ll_prior[current_model] = (accepted_ll[current_model] + accepted_prior[current_model]);
                    best_params[current_model] = acc_params[current_model];
                *stdout_ << "# " << output_root << ":\tBEST\tstep_num:\t" << current_step[current_model] << "\tcurrent_model:\t"<< this->model_id_mapping[current_model] << "\tLL:\t" << best_accepted_ll_prior[current_model];

                for (int ii = 0 ; ii < best_params[current_model].size(); ii++)
                {
                    *stdout_ << "\t" << best_params[current_model][ii];
                }
                *stdout_ << std::endl;
            }

            return 1;
        }


        virtual int rejected(){

            return 1;
        }

        virtual int rejected_output(){
        	/*
        	*stdout_ << "# " << output_root << ":\REJECTED\tstep_num:\t" << current_step[current_model] << "\tcurrent_model:\t"<< this->model_id_mapping[current_model];// << "\tLL:\t" << best_accepted_ll_prior[current_model];

        	for (int ii = 0 ; ii < prop_params[proposed_model].size(); ii++)
        	{
        		*stdout_ << "\t" << prop_params[proposed_model][ii];
        	}
        	*stdout_ << std::endl;
        	*/
            return 1;
        }



        virtual int post_update(const bool accept){


        	if (accept == true && move_type != 0) {
        		this->steps_since_jump = 0;
        	}
        	else {
        		this->steps_since_jump++;
        	}

            return 1;
        }

        virtual int chain_output(const std::string line_prefix, const bool accept)
        {
        	if (!supress_mcmc_foutput){
        		(*(output_chain[current_model])) << line_prefix << "ACCEPTED" << (accept ? 1 : 0) << "\tstep_num:\t" << current_step[this->current_model]  << "\tLL:\t" << accepted_ll[current_model]
																																													   << "\tlog_prior:\t" << accepted_prior[current_model];

        		for (int ii = 0 ; ii < acc_params[current_model].size(); ii++)
        		{
        			(*(output_chain[current_model])) << "\t" << acc_params[current_model][ii];
        		}
        		(*(output_chain[current_model])) << std::endl;
        	}

        	return 1;
        }

        virtual int post_update_output(const bool accept){
            //std::cout << "blah1:" << output_freq << "\t" << current_step << "\t" << burn_in << std::endl;
            if ((current_step[current_model] % (output_freq) == 0) && current_step[current_model] >= burn_in){
                //std::cout << "blah2" << std::endl;
                chain_output("",  accept);

				if (current_step[current_model] % (output_freq*100) == 0){
					*stdout_ << "# " << output_root << ":\taccepted_sumLL:\ttotal_steps:\t" << this->total_steps;
					for (unsigned int ii = 0; ii < num_models; ii++)
					{
						*stdout_ << "\t" << accepted_sum[ii];
					}
					*stdout_ << std::endl;
				}



            }
            return 1;
        }

        virtual int final_call(){
        	for (unsigned int ii = 0; ii < num_models; ii++)
        	{
        		sys_vec[ii]->set_joint_para_and_state(acc_params[ii]);
        		// *stdout_ << "\t" << accepted_sum[ii];
        	}
        	return 1;
        }


		void set_running_flag(bool val){
			boost::mutex::scoped_lock lock(is_running_mutex);
			is_running = val;
		}
		
		bool has_finished() {
			boost::mutex::scoped_lock lock(is_running_mutex);
			return !is_running;
		}

		int setup(){

			if (supress_foutput == true){
				supress_mcmc_foutput = true;
				supress_prob_foutput = true;
				supress_chkpt_foutput = true;
			}

			total_steps = 0;
			initialise();
			if (initialise_output_files() != 1){
				setup_ok = false;
				return -1;
			}
			chain_output("#\t", 1);
			setup_ok = true;
			return 1;
		}


        int operator()(  unsigned long steps_to_run ){
			
			{
			boost::mutex::scoped_lock lock(is_running_mutex);
			is_running = true;
			}

            //this->rand_gen = rand_ptr;

            //std::cout << "DEBUG: " << std::endl;
			/*
            if (sys_vec.size() != num_models){
                stderr_ << "ERROR: sys_vec.size()=" << sys_vec.size() << " is different from the required size from the initialisation parameters: " << num_models << std::endl;
                return -1;
            }
			*/
			if (setup_ok == false){
				*stderr_ << "ERROR: there was a problem with your setup - (e.g. are you overwriting exisiting files or did you not run setup()?)." << std::endl;
				return -1;
			}
			
			if (steps_to_run == std::numeric_limits< unsigned long >::max()){
				steps_to_run = std::numeric_limits< unsigned long >::max() - total_steps - 1;
			}

			const unsigned long start_step = total_steps+1;
			const unsigned long end_step = std::min(steps, total_steps+steps_to_run);


            for  (total_steps = start_step; total_steps <= end_step; total_steps++){
                //std::cout << "DEBUG: 1: step: " << total_steps << std::endl;
                //switch_model_move();
                proposal();
                //std::cout << "DEBUG: 2 " << current_model << "->" << proposed_model << std::endl;
                try {
                    //std::cout << "DEBUG: 3" << std::endl;
                    evaluate_ll( proposed_model);
                    //std::cout << "DEBUG: 4" << std::endl;
                }
                catch (const std::exception& e){
                    *stderr_ << "# WARNING: exception caught on evaluate_ll():\t" << e.what() << "# " << std::endl;
                    this_sum[current_model] =  std::numeric_limits< double>::lowest();
                }
                //std::cout << "DEBUG: 5" << std::endl;
                const bool accept = to_accept();
                //std::cout << "DEBUG: 6" << std::endl;
                post_metropolis_eval( accept);
                //std::cout << "DEBUG: 7" << std::endl;
                post_metropolis_eval_output( accept);
                //std::cout << "DEBUG: 8" << std::endl;

                if (accept){
                    //std::cout << "DEBUG: 9" << std::endl;
                    accepted();
                    //std::cout << "DEBUG: 10" << std::endl;
                    accepted_output( );
                    //std::cout << "DEBUG: 11" << std::endl;
                }
                else {
                    //std::cout << "DEBUG: 12" << std::endl;
                    rejected();
                    //std::cout << "DEBUG: 13" << std::endl;
                    rejected_output( );
                    //std::cout << "DEBUG: 14" << std::endl;
                }
                //std::cout << "DEBUG: 15" << std::endl;
                post_update( accept);
                //std::cout << "DEBUG: 16" << std::endl;
                post_update_output( accept);
                //std::cout << "DEBUG: 17" << std::endl;
                current_step[current_model]++;
            }

			{
			boost::mutex::scoped_lock lock(is_running_mutex);
			is_running = false;
			}

			final_call();

            return 1;
        }




    };


    template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
    class adaptive_mcmc_interface_v2: public mcmc_interface_v2<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_>  {

    public:


        typedef std::vector<system_ptr_> system_ptr_vector;
        typedef std::vector<ll_functor_ptr_> ll_functor_ptr_vector;
        typedef std::vector<prior_ptr_> prior_ptr_vector;
        typedef std::vector<proposal_ptr_> proposal_ptr_vector;

        typedef std::vector<double> double_vec;
        typedef std::vector<unsigned long> ulong_vec;

        bool greedy_start;
        ulong_vec adaptation_step;
        unsigned long min_adaptation_step;

        adaptive_mcmc_interface_v2( MTRand_shared_ptr rand_ptr, std::ostream& stdout__, std::ostream& stderr__, uint_vec model_id_mapping_,
								   system_ptr_vector sys_vec_, ll_functor_ptr_vector get_LL_vec_, prior_ptr_vector log_prior_vec_, proposal_ptr_vector mover_vec_, std::string output_root_,
                       unsigned long steps_, unsigned long burn_in_, unsigned long output_freq_, const double inv_temp_) : mcmc_interface_v2<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_>(
                       rand_ptr, stdout__, stderr__, model_id_mapping_, sys_vec_, get_LL_vec_, log_prior_vec_, mover_vec_, output_root_,
                       steps_, burn_in_, output_freq_, inv_temp_), greedy_start(true), adaptation_step(this->num_models, 1), min_adaptation_step(10)  {
            //std::cout << this->burn_in << std::endl;
        }



        virtual int post_update( const bool accept) override{
            if (this->move_type == 0){//if (this->last_model == this->proposed_model ){
                        if (this->greedy_start && this->current_step[this->current_model] < this->burn_in)
                        {
                            if (accept)
                            {
                                this->mover_vec[this->proposed_model]->adapt(std::max(adaptation_step[this->proposed_model], min_adaptation_step), this->acceptance_probability, this->acc_params[this->proposed_model], accept);
                                adaptation_step[this->proposed_model]++;
                            }
                        }
                        else {
                            this->mover_vec[this->proposed_model]->adapt(std::max(adaptation_step[this->proposed_model], min_adaptation_step), this->acceptance_probability, this->acc_params[this->proposed_model], accept);
                            adaptation_step[this->proposed_model]++;
                        }
            }

            /**********************************/
            // record population mean and covar
            /**********************************/

            this->mover_vec[this->current_model]->update_para_covar(this->covar_calc_step[this->current_model], this->acc_params[this->current_model]);
            this->covar_calc_step[this->current_model]++;
            /**********************************/
            /**********************************/

            return 1;
        }

        virtual int post_update_output(const bool accept) override {
        	//std::cout << "blah1:" << output_freq << "\t" << current_step << "\t" << burn_in << std::endl;
        	if ((this->current_step[this->current_model] % (this->output_freq) == 0) && this->current_step[this->current_model] >= this->burn_in){
        		//std::cout << "blah2" << std::endl;
        		this->chain_output("",  accept);

        		if (this->current_step[this->current_model] % (this->output_freq*100) == 0){
        			*(this->stdout_) << "# " << this->output_root << ":\taccepted_sumLL:\ttotal_steps:\t" << this->total_steps;
        			for (unsigned int ii = 0; ii < this->num_models; ii++)
        			{
        				*(this->stdout_) << "\t" << this->accepted_sum[ii];
        			}
        			*(this->stdout_) << std::endl;
        		}

        		if ((this->current_step[this->current_model] % (this->output_freq*10) == 0)
        				&& (!this->supress_foutput) && (!this->supress_mcmc_foutput)){
        			const std::string chk_covar_filename = this->output_root
        					+ "_model_" + std::to_string(this->model_id_mapping[this->current_model]) + "_covar.dat";



        			std::ofstream chk_covar_file(chk_covar_filename.c_str(), std::ios::out| std::ios::trunc);
        			chk_covar_file << this->mover_vec[this->current_model]->get_covar() << std::endl;
        			chk_covar_file.close();
        		}



        	}
        	return 1;
        }





    };




} // namespace MODEL_SEL
} // namespace PTR
} // namespace MCMC
#endif /* mcmc_ptr_mod_sel_algorithms_h */
