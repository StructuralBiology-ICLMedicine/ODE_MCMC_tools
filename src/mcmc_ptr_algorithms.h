//
//  mcmc_algorithms.h
//  ode_mcmc_code
//
//  Created by jmacdon on 11/09/2016.
//  Copyright © 2016 James T. MacDonald. All rights reserved.
//

#ifndef mcmc_algorithms_h
#define mcmc_algorithms_h




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






namespace MCMC{
namespace PTR{


    typedef std::shared_ptr<std::ofstream> ofstream_shared_ptr;
    typedef std::shared_ptr<MTRand> MTRand_shared_ptr;


    template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
    class mcmc_interface{

    public:



    public:

        MTRand_shared_ptr rand_gen;

        ll_functor_ptr_ get_LL;
        prior_ptr_ log_prior;
        proposal_ptr_ mover;


        std::string output_root;
        typename system_ptr_::element_type::vector_type acc_params;
        typename system_ptr_::element_type::vector_type prop_params;
		typename system_ptr_::element_type::vector_type best_params;

        unsigned long steps;
        unsigned long burn_in;
        unsigned long output_freq;
        double inv_temp;

        double accepted_ll, accepted_prior, accepted_sum;
        double best_accepted_ll_prior;

        double this_ll, this_prior, this_sum;

        unsigned long accepted_count;
		unsigned long current_step;

		unsigned long output_acc_rate_freq;
		unsigned long output_acc_rate_window;

		double acc_rate_total, acc_rate_rolling;
		std::deque<unsigned int> accepted_history;
		unsigned long accepted_history_count = 0;

		double acceptance_probability;


		ofstream_shared_ptr output_chain;

		bool supress_foutput;


    public:
        mcmc_interface(system_ptr_ sys, ll_functor_ptr_ get_LL_, prior_ptr_ log_prior_, proposal_ptr_ mover_, std::string output_root_,
                       unsigned long steps_, unsigned long burn_in_, unsigned long output_freq_, const double inv_temp_) :
                           get_LL(get_LL_), log_prior(log_prior_), mover(mover_),
                           output_root(output_root_), acc_params(), prop_params(), best_params(),
                           steps(steps_), burn_in(burn_in_), output_freq(output_freq_), inv_temp(inv_temp_),
                           accepted_ll(0), accepted_prior(0), accepted_sum(0), best_accepted_ll_prior(0),
                           this_ll(0), this_prior(0), this_sum(0),
                           accepted_count(0), current_step(0),
                           output_acc_rate_freq(1000), output_acc_rate_window(1000), acc_rate_total(0), acc_rate_rolling(0), accepted_history(output_acc_rate_window, 0), accepted_history_count(0),
                           acceptance_probability(0), supress_foutput(false)    {
        }


        virtual int evaluate_ll(system_ptr_ sys){
            /*
            std::cout << "DEBUG_LL: 0" << std::endl;
            for (int ii = 0 ; ii < sys->get_joint_para_and_state().size(); ii++)
            {
                    std::cout << "\t" << sys->get_joint_para_and_state()[ii];
            }
            */
            //std::cout << std::endl;
            const bool bounds_ok = this->mover->check_bounds(this->prop_params);

            if (bounds_ok){
                this_prior = log_prior->operator()(sys->get_joint_para_and_state());
                //std::cout << "DEBUG_LL: 1" << std::endl;
                this_ll = get_LL->operator()(sys);
                //std::cout << "DEBUG_LL: 2" << std::endl;

                //std::cout << "DEBUG_LL: 3" << std::endl;
                this_sum = this_ll + this_prior;
                //std::cout << "DEBUG_LL: 4" << std::endl;
                return 1;
            }
            else {
                this_prior = std::numeric_limits<double>::lowest();
                this_sum = std::numeric_limits<double>::lowest();
                return 1;
            }
        }

        virtual int initialise(system_ptr_ sys){
            current_step = 0;
            acc_params = sys->get_joint_para_and_state();
            evaluate_ll(sys);
            accepted_ll = this_ll;
            accepted_prior = this_prior;
            accepted_sum = this_sum;

            best_accepted_ll_prior = (accepted_ll + accepted_prior);
            //acc_params = sys->get_joint_para_and_state(); // get initial conditions and parameters
            best_params = acc_params;
            return 1;
        }

        virtual int initialise_output_files(system_ptr_ sys, std::ostream& stdout_, std::ostream& stderr_){

                output_chain = 0;

                if (!this->supress_foutput)
                {
                    std::string output_chain_filename = output_root + "_std_mcmc_chain.dat";
                    boost::filesystem::path p (output_chain_filename);

                    if (!boost::filesystem::exists(p))
                    {
                        output_chain = ofstream_shared_ptr(new std::ofstream(output_chain_filename.c_str(), std::ios::out));
                    }
                    else
                    {
                        std::cerr << "ERROR: overwriting file: " << output_chain_filename << std::endl;
                        return 0;
                    }


                    if (output_chain->is_open())
                    {
                        std::cout << "# writing MCMC chain to: " << output_chain_filename << std::endl;
                        return 1;
                    }

                    std::cerr << "ERROR: couldn't open file: " << output_chain_filename << std::endl;

                    return 0;
                }
                else {
                    std::cout << "# mcmc_interface: Suppressing file output" << std::endl;
                    return 1;
                }
                return 0;
        }

        virtual int proposal(system_ptr_ sys){
            prop_params = mover->operator()(acc_params);
            sys->set_joint_para_and_state(prop_params);
            return 1;
        }

        virtual bool to_accept(system_ptr_ sys){

			
			acceptance_probability = 0;
            const bool bounds_ok = this->mover->check_bounds(this->prop_params);
            if (bounds_ok == false){
                //acceptance_probability = 0;
                return false;
            }

            bool accept = false;
			if (!std::isinf(this_sum) && !std::isnan(this_sum))
            {
                if (this_sum > accepted_sum)
                {
                    accept = true;
                    acceptance_probability = 1;
                    //cout << "accepted_better:\t" << "\told_ll:\t" << acc_ll << "\tnew_ll:\t" << this_ll << endl;
                }
                else
                {
                    double p_acc = exp(inv_temp*(this_sum - accepted_sum));
                    acceptance_probability = p_acc;
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
            }
			else {
				std::cerr << "# WARNING: this_sum[this->proposed_model] is nan or inf" << std::endl;
				
				std::cerr << "# prop_params:\t";
				for (int ii = 0 ; ii < this->prop_params.size(); ii++)
				{
					std::cerr << "\t" << this->prop_params[ii];
				}
				std::cerr << std::endl;
			}

            return accept;
        }

        virtual int post_metropolis_eval(system_ptr_ sys, const bool accept){
            accepted_history.push_back(accept ? 1 : 0);
            accepted_history_count += (accept ? 1 : 0);
            accepted_history_count -= accepted_history.front();
            accepted_history.pop_front();
            acc_rate_rolling = static_cast<double>(accepted_history_count)/ static_cast<double>(output_acc_rate_window);
            accepted_count += accept ? 1 : 0;
            acc_rate_total = static_cast<double>(accepted_count) / static_cast<double>(current_step);
            return 1;
        }

        virtual int post_metropolis_eval_output(system_ptr_ sys, const bool accept, std::ostream& stdout_, std::ostream& stderr_){
            if (current_step % output_acc_rate_freq == 0){
                stdout_ << "#step_num:\t" << current_step << "\taccepted_num:\t" << accepted_count << "\tAcceptance_rate:\t" << acc_rate_total
                        << "\tAcceptance_rate_running:\t" << acc_rate_rolling
                        << std::endl;
            }
            return 1;
        }

        virtual int accepted(system_ptr_ sys){
            accepted_ll = this_ll;
            accepted_prior = this_prior;
            accepted_sum = this_sum;
            acc_params = prop_params;//sys->get_joint_para_and_state();
            return 1;
        }

        virtual int accepted_output(system_ptr_ sys, std::ostream& stdout_, std::ostream& stderr_){
            if ( (accepted_ll + accepted_prior) > best_accepted_ll_prior){
                    best_accepted_ll_prior = (accepted_ll + accepted_prior);
                    best_params = acc_params;
                stdout_ << "# BEST\tstep_num:\t" << current_step << "\tLL:\t" << best_accepted_ll_prior;

                for (int ii = 0 ; ii < best_params.size(); ii++)
                {
                    stdout_ << "\t" << best_params[ii];
                }
                stdout_ << std::endl;
            }
            return 1;
        }


        virtual int rejected(system_ptr_ sys){
            return 1;
        }

        virtual int rejected_output(system_ptr_ sys, std::ostream& stdout_, std::ostream& stderr_){
            return 1;
        }



        virtual int post_update(system_ptr_ sys, const bool accept){
            return 1;
        }

        virtual int chain_output(const std::string line_prefix, system_ptr_ sys, const bool accept, std::ostream& stdout_, std::ostream& stderr_)
        {
            if (!this->supress_foutput){
                (*output_chain) << line_prefix << "ACCEPTED" << (accept ? 1 : 0) << "\tstep_num:\t" << current_step << "\tLL:\t" << accepted_ll
                    << "\tlog_prior:\t" << accepted_prior;

                for (int ii = 0 ; ii < acc_params.size(); ii++)
                {
                    (*output_chain) << "\t" << acc_params[ii];
                }
                (*output_chain) << std::endl;

                return 1;
            }
            return 1;
        }

        virtual int post_update_output(system_ptr_ sys, const bool accept, std::ostream& stdout_, std::ostream& stderr_){
            //std::cout << "blah1:" << output_freq << "\t" << current_step << "\t" << burn_in << std::endl;
            if ((current_step % output_freq == 0) && current_step >= burn_in){
                //std::cout << "blah2" << std::endl;
                chain_output("", sys, accept, stdout_, stderr_);
            }
            return 1;
        }






        int operator()(system_ptr_ sys, MTRand_shared_ptr rand_ptr, std::ostream& stdout_, std::ostream& stderr_){

            this->rand_gen = rand_ptr;

            //std::cout << "DEBUG: " << std::endl;

            initialise(sys);
            if (initialise_output_files(sys, stdout_, stderr_) != 1){
                return -1;
            }

            chain_output("#\t",sys, 1, stdout_, stderr_);


            for  (current_step = 1; current_step <= steps; current_step++){
                //std::cout << "DEBUG: 1" << std::endl;
                proposal(sys);
                //std::cout << "DEBUG: 2" << std::endl;
                try {
                    //std::cout << "DEBUG: 3" << std::endl;
                    evaluate_ll(sys);
                    //std::cout << "DEBUG: 4" << std::endl;
                }
                catch (const std::exception& e){
                    stderr_ << "# WARNING: exception caught on evaluate_ll():\t" << e.what() << "# " << std::endl;
                    this_sum =  std::numeric_limits< double>::lowest();
                }
                //std::cout << "DEBUG: 5" << std::endl;
                const bool accept = to_accept(sys);
                //std::cout << "DEBUG: 6" << std::endl;
                post_metropolis_eval(sys, accept);
                //std::cout << "DEBUG: 7" << std::endl;
                post_metropolis_eval_output(sys, accept, stdout_, stderr_);
                //std::cout << "DEBUG: 8" << std::endl;

                if (accept){
                    //std::cout << "DEBUG: 9" << std::endl;
                    accepted(sys);
                    //std::cout << "DEBUG: 10" << std::endl;
                    accepted_output(sys, stdout_, stderr_);
                    //std::cout << "DEBUG: 11" << std::endl;
                }
                else {
                    //std::cout << "DEBUG: 12" << std::endl;
                    rejected(sys);
                    //std::cout << "DEBUG: 13" << std::endl;
                    rejected_output(sys, stdout_, stderr_);
                    //std::cout << "DEBUG: 14" << std::endl;
                }
                //std::cout << "DEBUG: 15" << std::endl;
                post_update(sys, accept);
                //std::cout << "DEBUG: 16" << std::endl;
                post_update_output(sys, accept, stdout_, stderr_);
                //std::cout << "DEBUG: 17" << std::endl;

            }


            return 1;
        }




    };


    template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
    class adaptive_mcmc_interface: public mcmc_interface<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_>  {

    public:

        bool greedy_start;
        unsigned long adaptation_step;
        unsigned long min_adaptation_step;

        adaptive_mcmc_interface(system_ptr_ sys, ll_functor_ptr_ get_LL_, prior_ptr_ log_prior_, proposal_ptr_ mover_, std::string output_root_,
                       unsigned long steps_, unsigned long burn_in_, unsigned long output_freq_, const double inv_temp_) : mcmc_interface<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_>(
                       sys, get_LL_, log_prior_, mover_, output_root_,
                       steps_, burn_in_, output_freq_, inv_temp_), greedy_start(true), adaptation_step(1), min_adaptation_step(10)  {
            //std::cout << this->burn_in << std::endl;
        }



        virtual int post_update(system_ptr_ sys, const bool accept){
            if (this->greedy_start && this->current_step < this->burn_in){
                if (accept){
                    this->mover->adapt(std::max(this->accepted_count, min_adaptation_step), this->acceptance_probability, this->acc_params, accept);
                    adaptation_step++;
                }
            }
            else {
                this->mover->adapt(std::max(this->accepted_count, min_adaptation_step), this->acceptance_probability, this->acc_params, accept);
                adaptation_step++;
            }
            return 1;
        }


    };






} // namespace PTR
} // namespace MCMC
#endif /* mcmc_algorithms_h */
