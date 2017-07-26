//
//  mcmc_manager.h
//  ode_mcmc_code
//
//  Created by jmacdon on 29/12/2016.
//  Copyright © 2016 James. All rights reserved.
//

#ifndef mcmc_manager_h
#define mcmc_manager_h

#include "mcmc_prep_utils.h"
#include "mcmc_ptr_algorithms.h"
#include "mcmc_ptr_model_selection_algorithms_v2.h"

#include "thread_job_distributor.h"

#include "MersenneTwister.h"



namespace MCMC{


template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
class mcmc_manager {
	
	
public:
	
	
	
	typedef MCMC::functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_>  func_collection_type;
	typedef std::vector<func_collection_type>  func_collection_type_vec;
	
	
	unsigned long seed;
	std::shared_ptr<MTRand> rand_gen;
	std::string output_root;
	unsigned long steps;
	unsigned long burn_in;
	unsigned long output_freq;
	double inv_temp;





	unsigned long output_chkpoint_freq;
	double jump_move_probability;
	
	
	unsigned long initial_individual_burn_in;
	unsigned int num_threads;
	unsigned long job_steps_chunk;
	unsigned long exchange_freq;
	
	bool random_start;
	bool greedy_start;
	bool do_random_exchanges;
	

	std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > unbiased_obj_col;
	
	std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > all_sim_obj_col;
	
	std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > >  burnin_obj_col;
	
	typedef job_functor_interface base_job_type;

	typedef job_functor_run_chunk< std::shared_ptr < MCMC::PTR::MODEL_SEL::mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > job_type;

	typedef job_functor_do_exchange< std::shared_ptr < MCMC::PTR::MODEL_SEL::mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > exch_job_type;

	
	mcmc_manager(unsigned long seed_,
			std::string output_root_,
			unsigned long steps_,
			unsigned long burn_in_,
			unsigned long output_freq_,
			double inv_temp_,
			unsigned long output_chkpoint_freq_,
			double jump_move_probability_) :
	seed(seed_), rand_gen(new MTRand(seed)), output_root(output_root_), steps(steps_), burn_in(burn_in_),
	output_freq(output_freq_),  inv_temp(inv_temp_),
	output_chkpoint_freq(output_chkpoint_freq_), jump_move_probability(jump_move_probability_), initial_individual_burn_in(100000),
	num_threads(3), job_steps_chunk(1000), exchange_freq(0), random_start(true), greedy_start(false), do_random_exchanges(true){
		
	}
	

	void set_up_pt_system(system_ptr_ sys, ll_functor_ptr_ ll, prior_ptr_ prior, proposal_ptr_ prop,
			std::vector<double> beta_ladder, std::vector<typename system_ptr_::element_type::vector_type> para_vec){
		do_random_exchanges = false;
		unbiased_obj_col.clear();
		burnin_obj_col.clear();
		all_sim_obj_col.clear();

		for (unsigned int ii = 0; ii < beta_ladder.size(); ii++){
			const double this_beta = beta_ladder[ii];

			std::shared_ptr<MTRand> rand_ptr(new MTRand(seed+ii));
			std::shared_ptr<MTRand> rand_ptr2(new MTRand(seed+ii));

			std::vector<unsigned int> model_id_mapping_;
			model_id_mapping_.push_back(0);
			std::vector<system_ptr_> sys_vec_;
			system_ptr_ new_sys = sys->clone();
			new_sys->set_joint_para_and_state(para_vec[ii % para_vec.size()]);
			sys_vec_.push_back(new_sys);
			std::vector<ll_functor_ptr_> get_LL_vec_;
			get_LL_vec_.push_back(ll->clone());
			std::vector<prior_ptr_> log_prior_vec_;
			log_prior_vec_.push_back(prior->clone());
			std::vector<proposal_ptr_> mover_vec_;
			auto this_prop = prop->clone();
			this_prop->set_rand_gen(rand_ptr2);
			this_prop->scale_covar(1.0/this_beta);
			mover_vec_.push_back(this_prop);


			std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > mcmc_functor(
					new MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ >(
							rand_ptr, std::cout, std::cerr,
							model_id_mapping_, sys_vec_, get_LL_vec_, log_prior_vec_, mover_vec_,
							output_root + "_beta_" + std::to_string(ii), steps, burn_in, output_freq, this_beta));

			mcmc_functor->greedy_start = this->greedy_start;
			mcmc_functor->jump_move_probability = 0;


			unbiased_obj_col.push_back(mcmc_functor);
			all_sim_obj_col.push_back(mcmc_functor);
			seed++;

		}
		std::cout << "# mcmc_manager: created " << unbiased_obj_col.size() << " unbiased adaptive MCMC objects " << std::endl;

	}


	
	bool run(){
		
		
		bool setup_ok = true;
		std::cout << "# mcmc_manager: setting up individual MCMC objects " << std::endl;
		for (int ii = 0; ii < burnin_obj_col.size(); ii++){
			setup_ok = burnin_obj_col[ii]->setup() == 1 && setup_ok;
		}
		if (setup_ok == false){
			std::cerr << "ERROR: there was a problem with setting up the MCMC objects" << std::endl;
			return false;
		}
		
		
		//std::cout << "# mcmc_manager: doing initial burn in " << std::endl;


		this->do_initial_burn_in();

		
		std::cout << "# mcmc_manager: setting up MCMC objects " << std::endl;
		setup_ok = true;
		for (int ii = 0; ii < all_sim_obj_col.size(); ii++){
			setup_ok = all_sim_obj_col[ii]->setup() == 1 && setup_ok;
		}
		if (setup_ok == false){
			std::cerr << "ERROR: there was a problem with setting up the MCMC objects" << std::endl;
			return false;
		}
		
		this->main_loop();
		

		
		return true;
	}
	

	bool main_loop(){


		thread_job_distributor<base_job_type> jd;
		std::deque<std::shared_ptr<job_type> > mcmc_jobs;
		std::cout << "# mcmc_manager: entering the main loop of MCMC using thread_job_distributor()" << std::endl;
		for (int ii = 0; ii < all_sim_obj_col.size(); ii++){
			std::shared_ptr<job_type> job = std::shared_ptr<job_type>( new job_type());
			job->ptr_ = all_sim_obj_col[ii];
			job->num_steps = job_steps_chunk;
			mcmc_jobs.push_back(job);
		}

		std::deque<std::shared_ptr<exch_job_type> > exch_job_store = reserve_exchange_jobs(mcmc_jobs.size() / (unsigned int)2);
		std::deque<std::shared_ptr<exch_job_type> > exch_jobs;

		unsigned long steps_taken = 0;
		while (steps_taken < steps){

			std::deque<std::shared_ptr<base_job_type> > curr_jobs;
			for (int ii = 0 ; ii < mcmc_jobs.size(); ii++){
				curr_jobs.push_back(mcmc_jobs[ii]);
				curr_jobs.back()->reset_flags();
			}

			if (exchange_freq != 0){
				bool ok = create_exchange_jobs(mcmc_jobs,
							exch_jobs,
							exch_job_store);
				if (ok){
					for (int ii = 0 ; ii < exch_jobs.size(); ii++){
						curr_jobs.push_back(exch_jobs[ii]);
						curr_jobs.back()->reset_flags();
					}
				}
				else {
					std::cerr << "ERROR: some error occurred in create_exchange_jobs()" << std::endl;
				}
			}

			jd(num_threads, curr_jobs);
			steps_taken += job_steps_chunk;
			/*
			if ((exchange_freq != 0)){
				if (((steps_taken % exchange_freq) == 0) ){
					//std::cout << "# attempting exchanges at step: " << steps_taken << std::endl;
					do_exchanges();
				}
			} //
			*/

		}

		return true;
	}


	std::deque<std::shared_ptr<exch_job_type> > reserve_exchange_jobs(const unsigned int num){
		std::deque<std::shared_ptr<exch_job_type> > rtn_queue;

		for (int ii = 0 ; ii < num; ii++){
			std::shared_ptr<exch_job_type> job = std::shared_ptr<exch_job_type>(new exch_job_type());
			const unsigned long seed =  rand_gen->randInt();
			job->rand_gen = std::shared_ptr<MTRand>(new MTRand(seed));
			for (int jj = 0 ; jj < 1000; jj++){
				// burn in rand_gen
				job->rand_gen->randInt();
			}
			rtn_queue.push_back(job);
		}

		std::cout << "# reserved " << rtn_queue.size() << " exchange jobs" << std::endl;

		return rtn_queue;
	}


	//exch_job_type
	bool create_exchange_jobs(std::deque<std::shared_ptr<job_type> > mcmc_jobs,
			std::deque<std::shared_ptr<exch_job_type> >& exch_jobs,
			std::deque<std::shared_ptr<exch_job_type> >& exch_job_store){

		//std::cout << "DEBUG: create_exchange_jobs: " << exch_jobs.size() << " " << exch_job_store.size() << std::endl;

		// shuffle mcmc_jobs
		if (do_random_exchanges){
			for (unsigned int ii = mcmc_jobs.size() - 1; ii > 0; ii-- ){
				unsigned int jj = rand_gen->randInt(ii);
				std::swap(mcmc_jobs[ii], mcmc_jobs[jj]);
			}
		}

		// move all exch_jobs to exch_job_store
		const unsigned int to_copy = exch_jobs.size();
		for (int ii = 0; ii < to_copy; ii++){
			exch_job_store.push_back( exch_jobs.front() );
			exch_jobs.pop_front();
			//std::cout <<  ii << " DEBUG: create_exchange_jobs: eh " << exch_jobs.size() << " " << exch_job_store.size() << std::endl;
		}


		const unsigned int min_required_exch =  mcmc_jobs.size()  / (unsigned int)2;

		if (exch_job_store.size() < min_required_exch  ){
			std::cerr << "ERROR: you haven't reserved the correct number of exchange jobs, for "
					<< mcmc_jobs.size() << " MCMC jobs, you need " << min_required_exch << " exchange jobs reserved"
					<< ", not " << exch_job_store.size()
					<< " and " << exch_jobs.size() << " left in exch_jobs queue"
					<< std::endl;
			return false;
		}


		const unsigned int start = do_random_exchanges ? 0 : rand_gen->randInt(1);

		for (unsigned int ii = start; ii < mcmc_jobs.size(); ii+=2){
			if (ii+1 < mcmc_jobs.size()){
				//do_pair_exchange(objs[ii], objs[ii+1]);
				std::shared_ptr<exch_job_type> job = exch_job_store.front();
				exch_job_store.pop_front();
				job->mod_i = mcmc_jobs[ii]->ptr_;
				job->mod_j = mcmc_jobs[ii+1]->ptr_;
				job->clear_dependencies();
				job->add_dependency(*mcmc_jobs[ii]);
				job->add_dependency(*mcmc_jobs[ii+1]);
				exch_jobs.push_back(job);
			}
		}

		//std::cout << "DEBUG: create_exchange_jobs: end: " << exch_jobs.size() << " " << exch_job_store.size() << std::endl;

		return true;
	}



	bool do_initial_burn_in(){
		
		thread_job_distributor<job_type> jd;
		std::deque<std::shared_ptr<job_type> > jobs;
		std::cout << "# mcmc_manager: initial burning in of individual MCMC objects using thread_job_distributor()" << std::endl;
		for (int ii = 0; ii < burnin_obj_col.size(); ii++){
			//job_type job;
			std::shared_ptr<job_type> job = std::shared_ptr<job_type>( new job_type());
			job->ptr_ = burnin_obj_col[ii];
			job->num_steps = initial_individual_burn_in;
			jobs.push_back(job);
		}
		jd(num_threads, jobs);
		return true;
	}
	
	
};
	
}

#endif /* mcmc_manager_h */
