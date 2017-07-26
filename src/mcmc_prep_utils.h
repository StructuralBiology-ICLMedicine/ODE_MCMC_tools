//
//  mcmc_prep_utils.h
//  ode_mcmc_code
//
//  Created by jmacdon on 29/12/2016.
//  Copyright © 2016 James T. MacDonald. All rights reserved.
//

#ifndef mcmc_prep_utils_h
#define mcmc_prep_utils_h

#include <limits>

#include "mcmc_ptr_utils.h"
#include "mcmc_ptr_model_selection_algorithms_v2.h"
#include "mcmc_ptr_algorithms.h"
#include "MersenneTwister.h"
#include <map>
#include <algorithm>


namespace MCMC{



	


	template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
	std::vector< std::vector< functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > > equalise_vectors_by_repeat(std::vector<std::vector< functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > > &all){
		
		std::vector< std::vector< functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > > rtn_vec;
		
		std::cout << "DEBUG1" << std::endl;
		
		int max_len = 0;
		for (auto it = all.begin(); it != all.end(); it++ ){
			if (it->size() > max_len){
				max_len = it->size();
			}
		}
		
		std::cout << "DEBUG2" << std::endl;
		
		for (auto it = all.begin(); it != all.end(); it++ ){
			
			std::vector< functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > this_vec = (*it);
			int num_to_insert = max_len - it->size();
			const int orig_size = it->size();
			if (num_to_insert > 0){
				while (num_to_insert != 0){
					num_to_insert = std::min<int>(max_len - this_vec.size(), orig_size);
					if (num_to_insert > 0){
						this_vec.insert(this_vec.end(), this_vec.begin(), this_vec.begin() + num_to_insert );
					}
				}
			}
			rtn_vec.push_back(this_vec);
		}
		
		std::cout << "DEBUG3" << std::endl;
		return rtn_vec;
	}
	

	
	template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
	std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > > generate_unbiased_mcmc_objects( std::vector< functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> >  models, unsigned long seed, std::string output_root, unsigned long steps, unsigned long burn_in, unsigned long output_freq, double inv_temp){
		
		std::vector < std::shared_ptr <  MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> > > rtn_vec;
		
		for (auto it = models.begin(); it != models.end(); it++ ){
			
			// make sure everything is cloned
			functor_collection<system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_> this_model(*it);
			
			std::shared_ptr<MTRand> rand_ptr(new MTRand(seed));
			
			std::vector<unsigned int> model_id_mapping_;
			model_id_mapping_.push_back(this_model.model_id);
			std::vector<system_ptr_> sys_vec_;
			sys_vec_.push_back(this_model.sys);
			std::vector<ll_functor_ptr_> get_LL_vec_;
			get_LL_vec_.push_back(this_model.ll);
			std::vector<prior_ptr_> log_prior_vec_;
			log_prior_vec_.push_back(this_model.prior);
			std::vector<proposal_ptr_> mover_vec_;
			mover_vec_.push_back(this_model.proposal);
			
			
			//std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > single_mcmc_functor( new MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ >(rand_ptr, std::cout, std::cerr, model_id_mapping_, sys_vec_, get_LL_vec_, log_prior_vec_, mover_vec_, output_root + "_unbiased", steps, burn_in,  output_freq , inv_temp));
			
			std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > mcmc_functor( new MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ >(rand_ptr, std::cout, std::cerr, model_id_mapping_, sys_vec_, get_LL_vec_, log_prior_vec_, mover_vec_, output_root + "_unbiased", steps, burn_in, output_freq, inv_temp));
			
			mcmc_functor->greedy_start = false;
			mcmc_functor->jump_move_probability = 0;
			
			
			rtn_vec.push_back(mcmc_functor);
			seed++;
		}
		
		
		return rtn_vec;
	}

	
	template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_ >
	std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > generate_burnin_mcmc_objects(unsigned long seed, std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > >  models, std::string output_root, unsigned long steps, double inv_temp, const bool random_start){
		
		const unsigned long burn_in = 1;
		
		std::vector < std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > rtn_vec;
		
		
		for (unsigned int mod_num = 0; mod_num < models->sys_vec.size(); mod_num++){
			
			std::shared_ptr<MTRand> rand_ptr(new MTRand(seed+mod_num));
			
			if (random_start){
				auto joint = models->sys_vec[mod_num]->get_joint_para_and_state();
				joint = models->mover_vec[mod_num]->uniform_random(joint);
				models->sys_vec[mod_num]->set_joint_para_and_state(joint);
			}
			
			std::vector<unsigned int> model_id_mapping_;
			model_id_mapping_.push_back(models->get_mapping()[mod_num]);
			std::vector<system_ptr_> sys_vec_;
			sys_vec_.push_back(models->sys_vec[mod_num]);
			std::vector<ll_functor_ptr_> get_LL_vec_;
			get_LL_vec_.push_back(models->get_LL_vec[mod_num]);
			std::vector<prior_ptr_> log_prior_vec_;
			log_prior_vec_.push_back(models->log_prior_vec[mod_num]);
			std::vector<proposal_ptr_> mover_vec_;
			mover_vec_.push_back(models->mover_vec[mod_num]);
			
			std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > single_mcmc_functor( new MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ >(rand_ptr, std::cout, std::cerr, model_id_mapping_, sys_vec_, get_LL_vec_, log_prior_vec_, mover_vec_, "BURNIN" + output_root + "_model_" + std::to_string(models->get_mapping()[mod_num]) + "_" + std::to_string(seed+mod_num), steps, burn_in,  99999999 , inv_temp));
			
			single_mcmc_functor->greedy_start = true;//false;
			single_mcmc_functor->supress_foutput = true;
			single_mcmc_functor->jump_move_probability = 0;
			rtn_vec.push_back(single_mcmc_functor);
		}
		
		return rtn_vec;
	}

	
	template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_>
	std::map< unsigned int,  std::vector< std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > >group_by_current_model_id(std::vector< std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > obj_col){
		
		std::map< unsigned int, std::vector< std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > > > rtn_vec;
		
		for (auto it = obj_col.begin(); it != obj_col.end(); it++ ){
			const unsigned int key = (*it)->get_current_model_mapped_id();
			if (rtn_vec.find(key) == rtn_vec.end()){
				rtn_vec[key] = std::vector< std::shared_ptr < MCMC::PTR::MODEL_SEL::adaptive_mcmc_interface_v2< system_ptr_, ll_functor_ptr_, prior_ptr_, proposal_ptr_ > > >();
			}
			rtn_vec[key].push_back(*it);
		}
		
		
		
		return rtn_vec;
	}
	
	
} // namespace MCMC
#endif /* mcmc_prep_utils_h */
