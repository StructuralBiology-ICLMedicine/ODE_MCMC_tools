//
//  thread_job_distributor.h
//  ode_mcmc_code
//
//  Created by jmacdon on 29/12/2016.
//  Copyright © 2016 James T. MacDonald. All rights reserved.
//

#ifndef thread_job_distributor_h
#define thread_job_distributor_h


#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
//#include <boost/timer/timer.hpp>
#include <vector>
#include <deque>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <chrono>
#include "MersenneTwister.h"


class dummy_test_functor {
	
public:
	void operator()(){
		std::cout << "dummy_test_functor: running my thread" << std::endl;
	}
};


class job_functor_interface {

	class aux {

	private:
		boost::mutex is_running_mutex;
		boost::mutex has_finished_mutex;
		boost::mutex times_work_mutex;

		bool is_running_flag;
		bool has_finished_flag;

		double total_times, last_time;
		double total_workload, last_workload;

		//boost::timer::cpu_times total_times;

	public:

		aux() : is_running_flag(false), has_finished_flag(false), total_times(0), last_time(0), total_workload(0), last_workload(0) {
			//total_times.clear();
		}

		void set_running_flag(bool val){
			boost::mutex::scoped_lock lock1(is_running_mutex);
			is_running_flag = val;
		}

		void set_has_finished_flag(bool val){
			boost::mutex::scoped_lock lock1(has_finished_mutex);
			has_finished_flag = val;
		}

		bool get_running_flag(){
			boost::mutex::scoped_lock lock1(is_running_mutex);
			return is_running_flag;
		}

		bool get_has_finished_flag(){
			boost::mutex::scoped_lock lock1(has_finished_mutex);
			return has_finished_flag;
		}

		void add_runtime(double seconds){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			total_times += seconds;
			last_time = seconds;
			//total_workload += workload;
		}

		void add_workload(double workunits){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			total_workload += workunits;
			last_workload = workunits;
		}

		double get_total_runtime(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return total_times;
		}

		double get_total_workload(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return total_workload;
		}

		double get_avg_time_per_work(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return total_times / total_workload;
		}

		double get_last_runtime(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return last_time;
		}

		double get_last_workload(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return last_workload;
		}

		double get_last_time_per_work(){
			boost::mutex::scoped_lock lock1(times_work_mutex);
			return last_time / last_workload;
		}

	};

protected:

	std::string label;
	boost::shared_ptr<aux> aux_ptr_;
	std::shared_ptr<std::chrono::high_resolution_clock> timer_ptr;
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

	bool starting(){
		//std::cerr << "job_functor_template: HELLO" << std::endl;

		if (aux_ptr_->get_has_finished_flag()){
			std::cerr << "job_functor_interface: ERROR: this job has already been run" << std::endl;
			return false;
		}
		else if (aux_ptr_->get_running_flag()) {
			std::cerr << "job_functor_interface: ERROR: this job is already running" << std::endl;
			return false;
		}
		else if (this->check_dependencies() == false){
			std::cerr << "job_functor_interface: ERROR: this job has failed the dependencies check" << std::endl;
			return false;
		}

		aux_ptr_->set_running_flag(true);
		aux_ptr_->set_has_finished_flag(false);
		timer_ptr = std::shared_ptr<std::chrono::high_resolution_clock>(new std::chrono::high_resolution_clock);
		start_time = timer_ptr->now();
		return true;
	}

	bool finished(){
		//timer_ptr->stop();
		std::chrono::time_point<std::chrono::high_resolution_clock> end_time = timer_ptr->now();
		std::chrono::duration<double, std::ratio<1, 1> > dur_secs = end_time - start_time;
		const double run_time = dur_secs.count();
		//std::cout << "# thread took " << run_time << " seconds" << std::endl;
		//std::cout << "# " << boost::timer::format(timer_ptr->elapsed());
		aux_ptr_->add_runtime(run_time);
		aux_ptr_->set_running_flag(false);
		aux_ptr_->set_has_finished_flag(true);
		return true;
	}

	std::vector<  job_functor_interface  > dependencies;


public:

	job_functor_interface() : label("UNK"), aux_ptr_(new aux()){

	}

	virtual ~job_functor_interface() {
	}

	/*
	job_functor_interface(const job_functor_template<ptr_type> &old){

		aux_ptr_ = old.aux_ptr_;

	}

	job_functor_template<ptr_type>& operator=(const job_functor_template<ptr_type>& old) {

		aux_ptr_ = old.aux_ptr_;

		return *this;
	}
	*/

	virtual bool operator()(){
		//std::cout << "# job_functor_template: running my thread" << std::endl;

		if (this->starting()){
			std::cerr << "job_functor_interface: ERROR: this is just an interface class" << std::endl;
			this->finished();
			return true;
		}
		else {
			std::cerr << "job_functor_interface: ERROR: couldn't start thread job" << std::endl;
			return false;
		}

	}

	bool has_finished() const{
		return aux_ptr_->get_has_finished_flag();
	}

	void clear_dependencies(){
		dependencies.clear();
	}

	void add_dependency(job_functor_interface dep){
		//boost::shared_ptr< job_functor_interface > nptr( new job_functor_interface(dep));
		dependencies.push_back(dep);
	}

	// check for job dependencies - returns true if all dependencies return true for has_finished()
	bool check_dependencies() const {
		for (auto it = dependencies.begin(); it != dependencies.end(); it++){
			if (it->has_finished() == false){
				return false;
			}
		}
		return true;
	}

	void reset_flags(){
		aux_ptr_->set_running_flag(false);
		aux_ptr_->set_has_finished_flag(false);
	}

	void set_label(std::string lab_){
		this->label = lab_;
	}

	std::string get_label(){
		return this->label;
	}

};

template<typename ptr_type>
class job_functor_run_chunk : public job_functor_interface {

public:
	ptr_type ptr_;
	unsigned long num_steps;

	job_functor_run_chunk() : job_functor_interface(), ptr_(), num_steps(0) {

	}

	virtual ~job_functor_run_chunk(){
	}

	virtual bool operator()() override{
		//std::cout << "# job_functor_template: running my thread" << std::endl;

		if (this->starting()){
			ptr_->operator()(num_steps);
			aux_ptr_->add_workload(double(num_steps));
			this->finished();
			std::cout << "# job_functor_run_chunk job (" << this->label << ") took " << aux_ptr_->get_last_runtime() << " seconds for " <<   num_steps << " steps."
					<< " avg seconds/step = " << aux_ptr_->get_avg_time_per_work()  << std::endl;
			return true;
		}
		else {
			std::cerr << "job_functor_run_chunk: ERROR: couldn't start job" << std::endl;
			return false;
		}
	}
};

template<typename ptr_type>
class job_functor_do_exchange : public job_functor_interface {

public:
	ptr_type mod_i;
	ptr_type mod_j;
	std::shared_ptr<MTRand> rand_gen;
	//unsigned long num_steps;

	job_functor_do_exchange() : job_functor_interface(), mod_i(), mod_j(){

	}

	virtual ~job_functor_do_exchange(){
	}

	virtual bool operator()() override{
		//std::cout << "# job_functor_template: running my thread" << std::endl;

		if (this->starting()){


			if (mod_i->get_current_model_mapped_id() == mod_j->get_current_model_mapped_id()){

				//typename system_ptr_::element_type::vector_type
				typename ptr_type::element_type::system_ptr_vector::value_type::element_type::vector_type param_i = mod_i->get_current_model_acc_params();
				//typename system_ptr_::element_type::vector_type
				typename ptr_type::element_type::system_ptr_vector::value_type::element_type::vector_type param_j = mod_j->get_current_model_acc_params();

				/*
					std::cout << "param_i";
					for (int ii = 0 ; ii < param_i.size(); ii++)
					{
						std::cout << "\t" << param_i[ii];
					}
					std::cout << std::endl;

					std::cout << "param_j";
					for (int ii = 0 ; ii < param_j.size(); ii++)
					{
						std::cout << "\t" << param_j[ii];
					}
					std::cout << std::endl;
				 */


				const double beta_i = mod_i->get_inv_temp();
				const double beta_j = mod_j->get_inv_temp();
				const double sum_ll_i_param_i = mod_i->get_current_acc_sum_ll();
				const double sum_ll_j_param_j = mod_j->get_current_acc_sum_ll();
				const double sum_ll_i_param_j= mod_i->get_current_model_summed_ll(param_j);
				const double sum_ll_j_param_i= mod_j->get_current_model_summed_ll(param_i);

				const double sum_neg = (beta_i*sum_ll_i_param_i) + (beta_j*sum_ll_j_param_j);
				const double sum_pos = (beta_j*sum_ll_j_param_i) + (beta_i*sum_ll_i_param_j);
				const double diff = sum_pos - sum_neg;
				double prob_acc = std::exp(diff);
				bool accept = false;
				if (!std::isinf(diff) && !std::isnan(diff)){
					if (diff >= 0){
						prob_acc = 1;
						accept = true;
					}
					else {
						const double rand = this->rand_gen->rand();
						if (rand < prob_acc)
						{
							accept = true;
						}
					}
					if (accept){
						mod_i->set_current_model_acc_params(param_j);
						mod_j->set_current_model_acc_params(param_i);

						std::cout << "# EXCHANGE accepted model_id: " << mod_i->get_current_model_mapped_id()
											<< " old_sum: " << sum_neg
											<< " new_sum: " << sum_pos
											<< " diff: " << diff
											<< " prob_acc: " << prob_acc
											<< std::endl;

					}
					else {

						std::cout << "# EXCHANGE rejected model_id: " << mod_i->get_current_model_mapped_id()
											<< " old_sum: " << sum_neg
											<< " new_sum: " << sum_pos
											<< " diff: " << diff
											<< " prob_acc: " << prob_acc
											<< std::endl;

					}
				}
				aux_ptr_->add_workload(double(1));
			}

			this->finished();
			std::cout << "# job_functor_do_exchange job (" << this->label << ")  took " << aux_ptr_->get_last_runtime() << " seconds for " <<   1 << " exchange."
								<< " avg seconds/exchange = " << aux_ptr_->get_avg_time_per_work()  << std::endl;
			return true;
		}
		else {
			std::cerr << "job_functor_do_exchange: ERROR: couldn't start job" << std::endl;
			return false;
		}
	}
};


template <typename job_functor_type>
class wrapper_functor{
public:
	std::shared_ptr<job_functor_type> ptr_;

	void operator()(){
		ptr_->operator()();
	}

	bool check_dependencies(){
		return ptr_->check_dependencies();
	}

	bool has_finished(){
		return ptr_->has_finished();
	}

};

template<typename job_functor_type>
class thread_job_distributor {
	
public:
	
	//typedef std::vector<job_functor_type> job_functor_type_vector;
	typedef std::deque<std::shared_ptr<job_functor_type> > job_functor_type_queue;
	typedef std::deque<std::shared_ptr<wrapper_functor<job_functor_type> > > job_wrapper_functor_type_queue;
	typedef boost::thread thread_type;
	typedef std::shared_ptr<thread_type> thread_ptr_type;
	typedef std::map<thread_ptr_type, std::shared_ptr<wrapper_functor<job_functor_type> > > thread_ptr_job_map;
	
	unsigned long check_freq_ms;
	
	thread_job_distributor():check_freq_ms(1){
		
	}
	

	void operator()(const unsigned int num_threads, job_functor_type_queue bare_jobs){
		job_wrapper_functor_type_queue jobs;
		for (unsigned int ii = 0; ii < bare_jobs.size(); ii++){
			std::shared_ptr<wrapper_functor<job_functor_type> > job_wrap = std::shared_ptr<wrapper_functor<job_functor_type> >( new wrapper_functor<job_functor_type>());
			job_wrap->ptr_ = bare_jobs[ii];
			jobs.push_back(job_wrap);

			// make sure all running flags are set to true
			//jobs[ii].set_running_flag(true);

			// reset flags
			//jobs[ii]->reset_flags();
		}


		thread_ptr_job_map running_jobs;
		while (jobs.size() != 0){

			if (running_jobs.size() < num_threads && jobs.size() != 0 ){
				// fill up the jobs
				unsigned int num_jobs_to_add = std::min(jobs.size(), num_threads - running_jobs.size());
				unsigned int num_added = 0;

				auto curr_pos = jobs.begin();

				while (num_added < num_jobs_to_add && curr_pos != jobs.end()  && jobs.size() != 0){

					if ((*curr_pos)->check_dependencies()){
						thread_ptr_type new_thr_ptr(new thread_type(*(*curr_pos)));
						running_jobs[new_thr_ptr] = *curr_pos;
						curr_pos = jobs.erase(curr_pos);
						//jobs.pop_front();
						boost::this_thread::sleep( boost::posix_time::milliseconds(1) );
						num_added++;
					}
					else {
						curr_pos++;
					}

				}

				/*
				//std::cout << "#\n# adding: " << num_jobs_to_add << " jobs" << std::endl;
				for (int ii = 0; ii < num_jobs_to_add; ii++){
					thread_ptr_type new_thr_ptr(new thread_type(jobs.front()));
					running_jobs[new_thr_ptr] = jobs.front();
					jobs.pop_front();
					boost::this_thread::sleep( boost::posix_time::milliseconds(1) );
				}
				*/

			}

			boost::this_thread::sleep( boost::posix_time::milliseconds(check_freq_ms) );

			for (auto it = running_jobs.cbegin(); it != running_jobs.cend() /* not hoisted */; /* no increment */)
			{
				if (it->second->has_finished())
				{
					//std::cout << "#\n# joining then removing completed job" << std::endl;
					it->first->join();
					running_jobs.erase(it++);    // or "it = m.erase(it)" since C++11
				}
				else
				{
					++it;
				}
			}



			//std::cout << "tick...";
		}

		// join final remaining threads
		for (auto it = running_jobs.cbegin(); it != running_jobs.cend(); it++)
		{
			it->first->join();
		}


		/*
				std::cout << "D8" << std::endl;

				jobs[0].set_running_flag(true);
				thread_ptr_type ptr(new thread_type(jobs[0]));
				boost::this_thread::sleep( boost::posix_time::milliseconds(100) );
				std::cout << "D9" << std::endl;
				if (jobs[0].has_finished() == false){
					std::cout << "# STILL RUNNING" << std::endl;
				}
				else{
					std::cout << "# WTF: HAS FINISHED" << std::endl;
				}
				ptr->join();
				std::cout << "D10" << std::endl;
				if (jobs[0].has_finished() == true){
					std::cout << "# HAS FINISHED" << std::endl;
				}
		 */
	}

	void old_op(const unsigned int num_threads, job_functor_type_queue jobs){


		for (unsigned int ii = 0; ii < jobs.size(); ii++){
			// make sure all running flags are set to true
			//jobs[ii].set_running_flag(true);

			// reset flags
			jobs[ii].reset_flags();
		}
		

		thread_ptr_job_map running_jobs;
		while (jobs.size() != 0){
			
			if (running_jobs.size() < num_threads && jobs.size() != 0 ){
				// fill up the jobs
				unsigned int num_jobs_to_add = std::min(jobs.size(), num_threads - running_jobs.size());
				//std::cout << "#\n# adding: " << num_jobs_to_add << " jobs" << std::endl;
				for (int ii = 0; ii < num_jobs_to_add; ii++){
					thread_ptr_type new_thr_ptr(new thread_type(jobs.front()));
					running_jobs[new_thr_ptr] = jobs.front();
					jobs.pop_front();
					boost::this_thread::sleep( boost::posix_time::milliseconds(1) );
				}
				
			}
			
			for (auto it = running_jobs.cbegin(); it != running_jobs.cend() /* not hoisted */; /* no increment */)
			{
				if (it->second.has_finished())
				{
					//std::cout << "#\n# joining then removing completed job" << std::endl;
					it->first->join();
					running_jobs.erase(it++);    // or "it = m.erase(it)" since C++11
				}
				else
				{
					++it;
				}
			}
			
			
			boost::this_thread::sleep( boost::posix_time::milliseconds(check_freq_ms) );
			//std::cout << "tick...";
		}
		
		// join final remaining threads
		for (auto it = running_jobs.cbegin(); it != running_jobs.cend(); it++)
		{
			it->first->join();
		}
		
		
		/*
		std::cout << "D8" << std::endl;
		
		jobs[0].set_running_flag(true);
		thread_ptr_type ptr(new thread_type(jobs[0]));
		boost::this_thread::sleep( boost::posix_time::milliseconds(100) );
		std::cout << "D9" << std::endl;
		if (jobs[0].has_finished() == false){
			std::cout << "# STILL RUNNING" << std::endl;
		}
		else{
			std::cout << "# WTF: HAS FINISHED" << std::endl;
		}
		ptr->join();
		std::cout << "D10" << std::endl;
		if (jobs[0].has_finished() == true){
			std::cout << "# HAS FINISHED" << std::endl;
		}
		 */
	}

};


#endif /* thread_job_distributor_h */
