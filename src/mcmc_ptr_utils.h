#ifndef MCMC_PTR_UTILS_H_INCLUDED
#define MCMC_PTR_UTILS_H_INCLUDED

#include <tuple>
#include <memory>

#include "mcmc_interfaces.h"
#include "ode_likelihood_function.h"
#include "MersenneTwister.h"

namespace MCMC{

    typedef ode_model_system_interface model_sys_intf_type;
    typedef prior_functor_interface<MCMC::model_system_interface::vector_type> prior_intf_type;
    typedef proposal_functor_interface<model_sys_intf_type::vector_type, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > proposal_intf_type;
    typedef loglikelihood_functor_interface<model_sys_intf_type> ll_intf_type;

    typedef std::shared_ptr<model_sys_intf_type>  model_sys_intf_ptr_type;
    typedef std::shared_ptr<prior_intf_type>  prior_intf_ptr_type;
    typedef std::shared_ptr<proposal_intf_type>  proposal_intf_ptr_type;
    typedef std::shared_ptr<ll_intf_type>  ll_intf_ptr_type;
	
	
	template<typename system_ptr_, typename ll_functor_ptr_, typename prior_ptr_, typename proposal_ptr_ >
	class functor_collection {
	public:
		unsigned int model_id;
		unsigned int P1;
		unsigned int P2;
		system_ptr_ sys;
		ll_functor_ptr_ ll;
		prior_ptr_ prior;
		proposal_ptr_ proposal;
		unsigned int dims;
		
		functor_collection() :
		model_id(0), P1(0), P2(0), dims(2)
		{
			
		}
		
		functor_collection(const functor_collection& old) :
		model_id(old.model_id), P1(old.P1), P2(old.P2),
		sys(old.sys ? old.sys->clone() : system_ptr_() ), // check pointers are not null before cloning
		ll(old.ll ? old.ll->clone() : ll_functor_ptr_() ),
		prior(old.prior ? old.prior->clone() : prior_ptr_()),
		proposal(old.proposal ? old.proposal->clone() : proposal_ptr_() ), dims(old.dims)
		{
			
		}
		
		std::string get_descr() const{
			std::string rtn_string = "ERROR";
			if (dims == 2){
				rtn_string = "M" + std::to_string(model_id) + "P" + std::to_string(P1) + "P" + std::to_string(P2);
			}
			else if (dims == 1) {
				rtn_string = "M" + std::to_string(model_id) + "P" + std::to_string(P1);
			}
			else {
				//WTF;
			}
			return rtn_string;
		}
		
		
	};

    template<typename system_, typename loglikelihood_,  typename prior_, typename proposal_>
	std::tuple<unsigned int, unsigned int, model_sys_intf_ptr_type, ll_intf_ptr_type, prior_intf_ptr_type, proposal_intf_ptr_type > load_ode_model_set( MTRand::MTRand_shared_ptr rand_gen,
                                                                                                                  const std:: string P1_name,
                                                                                                                  const std:: string P2_name,
                                                                                                                  const std::string init_param_filename,
                                                                                                                  const std::string ll_data_list_file,
                                                                                                                  const std::string parameter_bounds_file,
                                                                                                                  const std::string proposal_mvn_sigma_filename,
                                                                                                                  const std::string proposal_mvn_mapping_filename,
                                                                                                                  const std::string prior_mvn_mean_filename,
                                                                                                                  const std::string prior_mvn_sigma_filename,
                                                                                                                  const std::string prior_mvn_mapping_filename){


        std::shared_ptr<system_> sys(new system_());

        const unsigned int P1 = sys->get_joint_para_and_state_index(P1_name);
        const unsigned int P2 = sys->get_joint_para_and_state_index(P2_name);

        if (P1 == NAME_NOT_FOUND || P2 == NAME_NOT_FOUND){
            std::cerr << "ERROR: one or both parameter labels not found: " << P1_name << " " << P2_name << std::endl;
            throw std::exception();
        }
        else {
            std::cout << "# P1=" << P1_name << "=" << P1
                << " P2=" << P2_name << "="<< P2 << std::endl;
        }



        std::shared_ptr<loglikelihood_ > ll_functor(new loglikelihood_());
        std::cout << "# parsing input data list file" << std::endl;
        bool parseOK = ODE::parse_MVN_list_file<system_>(ll_data_list_file, ll_functor->conds, ll_functor->times);
        if (parseOK == false){
            std::cerr << "ERROR: problem parsing file: " << ll_data_list_file << "\n Depending on the error, you probably need to check your file paths" << std::endl;
            exit(EXIT_FAILURE);
        }

        std::cout << "# number of time steps: " << ll_functor->times.size() << std::endl;




        typename system_::vector_type joint_state = sys->new_joint_para_and_state_type();
        std::cout << "# parsing parameters from file: " << init_param_filename << std::endl;
        ODE::parse_initial_params(init_param_filename, joint_state);
        sys->set_joint_para_and_state(joint_state);


        std::cout << "# reading parameter bounds for move_functor: " << parameter_bounds_file << std::endl;
        proposal_intf_ptr_type move_functor(new  proposal_(sys, rand_gen));
        move_functor->loadBounds(*sys, parameter_bounds_file);

        if (proposal_mvn_sigma_filename.compare("default") != 0 && proposal_mvn_mapping_filename.compare("default") != 0){
            std::cout << "# loading MVN proposal sigma from command line" << std::endl;
            //use_RAM ? ram_move_functor->use_mvn_mover = true : move_functor->use_mvn_mover = true;
            std::cout << "# reading MVN sigma matrix for mover_functor: " <<  proposal_mvn_sigma_filename << std::endl;
            move_functor->loadSigma(proposal_mvn_sigma_filename);
            std::cout << "# loading move_mapping for mover_functor: " <<  proposal_mvn_mapping_filename << std::endl;
            move_functor->load_mapping(*sys, proposal_mvn_mapping_filename);
        }
        else {
            // uninf prior
            std::cout << "# using parameter bounds file to set up initial default MVN sigma matrix for move_functor" << std::endl;
        }


        prior_intf_ptr_type prior_ptr(new ODE::uninf_prior_functor<typename system_::vector_type>());
        if (prior_mvn_sigma_filename.compare("default") != 0 && prior_mvn_mapping_filename.compare("default") != 0 && prior_mvn_mean_filename.compare("default") != 0){
            std::cout << "# parsing MVN prior files" << std::endl;
            Eigen::Matrix<double, Eigen::Dynamic, 1> prior_mean = readDynamicMatrix(prior_mvn_mean_filename);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> prior_sigma = readDynamicMatrix(prior_mvn_sigma_filename);
            std::shared_ptr<ODE::multigauss_prior_functor<model_sys_intf_type> > prior_functor(new ODE::multigauss_prior_functor<model_sys_intf_type>() );
            prior_ptr = prior_functor;
            prior_functor->load_mapping(*sys, prior_mvn_mapping_filename);

            prior_functor->set_covar(prior_sigma);
            prior_functor->set_mean(prior_mean);
            bool priorOK = prior_functor->verify();
            if (priorOK == false){
                throw std::exception();
            }

        }
        else {
            std::cout << "# using default uninformative prior" << std::endl;
        }



        return std::make_tuple(P1, P2, sys, ll_functor, prior_ptr, move_functor);
    }


}
#endif // MCMC_PTR_UTILS_H_INCLUDED
