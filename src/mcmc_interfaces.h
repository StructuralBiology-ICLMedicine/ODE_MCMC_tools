#ifndef MCMC_INTERFACES_H_INCLUDED
#define MCMC_INTERFACES_H_INCLUDED

#include <memory>
#include <vector>
#include "common_model_defs.h"

#include <eigen3/Eigen/Dense>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "MersenneTwister.h"


namespace MCMC{


class model_system_interface{

public:
    typedef boost::numeric::ublas::vector< double> vector_type;
    typedef boost::numeric::ublas::matrix< double > state_matrix_type;


    virtual unsigned int get_dims() const = 0;

    virtual vector_type new_joint_para_and_state_type() const = 0;
    virtual std::vector<std::string> get_joint_para_and_state_name_vec() const = 0;
    virtual unsigned int get_joint_para_and_state_index(const std::string name) const = 0;
    virtual vector_type get_joint_para_and_state()  = 0;
    virtual void set_joint_para_and_state(vector_type vec)  = 0;
	
	virtual std::shared_ptr< model_system_interface > clone() = 0;

};
typedef std::shared_ptr<model_system_interface> model_system_interface_shared_ptr;
typedef std::shared_ptr<const model_system_interface> const_model_system_interface_shared_ptr;


class ode_model_system_interface{

public:
    typedef boost::numeric::ublas::vector< double> vector_type;
    typedef boost::numeric::ublas::matrix< double > state_matrix_type;


    virtual unsigned int get_dims() const = 0;
    virtual vector_type new_state_type() const = 0;
    virtual state_matrix_type new_state_matrix_type() const = 0;

    virtual vector_type get_state() const = 0;
    virtual vector_type get_dependent_vars(vector_type x_)  = 0;

    virtual vector_type new_joint_para_and_state_type() const = 0;
    virtual std::vector<std::string> get_joint_para_and_state_name_vec() const = 0;
    virtual unsigned int get_joint_para_and_state_index(const std::string name) const = 0;
    virtual vector_type get_joint_para_and_state()  = 0;
    virtual void set_joint_para_and_state(vector_type vec)  = 0;

    virtual void initial_assignments() = 0;


    virtual void operator()( const vector_type &x_ , vector_type &dx , double t ) = 0;
	
	virtual std::shared_ptr< ode_model_system_interface > clone() = 0;

};
typedef std::shared_ptr<ode_model_system_interface> ode_model_system_interface_shared_ptr;
typedef std::shared_ptr<const ode_model_system_interface> const_ode_model_system_interface_shared_ptr;



template<typename system_ptr_>
class ode_model_system_wrapper{

private:
    system_ptr_ sys_ptr;

public:

    ode_model_system_wrapper(system_ptr_ ptr_) : sys_ptr(ptr_){
    }

    void operator()( const typename system_ptr_::element_type::vector_type &x_ , typename system_ptr_::element_type::vector_type &dx , double t ){
        sys_ptr->operator()(x_, dx, t);
    }

};

template<typename system_>
class loglikelihood_functor_interface {
public:

    virtual double operator()( std::shared_ptr<  system_> sys_ptr ) = 0;
	
	virtual std::shared_ptr< loglikelihood_functor_interface<system_> > clone() = 0;
	


};

template <typename para_type_>
class prior_functor_interface{
public:
	
	virtual std::shared_ptr< prior_functor_interface<para_type_> > clone() = 0;
	
    virtual double operator()(const para_type_&  p) = 0;
};

template <typename vector_type, typename covar_matrix_type>
class proposal_functor_interface {
public:
    virtual vector_type operator()( vector_type para_vec) = 0;
	virtual vector_type uniform_random( vector_type para_vec) = 0;

	virtual bool adapt(const unsigned long step, const double acceptance_prob, vector_type para_vec, const bool was_accepted) = 0;
	virtual bool update_para_covar(const unsigned long step, vector_type para_vec) = 0;

    virtual covar_matrix_type get_covar() const = 0;
    virtual covar_matrix_type get_para_covar() const = 0;
    virtual covar_matrix_type get_S_matrix() const = 0;
    virtual covar_matrix_type get_B_matrix() const = 0;
    virtual typename Eigen::Matrix<double,Eigen::Dynamic,1> get_mean() const = 0;
    virtual bool set_para_mean_vec(vector_type para_vec) = 0;
    virtual Eigen::Matrix<double,Eigen::Dynamic, 1> get_mapped_parameters( vector_type para_vec) const = 0;
    virtual vector_type map_changes_back( vector_type  joint, const Eigen::Matrix<double,Eigen::Dynamic,1> move_ ) const = 0;
    virtual void loadSigma(std::string filename) = 0;
    virtual void setSigma(const covar_matrix_type& sigma_) = 0;
    virtual bool set_up_default_covar(const double scale_ = 1) = 0;
    virtual bool scale_covar(const double scale_) = 0;
    virtual bool load_mapping(const ode_model_system_interface &sys, const std::string filename) = 0;
	//! get moveable parameter indices:
	virtual std::vector<unsigned int> get_mapping() const = 0;
    virtual bool loadBounds(const ode_model_system_interface &sys, std::string filename) = 0;
	virtual bool setBounds(unsigned int joint_para_state_index, const double lower_bound_, const double upper_bound_, const double stddev_) = 0;
    virtual std::ostream&  output_mapping_names(std::ostream& output) const = 0;
    virtual bool check_bounds( vector_type  joint, const bool verbose = false) = 0;
    virtual vector_type force_in_bounds(vector_type  joint, const bool verbose = false) = 0;
    virtual typename Eigen::Matrix<double,Eigen::Dynamic,1> draw_test_sample() = 0;


    virtual double get_lower_bound(unsigned int p) const = 0;
    virtual double get_upper_bound(unsigned int p) const = 0;
	
	virtual void set_rand_gen(MTRand::MTRand_shared_ptr rand_ptr) = 0;

    virtual double get_proposal_log_prob(vector_type para_vec_start, vector_type para_vec_end) const{
        return 0;
    }

    virtual unsigned int get_n_dims() const = 0;
	
	virtual std::shared_ptr<proposal_functor_interface<vector_type, covar_matrix_type> > clone() = 0;

};


}
#endif // MCMC_INTERFACES_H_INCLUDED
