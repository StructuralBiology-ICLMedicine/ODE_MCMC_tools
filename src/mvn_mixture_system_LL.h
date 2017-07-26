#ifndef MVN_MIXTURE_SYSTEM_LL_H_INCLUDED
#define MVN_MIXTURE_SYSTEM_LL_H_INCLUDED


#include <string>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "common_model_defs.h"
#include <limits>

#include "multivariate_normal.h"
#include "mcmc_interfaces.h"




//template <unsigned int N>
class mvn_mixture_model_functor : public MCMC::model_system_interface {

private:


public:
	
	virtual std::shared_ptr< MCMC::model_system_interface > clone() override{
		std::shared_ptr< MCMC::model_system_interface > new_ptr(new mvn_mixture_model_functor(*this));
		return new_ptr;
	}

    const unsigned int dims;

    unsigned int get_dims() const override{
        return dims;
    }

    //! joint type for use in sampling both parameters and initial conditions
    typedef boost::numeric::ublas::vector< double   > joint_para_and_state_type;
    joint_para_and_state_type new_joint_para_and_state_type() const override{
        return joint_para_and_state_type(dims, 0);
    }



    unsigned int get_joint_para_x_y_state_index(const std::string name) const {

        const unsigned int index =  boost::lexical_cast<unsigned int>(name);

        if (index < dims){
            return index;
        }

        return NAME_NOT_FOUND;
    }

    virtual unsigned int get_joint_para_and_state_index(const std::string name) const override{
        return get_joint_para_x_y_state_index(name);
    }

    std::vector<std::string> get_joint_para_and_state_name_vec() const override{
        std::vector<std::string> name_vec(0);
        for (unsigned int ii = 0; ii<dims; ii++){
            name_vec.push_back(std::to_string(ii));
        }
        return name_vec;
    }

    joint_para_and_state_type params;

    mvn_mixture_model_functor(const unsigned int dims_) : dims(dims_) {
        params = new_joint_para_and_state_type();
    }

    mvn_mixture_model_functor(const unsigned int dims_, joint_para_and_state_type joint_) : dims(dims_), params(joint_) {

    }

    joint_para_and_state_type get_joint_para_and_state()  override{
        return params;
    }


    void set_joint_para_and_state(joint_para_and_state_type joint_) override{
        params = joint_;
    }




};


typedef std::vector< MultivariateNormalPDF > MultivariateNormalPDF_vec;
typedef std::vector<double> double_vec;
//typedef std::vector<unsigned int> uint_vec;

template<typename system_>
class mvn_mixture_loglikelihood_functor : public MCMC::loglikelihood_functor_interface<system_> {
    private:

    public:

        const unsigned int dims;

        MultivariateNormalPDF_vec mvns;
        double_vec weights;

	
	virtual std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > clone() override{
		std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > new_ptr(new mvn_mixture_loglikelihood_functor<system_>(*this));
		return new_ptr;
	}


    mvn_mixture_loglikelihood_functor(system_ sys) : dims(sys.get_dims()){
		
    }


    double operator()(system_ sys)
    {

        const typename system_::joint_para_and_state_type  params = sys.get_joint_para_and_state();
        Eigen::Matrix<double,Eigen::Dynamic,1> params_vec(dims, 1) ;//
        //Eigen::Matrix<double,dims,1> params_vec;

        for (unsigned int ii = 0; ii< dims; ii++)
        {
            params_vec(ii) = params[ii];
        }

        double sum_prob = std::numeric_limits< double>::min() * 100;//0;
        //typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };


        for (unsigned int ii = 0; ii < mvns.size(); ii++)
        {
            sum_prob += weights[ii] * mvns[ii].getPdf(params_vec);
        }

        return std::log(sum_prob);
    }



    virtual double operator()( std::shared_ptr< system_> sys_ptr ) override{
        const typename system_::joint_para_and_state_type  params = sys_ptr->get_joint_para_and_state();
        Eigen::Matrix<double,Eigen::Dynamic,1> params_vec(dims, 1) ;//
        //Eigen::Matrix<double,dims,1> params_vec;

        for (unsigned int ii = 0; ii< dims; ii++)
        {
            params_vec(ii) = params[ii];
        }

        double sum_prob = std::numeric_limits< double>::min() * 100;//0;
        //typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };


        for (unsigned int ii = 0; ii < mvns.size(); ii++)
        {
            sum_prob += weights[ii] * mvns[ii].getPdf(params_vec);
        }

        return std::log(sum_prob);
    }


};


template<typename system_>
class mvn_mixture_loglikelihood_ptr_functor : public MCMC::loglikelihood_functor_interface<system_> {
    private:

    public:

        const unsigned int dims;

        MultivariateNormalPDF_vec mvns;
        double_vec weights;
	
	
	virtual std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > clone() override{
		std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > new_ptr(new mvn_mixture_loglikelihood_ptr_functor<system_>(*this));
		return new_ptr;
	}



    mvn_mixture_loglikelihood_ptr_functor(std::shared_ptr< system_> sys_ptr) : dims(sys_ptr->get_dims()){
    }



    //virtual double operator()( std::shared_ptr<  system_> sys_ptr ){}

    virtual double operator()( std::shared_ptr< system_> sys_ptr ) override{
        const typename system_::vector_type  params = sys_ptr->get_joint_para_and_state();
        Eigen::Matrix<double,Eigen::Dynamic,1> params_vec(dims, 1) ;//
        //Eigen::Matrix<double,dims,1> params_vec;

        for (unsigned int ii = 0; ii< dims; ii++)
        {
            params_vec(ii) = params[ii];
        }

        double sum_prob = std::numeric_limits< double>::min() * 100;//0;
        //typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };


        for (unsigned int ii = 0; ii < mvns.size(); ii++)
        {
            sum_prob += weights[ii] * mvns[ii].getPdf(params_vec);
        }

        return std::log(sum_prob);
    }


};

template<typename system_>
class flat_loglikelihood_ptr_functor : public MCMC::loglikelihood_functor_interface<system_> {
    private:

    public:

        const unsigned int dims;

        double val;
        //MultivariateNormalPDF_vec mvns;
        //double_vec weights;
	
	virtual std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > clone() override{
		std::shared_ptr< MCMC::loglikelihood_functor_interface<system_> > new_ptr(new flat_loglikelihood_ptr_functor<system_>(*this));
		return new_ptr;
	}



    flat_loglikelihood_ptr_functor(std::shared_ptr< system_> sys_ptr) : dims(sys_ptr->get_dims()){
    }



    //virtual double operator()( std::shared_ptr<  system_> sys_ptr ){}

    virtual double operator()( std::shared_ptr< system_> sys_ptr ) override{
        return val;
    }


};



#endif // MVN_MIXTURE_SYSTEM_LL_H_INCLUDED
