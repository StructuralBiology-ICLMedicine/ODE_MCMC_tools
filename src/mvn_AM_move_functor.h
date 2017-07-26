#ifndef MVN_MOVE_FUNCTOR_INCLUDED
#define MVN_MOVE_FUNCTOR_INCLUDED

#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <limits>

// #include "eigenmvn.h"
#include "multivariate_normal_rand.h"
#include "eigen_matrix_utils.h"

#include "common_variables.h"
#include "mcmc_interfaces.h"





template<typename system_>
class AM_move_parameters_functor : public MCMC::proposal_functor_interface<typename system_::vector_type, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > {
	
private:
	MTRand::MTRand_shared_ptr rand_gen;
	
public:
	
		typedef std::shared_ptr< typename MCMC::proposal_functor_interface<typename system_::vector_type, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > > move_shared_ptr;
	
        typename system_::vector_type stddev, low, high;
        std::vector<std::string> para_name_vec;
        Eigen::Matrix<double,Eigen::Dynamic,1> mean , para_mean_vec;
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar, para_covar;
        Eigen::multivariate_normal_rand norm_cholesk;

        std::vector<unsigned int> mapping;
        bool use_mvn_mover;

        double covar_scale_factor;
	
	double ram_gamma;
	unsigned long  max_adapt_step;


		virtual move_shared_ptr clone() override {
			move_shared_ptr new_ptr(new AM_move_parameters_functor<system_>(*this));
			return new_ptr;
		}
	
	AM_move_parameters_functor(const AM_move_parameters_functor<system_> &old)
	: rand_gen(new MTRand(*(old.rand_gen))),
	stddev(old.stddev), low(old.low), high(old.high), para_name_vec(old.para_name_vec),
	mean(old.mean), para_mean_vec(old.para_mean_vec), covar(old.covar), para_covar(old.para_covar),
	norm_cholesk(old.norm_cholesk), mapping(old.mapping), use_mvn_mover(old.use_mvn_mover),
	covar_scale_factor(old.covar_scale_factor), ram_gamma(old.ram_gamma), max_adapt_step(old.max_adapt_step)
	{
		this->set_rand_gen(rand_gen);
	}

	AM_move_parameters_functor(std::shared_ptr<system_> sys_ptr, MTRand::MTRand_shared_ptr rand_gen_) :
		rand_gen(rand_gen_),
		stddev(sys_ptr->new_joint_para_and_state_type()), low(sys_ptr->new_joint_para_and_state_type()), high(sys_ptr->new_joint_para_and_state_type()),
		para_name_vec(sys_ptr->get_joint_para_and_state_name_vec()),
		mean(Eigen::Matrix<double, 1, 1>::Zero()), para_mean_vec(Eigen::Matrix<double, 1, 1>::Zero()), covar(Eigen::Matrix<double, 1, 1>::Zero()), norm_cholesk(mean, covar, false, rand_gen){

		use_mvn_mover = true;
		covar_scale_factor = 1;

		ram_gamma = 0.5;
		max_adapt_step = 1000;

		// default to allow any parameter value but no movement
		//stddev.fill(0);
		std::fill(stddev.begin(), stddev.end(), 0);
		//low.fill(std::numeric_limits< typename system_::vector_type::value_type>::lowest());
		std::fill(low.begin(), low.end(), std::numeric_limits< typename system_::vector_type::value_type>::lowest());
		//high.fill(std::numeric_limits< typename system_::vector_type::value_type>::max());
		std::fill(high.begin(), high.end(), std::numeric_limits< typename system_::vector_type::value_type>::max());
	}

    virtual double get_lower_bound(unsigned int p) const override{
        return low[p];
    }
    virtual double get_upper_bound(unsigned int p) const override{
        return high[p];
    }


	double get_gaussian_sample(double value, double stddev, double low, double high){
		bool ok = false;
		double new_value = value;
		while (ok == false){
			new_value = rand_gen->randNorm(value, stddev);
			if (new_value >= low && new_value <= high) ok = true;
		}
		return new_value;
	}



        bool check_bounds(typename system_::vector_type  joint, const bool verbose = false) override{
            for  (unsigned int ii = 0; ii < joint.size(); ii++){
                if (joint[ii] < low[ii] ){
					if (verbose) std::cout << "WARNING: initial parameter out of bounds - too low. param_index:\t" << ii << "\tvalue:\t" << joint[ii] << std::endl;
                    return false;
                }
                else if (joint[ii] >= high[ii]){
                    if (verbose) std::cout << "WARNING: initial parameter out of bounds - too high. param_index:\t" << ii << "\tvalue:\t" << joint[ii] << std::endl;
                    return false;
                }
            }
            //cout << "OK" << endl;
            return true;
        }

        virtual typename system_::vector_type force_in_bounds(typename system_::vector_type  joint, const bool verbose = false ) override {
        	for  (unsigned int ii = 0; ii < joint.size(); ii++){
        		if (joint[ii] < low[ii] ){
        			const double new_val = low[ii];
        			if (verbose) std::cout << "# changing out-of-bounds parameter:\t" << ii << "\tfrom\t" << joint[ii] << "\tto\t" << new_val << std::endl;
        			joint[ii] = new_val;
        		}
        		else if (joint[ii] >= high[ii]){
        			const double new_val = std::nextafter(high[ii], low[ii]);
        			if (verbose) std::cout << "# changing out-of-bounds parameter:\t" << ii << "\tfrom\t" << joint[ii] << "\tto\t" << new_val << std::endl;
        			joint[ii] = new_val;//std::nextafter(high[ii], low[ii]);//high[ii] - (std::numeric_limits< double>::min() * 10000); // this may not work well
        		}
        	}
        	return joint;
        }

        typename system_::vector_type operator()(typename system_::vector_type  joint) override{



            typename system_::vector_type  orig_joint = joint;
            //bool bounds_ok = false;

            //while (bounds_ok == false){
                if (use_mvn_mover){
                    joint = move_MultiGauss(orig_joint); // move_one(orig_joint);
                }
                else {
                    joint = move_one(orig_joint);
                }
                //bounds_ok = check_bounds(joint);
            //}


            return joint;


            //return move_one(orig_p);

        }
	
	virtual typename system_::vector_type uniform_random( typename system_::vector_type joint) override{
		for (unsigned int ii = 0; ii < joint.size(); ii++){
			if (stddev[ii] != 0 ){
				joint[ii] = low[ii] + rand_gen->randExc(high[ii]);//get_gaussian_sample(joint[ii], stddev[ii], low[ii], high[ii]);
			}
		}
		
		check_bounds(joint, true);
		
		return joint;
		
	}

        typename system_::vector_type map_changes_back( typename system_::vector_type  joint, const Eigen::Matrix<double,Eigen::Dynamic,1> move_ ) const override{
            for (unsigned int ii = 0; ii < move_.rows(); ii++){
                    joint[mapping[ii]] = move_(ii);
            }
            return joint;
        }
	
		virtual std::vector<unsigned int> get_mapping() const override {
			return mapping;
		}	

        typename system_::vector_type move_all(typename system_::vector_type  joint){
                for (unsigned int ii = 0; ii < joint.size(); ii++){
                        if (stddev[ii] != 0 ){
                                //cout << ii << "\t" << p[ii] << " " << stddev[ii] << " " << low[ii] << " " << high[ii] << endl;
                                joint[ii] = get_gaussian_sample(joint[ii], stddev[ii], low[ii], high[ii]);
                        }
                }
                return joint;
        }

        typename system_::vector_type move_MultiGauss(typename system_::vector_type  joint){

            //bmeg_functor::para_type  p = orig_p;

            const double scale = 1;//0.1;

            Eigen::Matrix<double,Eigen::Dynamic,1> move_ = this->norm_cholesk.samples(1);


            for (unsigned int ii = 0; ii < move_.rows(); ii++){

                    joint[mapping[ii]] += scale * move_(ii);

            }


            return joint;
        }


        typename Eigen::Matrix<double,Eigen::Dynamic,1> draw_test_sample() override{

            Eigen::Matrix<double,Eigen::Dynamic,1> move_ = this->norm_cholesk.samples(1);

            return move_;

        }

        typename system_::vector_type move_one(typename system_::vector_type  joint){

            const double scale = 1;

            bool found_suitable =false;
            int para = 0; //rand_gen.randInt(p.size());
            while (found_suitable == false){
                para = rand_gen->randInt(joint.size()-1);
                if (stddev[para] != 0 ){
                    found_suitable = true;
                }
            }

            joint[para] = this->get_gaussian_sample(joint[para], scale * stddev[para], low[para], high[para]);
            //cout << para << "\t" << joint[para] << endl;

            return joint;
        }

        void loadSigma(std::string filename) override{

            covar = readDynamicMatrix(filename);
            mean = Eigen::MatrixXd::Zero(covar.cols(), 1);
            para_mean_vec = Eigen::MatrixXd::Zero(covar.cols(), 1);

            // input file  is assumed to contain actual proposal covariance (i.e. scaled already)
            covar_scale_factor = (2.38*2.38)/double(covar.cols());
            covar = covar / covar_scale_factor;

            norm_cholesk = Eigen::multivariate_normal_rand(mean, covar_scale_factor * covar, false, rand_gen);




        }

        void setSigma(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& sigma_) override{

        	covar = sigma_;
        	mean = Eigen::MatrixXd::Zero(covar.cols(), 1);
        	para_mean_vec = Eigen::MatrixXd::Zero(covar.cols(), 1);

        	// input file  is assumed to contain actual proposal covariance (i.e. scaled already)
        	covar_scale_factor = (2.38*2.38)/double(covar.cols());
        	covar = covar / covar_scale_factor;

        	norm_cholesk = Eigen::multivariate_normal_rand(mean, covar_scale_factor * covar, false, rand_gen);

        }

    bool load_mapping(const MCMC::ode_model_system_interface& sys, const std::string filename) override
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

                    std::cout << "# Reading AM_move_parameters_functor mapping:\t" << name << "\t" << lineNum << "\t" << joint_index << std::endl;
                    if (joint_index == NAME_NOT_FOUND)
                    {
                        throw std::exception();
                    }
                    mapping.push_back(joint_index);
                    lineNum++;
                }
            }
        }
        covar_scale_factor = (2.38*2.38)/double(mapping.size());
        return true;
    }

    virtual bool setBounds(unsigned int joint_para_state_index, const double lower_bound_, const double upper_bound_, const double stddev_) override{


        if (joint_para_state_index > low.size()){
            std::cerr << "# ERROR: setBounds(): index is out of bounds: " << joint_para_state_index << " > " << low.size() << std::endl;
            return false;
        }

        if (lower_bound_ > upper_bound_){
            std::cerr << "# WARNING: setBounds(): lower bound is higher than the upper bound, this will not work I'm afraid" << std::endl;
            return false;
        }


        stddev[joint_para_state_index] = stddev_;
        low[joint_para_state_index] = lower_bound_;
        high[joint_para_state_index] = upper_bound_;

        return true;
    }

    bool loadBounds(const MCMC::ode_model_system_interface &sys, std::string filename) override{

        std::ifstream input(filename.c_str(), std::ios::in);

        string_vector SplitVec;
        std::string lineStr;
        unsigned long length;// lineNum = 0;
        unsigned int num_nonzero = 0;

        while ( !input.eof() ){
            getline(input, lineStr);


            length = lineStr.length();


            if (length > 0){
                boost::split( SplitVec, lineStr, boost::is_any_of("\t") );
                if (SplitVec.size() >= 4 && lineStr[0] != '#'){
                    const std::string name = SplitVec[0];
                    const double this_stddev = boost::lexical_cast<double>(SplitVec[1]);
                    const double this_lower = boost::lexical_cast<double>(SplitVec[2]);
                    const double this_upper = boost::lexical_cast<double>(SplitVec[3]);

                    const unsigned int joint_index = sys.get_joint_para_and_state_index(name);//system_::get_joint_para_x_y_state_index(name);
                    stddev[joint_index] = this_stddev;
                    low[joint_index] = this_lower;
                    high[joint_index] = this_upper;
                    std::cout << "# Reading AM_move_parameters_functor parameter bounds:\t" << name << "\t" << this_lower << "\t" << this_upper << std::endl;
                    if (this_lower > this_upper){
                            std::cerr << "# WARNING: lower bound is higher than the upper bound, this will not work I'm afraid" << std::endl;

                    }

                    if (this_stddev > 0){
                        num_nonzero++;
                    }

                }
            }
        }

        std::cout << "# AM_move_parameters_functor: setting up default covariance matrix from bounds file, " << num_nonzero << " non-zero standard deviations found" << std::endl;
        covar = Eigen::MatrixXd::Zero(num_nonzero, num_nonzero);
        mean = Eigen::MatrixXd::Zero(covar.cols(), 1);
        para_mean_vec = Eigen::MatrixXd::Zero(covar.cols(), 1);
        mapping.clear();
        unsigned int pos = 0;
        for (unsigned int joint_index = 0; joint_index < stddev.size(); joint_index++){
                if (stddev[joint_index] > 0){
                    mapping.push_back(joint_index);
                    covar(pos,pos) = stddev[joint_index]*stddev[joint_index];
                    pos++;
                }
        }
        covar_scale_factor = (2.38*2.38)/double(mapping.size());
        norm_cholesk = Eigen::multivariate_normal_rand(mean, covar_scale_factor*covar, false, rand_gen);


        return true;
    }

    virtual bool set_up_default_covar(const double scale_ = 1) override{



        unsigned int num_nonzero = 0;
        for (unsigned int joint_index = 0; joint_index < stddev.size(); joint_index++){
            if (stddev[joint_index] > 0){
                num_nonzero++;
            }

        }

        std::cout << "# AM_move_parameters_functor: setting up default covariance matrix from bounds file, " << num_nonzero << " non-zero standard deviations found" << std::endl;

        covar = Eigen::MatrixXd::Zero(num_nonzero, num_nonzero);
        mean = Eigen::MatrixXd::Zero(covar.cols(), 1);
        para_mean_vec = Eigen::MatrixXd::Zero(covar.cols(), 1);
        mapping.clear();
        unsigned int pos = 0;
        for (unsigned int joint_index = 0; joint_index < stddev.size(); joint_index++){
                if (stddev[joint_index] > 0){
                    mapping.push_back(joint_index);
                    covar(pos,pos) = scale_ * stddev[joint_index]*stddev[joint_index];
                    pos++;
                }
        }
        covar_scale_factor = (2.38*2.38)/double(mapping.size());
        norm_cholesk = Eigen::multivariate_normal_rand(mean, covar_scale_factor*covar, false, rand_gen);

        return true;
    }




    Eigen::Matrix<double,Eigen::Dynamic, 1> get_mapped_parameters(typename system_::vector_type para_vec) const override {
        Eigen::Matrix<double,Eigen::Dynamic, 1> rtn_vec = Eigen::MatrixXd::Zero(mapping.size(), 1);
        for (unsigned int ii = 0; ii < mapping.size(); ii++){
            rtn_vec(ii) = para_vec[mapping[ii]];
        }
        return rtn_vec;
    }


    bool set_para_mean_vec(typename system_::vector_type para_vec) override{

        for (unsigned int ii = 0; ii < mapping.size(); ii++){
            para_mean_vec(ii) = para_vec[mapping[ii]];
        }

        return true;
    }


    bool update_para_covar(const unsigned long step,  Eigen::MatrixXd new_para_vec){

    	const double para_gamma = 1;
    	const unsigned long min_step = 100;

    	const double epsilon = std::numeric_limits< double>::min()*1000;
    	Eigen::MatrixXd identity = epsilon * Eigen::MatrixXd::Identity(para_covar.cols(), para_covar.rows());
    	//const double fstep = double(step);
    	const double fstep = std::pow(double(std::min(std::max(step,(unsigned long)min_step), max_adapt_step)), para_gamma);

    	para_mean_vec = para_mean_vec + (1/(fstep)) * (new_para_vec - para_mean_vec);
    	para_covar = para_covar + (1/(fstep+1)) * ( ((new_para_vec - para_mean_vec) * (new_para_vec - para_mean_vec).transpose()) - para_covar + identity);

    	return true;
    }

    virtual bool update_para_covar(const unsigned long step, typename system_::vector_type para_vec) override {

    	//std::cerr << "ERROR!!!!!!!!!!!!: Don't this this mover use RAM mover instead!!!" << std::endl;

    	Eigen::MatrixXd new_para_vec = Eigen::MatrixXd::Zero(para_vec.size(), 1);
    	for (unsigned int ii = 0; ii < mapping.size(); ii++){
    		new_para_vec(ii) = para_vec[mapping[ii]];
    	}

    	return false;//this->update_para_covar(step, new_para_vec);
    }

    bool update_proposal_covar_AM(const unsigned long step, typename system_::vector_type para_vec){

        const double epsilon = std::numeric_limits< double>::min();
        //const double fstep = double(step);
		const double fstep = std::pow(double(std::min(std::max(step,(unsigned long)1), max_adapt_step)), ram_gamma);
        Eigen::MatrixXd new_para_vec = Eigen::MatrixXd::Zero(covar.cols(), 1);
        Eigen::MatrixXd identity = epsilon * Eigen::MatrixXd::Identity(covar.cols(), covar.rows());
        for (unsigned int ii = 0; ii < mapping.size(); ii++){
            new_para_vec(ii) = para_vec[mapping[ii]];
        }
        // update mean
        para_mean_vec = para_mean_vec + (1/(fstep+1)) * (new_para_vec - para_mean_vec);

        //update covar
        covar = covar + (1/(fstep+1)) * ( ((new_para_vec - para_mean_vec) * (new_para_vec - para_mean_vec).transpose()) - covar + identity);

        norm_cholesk.setCovar(covar_scale_factor * covar);

        return true;
    }

    bool adapt(const unsigned long step, const double acceptance_prob, typename system_::vector_type para_vec, const bool was_accepted) override{
        return this->update_proposal_covar_AM(step, para_vec);
    }


    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> get_covar()const override{
        // return proposal covariance (i.e. scaled covariance)
        return covar_scale_factor * covar;
    }

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> get_para_covar()const override{
    	Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B =  get_B_matrix();
    	return B * B.transpose();
    }

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> get_S_matrix()const override{
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > _eigenSolver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >(covar);
        Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> S = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
        return S;
    }

    // TODO update this to para_covar
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> get_B_matrix()const override{
    	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > _eigenSolver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> >(covar);
    	Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> B = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
    	return B;
    }

    Eigen::Matrix<double,Eigen::Dynamic, 1> get_mean() const override{
        return para_mean_vec;
    }


    virtual bool scale_covar(const double scale_factor_) override{
        covar = scale_factor_ * covar;
        norm_cholesk.setCovar(covar_scale_factor * covar);

        return true;
    }

    /*
    virtual bool scale_covar(const double scale_ ) override {
    	covar = scale_ * covar;
    }*/

    std::ostream&  output_mapping_names(std::ostream& output) const override{
        //const std::vector<std::string> name_vec = sys.get_joint_para_and_state_name_vec();
        std::vector<unsigned int>::const_iterator iter;
        for (iter = mapping.begin(); iter != mapping.end(); iter++){
                output << para_name_vec[*iter] << std::endl;
        }
        return output;
    }
	
	virtual void set_rand_gen(MTRand::MTRand_shared_ptr rand_ptr) override{
		rand_gen = rand_ptr;
		norm_cholesk.set_rand_gen(rand_gen);
	}

	virtual unsigned int get_n_dims() const override {
		return mapping.size();
	}



};



template <typename system_>
class init_conds_functor
{
public:
    typename system_::state_type operator()(system_ sys)
    {
        //typename system_::state_type init_conds = { 0.0 , 0.0 , sys.p[9], 0.0 };
        return sys.x;
    }
};



#endif // MVN_MOVE_FUNCTOR_INCLUDED
