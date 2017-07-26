/*
 * minimiser.h
 *
 *  Created on: 8 Feb 2017
 *      Author: jmacdona
 */

#ifndef SRC_MINIMISER_H_
#define SRC_MINIMISER_H_

#include "ceres/ceres.h"
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Dense>





template<typename model_sys_intf_ptr_type, typename  ll_intf_ptr_type, typename  prior_intf_ptr_type, typename proposal_intf_ptr_type>
class minfunctor : public ceres::FirstOrderFunction {
public:
	virtual ~minfunctor(){

	}
	minfunctor(model_sys_intf_ptr_type sys_,
			ll_intf_ptr_type ll_,
			prior_intf_ptr_type prior_,
			proposal_intf_ptr_type prop_) : sys(sys_), ll(ll_), prior(prior_), prop(prop_){

		mapping = prop_->get_mapping();

	}


	virtual bool Evaluate(const double* parameters,
			double* cost,
			double* gradient) const{

		//typename model_sys_intf_ptr_type::element_type::vector_type intial_params = sys->get_joint_para_and_state();
		map_params(parameters);
		cost[0] = evaluate();

		calc_gradient(parameters, gradient);

		return true;
	}

	virtual int NumParameters() const{
		return mapping.size();
	}

	void initialise_params(double* parameters) {
		typename model_sys_intf_ptr_type::element_type::vector_type params = sys->get_joint_para_and_state();
		for (unsigned ii = 0; ii < mapping.size(); ii++){
			parameters[ii] = params[mapping[ii]];
		}
	}


public:

	void map_params(const double* parameters) const{
		typename model_sys_intf_ptr_type::element_type::vector_type params = sys->get_joint_para_and_state();
		for (unsigned ii = 0; ii < mapping.size(); ii++){
			params[mapping[ii]] = parameters[ii];
			//std::cout << ii << "\t" << params[mapping[ii]] << "\n";
		}
		params = prop->force_in_bounds(params);
		/*
		for (unsigned ii = 0; ii < mapping.size(); ii++){
			std::cout << ii << "\t" << params[mapping[ii]] << "\n";
		}
		*/
		sys->set_joint_para_and_state(params);
	}

	double evaluate() const{
		const double ll_val = ll->operator()(sys);
		const double prior_val =  prior->operator()(sys->get_joint_para_and_state());
		return -(ll_val + prior_val);
	}

	void calc_gradient(const double* parameters,
			double* gradient) const{
		double lowparams[mapping.size()];
		double highparams[mapping.size()];
		std::vector<double> lowvals(mapping.size()), highvals(mapping.size());//, exact_h(mapping.size());

		const double near_0 = 0.01;
		const double standard_h = near_0 * std::cbrt(std::nextafter(near_0, INFINITY) - near_0); //const double standard_h = 0.000001;

		double temp = 0;
		for (unsigned ii = 0; ii < mapping.size(); ii++){

			std::copy(parameters, parameters + mapping.size(), lowparams );
			std::copy(parameters, parameters + mapping.size(), highparams );
			const double this_err = std::nextafter(parameters[ii], INFINITY) - parameters[ii];
			const double h = std::fabs(parameters[ii]) < near_0 ? standard_h : (parameters[ii] * std::cbrt(this_err));//;std::sqrt(this_err)) ;

			//std::cout << "x\t" << parameters[ii] << "\t" << h << std::endl;

			temp = lowparams[ii] - h;
			do_nothing(temp);
			double this_neg_h = temp - lowparams[ii];
			lowparams[ii] = lowparams[ii] + this_neg_h;
			map_params(lowparams);
			double lowval = evaluate();

			temp = highparams[ii] + h;
			do_nothing(temp);
			double this_pos_h = temp - highparams[ii];
			highparams[ii] = highparams[ii] + this_pos_h;
			map_params(highparams);
			double highval = evaluate();

			const double exact_h = this_pos_h - this_neg_h;
			const double this_grad = (highval - lowval) / exact_h;
			gradient[ii] = this_grad;
			//std::cout << "GRAD:\t" << ii << "\t" << this_grad << std::endl;

		}
	}
	
	Eigen::MatrixXd calculate_hessian(){
		
		Eigen::MatrixXd rtn_mat(mapping.size(),mapping.size());
		
		typename model_sys_intf_ptr_type::element_type::vector_type intial_params = sys->get_joint_para_and_state();
		const double central_val = evaluate();
		
		double temp = 0;
		const double near_0 = 0.01;
		const double standard_h = near_0 * std::cbrt(std::nextafter(near_0, INFINITY) - near_0);
		for (unsigned x = 0; x < mapping.size(); x++){
			const double x_val = intial_params[mapping[x]];
			const double x_err = std::nextafter(x_val, INFINITY) - x_val;
			const double h = std::fabs(x_val) < near_0 ? standard_h : (x_val * std::cbrt(x_err));
			
			temp = x_val - h;
			do_nothing(temp);
			const double this_neg_h = temp - x_val;
			const double x_val_low = x_val + this_neg_h;
			
			temp = x_val + h;
			do_nothing(temp);
			double this_pos_h = temp - x_val;
			const double x_val_high = x_val + this_pos_h;
			
			const double exact_2h = this_pos_h - this_neg_h;
			
			
			// do double df_dxdx here
			{
				typename model_sys_intf_ptr_type::element_type::vector_type temp_params = intial_params;
				temp_params[mapping[x]] = x_val_high;
				sys->set_joint_para_and_state(temp_params);
				const double f_px = evaluate();
				
				temp_params = intial_params;
				temp_params[mapping[x]] = x_val_low;
				sys->set_joint_para_and_state(temp_params);
				const double f_mx = evaluate();
				
				const double df_dxdx = (f_px - (2*central_val) + f_mx) / std::pow(exact_2h/2.0,2);
				rtn_mat(x,x) = df_dxdx;
				
			}
			
			
			for (unsigned y = x+1; y < mapping.size(); y++){
				const double y_val = intial_params[mapping[y]];
				const double y_err = std::nextafter(y_val, INFINITY) - y_val;
				const double k = std::fabs(y_val) < near_0 ? standard_h : (y_val * std::cbrt(y_err));
				
				temp = y_val - k;
				do_nothing(temp);
				const double this_neg_k = temp - y_val;
				const double y_val_low = y_val + this_neg_k;
				
				temp = y_val + k;
				do_nothing(temp);
				double this_pos_k = temp - y_val;
				const double y_val_high = y_val + this_pos_k;
				
				const double exact_2k = this_pos_k - this_neg_k;
				
				typename model_sys_intf_ptr_type::element_type::vector_type temp_params = intial_params;
				temp_params[mapping[x]] = x_val_high;
				temp_params[mapping[y]] = y_val_high;
				sys->set_joint_para_and_state(temp_params);
				const double f_px_py = evaluate();
				
				temp_params = intial_params;
				temp_params[mapping[x]] = x_val_high;
				temp_params[mapping[y]] = y_val_low;
				sys->set_joint_para_and_state(temp_params);
				const double f_px_my = evaluate();
				
				temp_params = intial_params;
				temp_params[mapping[x]] = x_val_low;
				temp_params[mapping[y]] = y_val_high;
				sys->set_joint_para_and_state(temp_params);
				const double f_mx_py = evaluate();
				
				temp_params = intial_params;
				temp_params[mapping[x]] = x_val_low;
				temp_params[mapping[y]] = y_val_low;
				sys->set_joint_para_and_state(temp_params);
				const double f_mx_my = evaluate();
				
				const double df_dxdy = (f_px_py - f_px_my - f_mx_py + f_mx_my) / (exact_2h * exact_2k);
				rtn_mat(x,y) = df_dxdy;
				rtn_mat(y,x) = df_dxdy;
				
				
				
			}
		}
		
		sys->set_joint_para_and_state(intial_params);
		return rtn_mat;
	}

	model_sys_intf_ptr_type sys;
	ll_intf_ptr_type ll;
	prior_intf_ptr_type prior;
	proposal_intf_ptr_type prop;

	std::vector<unsigned int> mapping;

	void do_nothing(double val) const{}

};

template<typename model_sys_intf_ptr_type, typename  ll_intf_ptr_type, typename  prior_intf_ptr_type, typename proposal_intf_ptr_type>
bool minimise(model_sys_intf_ptr_type sys, ll_intf_ptr_type ll,
		prior_intf_ptr_type prior,
		proposal_intf_ptr_type prop, const unsigned long steps){
	//std::cout << "# USING ceres minimiser!!!" << std::endl;

	minfunctor< model_sys_intf_ptr_type,
	ll_intf_ptr_type,
	prior_intf_ptr_type,
	proposal_intf_ptr_type> *funct = new minfunctor< model_sys_intf_ptr_type,
	ll_intf_ptr_type,
	prior_intf_ptr_type,
	proposal_intf_ptr_type>(sys, ll, prior, prop);

	unsigned int param_size = funct->NumParameters();
	double *parameters = new double[param_size];
	funct->initialise_params(parameters);


	ceres::GradientProblemSolver::Options options;
	options.logging_type = ceres::SILENT;
	options.minimizer_progress_to_stdout = false;//true;
	options.gradient_tolerance = 1e-12;	//this->tol_; // 1e-10 is default
	options.min_line_search_step_size = 1e-9;//0.00001;//1e-12;// 0.00001;//step_size; // 1e-9 is default
	options.max_num_iterations = steps;
	options.line_search_direction_type = ceres::LBFGS;
	options.function_tolerance = 1e-6;//1e-10; // 1e-6 is default
	options.parameter_tolerance = 1e-8;//1e-12;//1e-8 is default

	ceres::GradientProblemSolver::Summary summary;

	//ceres::GradientProblem takes ownership of funct so assume it deletes afterwards
	ceres::GradientProblem problem(funct);
	ceres::Solve(options, problem, parameters, &summary);
	//std::cout << summary.FullReport() << "\n";
	std::cout << "# minimiser:\t" << summary.BriefReport() << "\n";
	//std::cout << "Initial x: " << -1.2 << " y: " << 1.0 << "\n";
	//std::cout << "Final   x: " << parameters[0] << " y: " << parameters[1] << "\n";
	

	/*
	std::cout << "# calculating Hessian" << std::endl;
	Eigen::MatrixXd hessian = funct->calculate_hessian();
	std::cout << hessian << std::endl;
	std::cout << "# hessian determinant\t" << hessian.determinant() << std::endl;
	std::cout << "\n# calculating covariance from Hessian" << std::endl;
	Eigen::MatrixXd covar = hessian.inverse();
	std::cout << covar << "\n" << std::endl;
	std::cout << "# covar determinant\t" << covar.determinant() << std::endl;
	*/


	delete [] parameters;//free(parameters);

	return true;

}







#endif /* SRC_MINIMISER_H_ */
