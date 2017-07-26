#ifndef MODEL_NUMERICAL_JACOBIAN_H_INCLUDED
#define MODEL_NUMERICAL_JACOBIAN_H_INCLUDED

#include <string>
#include <cmath>


template<typename system_>
class model_jacobian
{

private:
    void do_nothing(double val){
    }

public:

    system_ sys;

    model_jacobian(system_ sys_) : sys(sys_)
    {
    }

    void operator()(const typename system_::state_type &x_ , typename system_::state_matrix_type &jacobi, double t, typename system_::state_type &dfdt){


        const double near_0 = 0.01;
        const double standard_h = near_0 * std::cbrt(std::nextafter(near_0, INFINITY) - near_0); //0.000001;
        jacobi = sys.new_state_matrix_type();

        //! dfdt is assumed to be zero at the moment but should calculate it
        dfdt = sys.new_state_type();               // zeros

        typename system_::state_matrix_type low = sys.new_state_matrix_type();
        typename system_::state_matrix_type high = sys.new_state_matrix_type();
        typename system_::state_matrix_type exact_h = sys.new_state_matrix_type();

        double temp = 0;

        for (int ii = 0; ii < x_.size(); ii++ ){

        	const double this_err = std::nextafter(x_[ii], INFINITY) - x_[ii];
        	const double h = std::fabs(x_[ii]) < near_0 ? standard_h : (x_[ii] * std::cbrt(this_err));

            typename system_::state_type this_state_low = x_;
            typename system_::state_type dxdt_low = sys.new_state_type();
            typename system_::state_type this_state_high = x_;
            typename system_::state_type dxdt_high = sys.new_state_type();

            temp = this_state_low[ii] - h;
            do_nothing(temp);
            double this_neg_h = temp - this_state_low[ii];
            this_state_low[ii] = this_state_low[ii] + this_neg_h;

            temp = this_state_high[ii] + h;
            do_nothing(temp);
            double this_pos_h = temp - this_state_high[ii];
            this_state_high[ii] = this_state_high[ii] + this_pos_h;

            sys(this_state_low, dxdt_low, t  );
            sys(this_state_high, dxdt_high, t  );
            column(low, ii) = dxdt_low;
            column(high, ii) = dxdt_high;
            column(exact_h, ii) = typename system_::state_type( x_.size(),this_pos_h - this_neg_h);
        }

        //jacobi = (high - low) / (2.0 * h);
        jacobi = element_div((high - low), exact_h);


        /*
        std::cout << "OUTPUTTING JACOBIAN\n";
        for (int ii = 0; ii < x_.size(); ii++ ){
            for (int jj = 0; jj < x_.size(); jj++ ){
                std::cout << jacobi(ii,jj) << "\t";
            }
            std::cout << std::endl;
        }
        */



    }

};


template<typename system_ptr_>
class model_jacobian_wrapper
{

private:
    void do_nothing(double val){
    }

public:

    system_ptr_ sys_ptr;

    model_jacobian_wrapper(system_ptr_ sys_) : sys_ptr(sys_)
    {
    }

    void operator()(const typename system_ptr_::element_type::vector_type &x_ , typename system_ptr_::element_type::state_matrix_type &jacobi, double t, typename system_ptr_::element_type::vector_type &dfdt){
        //const double h = 0.00001;
        const double near_0 = 0.01;
        const double standard_h = near_0 * std::cbrt(std::nextafter(near_0, INFINITY) - near_0); //0.000001;
        jacobi = sys_ptr->new_state_matrix_type();

        //! dfdt is assumed to be zero at the moment but should calculate it
        dfdt = sys_ptr->new_state_type();               // zeros

        typename system_ptr_::element_type::state_matrix_type low = sys_ptr->new_state_matrix_type();
        typename system_ptr_::element_type::state_matrix_type high = sys_ptr->new_state_matrix_type();
        typename system_ptr_::element_type::state_matrix_type exact_h = sys_ptr->new_state_matrix_type();

        double temp = 0;

        for (int ii = 0; ii < x_.size(); ii++ ){

        	const double this_err = std::nextafter(x_[ii], INFINITY) - x_[ii];
        	const double h = std::fabs(x_[ii]) < near_0 ? standard_h : (x_[ii] * std::cbrt(this_err));

            typename system_ptr_::element_type::vector_type this_state_low = x_;
            typename system_ptr_::element_type::vector_type dxdt_low = sys_ptr->new_state_type();
            typename system_ptr_::element_type::vector_type this_state_high = x_;
            typename system_ptr_::element_type::vector_type dxdt_high = sys_ptr->new_state_type();

            temp = this_state_low[ii] - h;
            do_nothing(temp);
            double this_neg_h = temp - this_state_low[ii];
            this_state_low[ii] = this_state_low[ii] + this_neg_h;

            temp = this_state_high[ii] + h;
            do_nothing(temp);
            double this_pos_h = temp - this_state_high[ii];
            this_state_high[ii] = this_state_high[ii] + this_pos_h;

            sys_ptr->operator()(this_state_low, dxdt_low, t  );
            sys_ptr->operator()(this_state_high, dxdt_high, t  );

            /*
            std::cout << "DEBUG:this_state_high:";
            for (int blah_ = 0; blah_ < this_state_high.size(); blah_++){
            	std::cout << "\t" << this_state_high[blah_];
            }
            std::cout << std::endl;
            */

            column(low, ii) = dxdt_low;
            column(high, ii) = dxdt_high;
            column(exact_h, ii) = typename system_ptr_::element_type::vector_type( x_.size(),this_pos_h - this_neg_h);
        }

        //jacobi = (high - low) / (2.0 * h);
        jacobi = element_div((high - low), exact_h);


        /*
        std::cout << "OUTPUTTING JACOBIAN\n";
        for (int ii = 0; ii < x_.size(); ii++ ){
            for (int jj = 0; jj < x_.size(); jj++ ){
                std::cout << jacobi(ii,jj) << "\t";
            }
            std::cout << std::endl;
        }
        */



    }

};


#endif // MODEL_NUMERICAL_JACOBIAN_H_INCLUDED
