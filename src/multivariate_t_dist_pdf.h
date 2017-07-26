#ifndef MULTIVARIATE_T_DIST_PDF_H_INCLUDED
#define MULTIVARIATE_T_DIST_PDF_H_INCLUDED

#ifndef M_PI
#define M_PI REAL(3.1415926535897932384626433832795029)
#endif

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include<boost/algorithm/string.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/algorithm/string/split.hpp>
#include <vector>
#include  <cmath>

#define MAXBUFSIZE  ((int) 1e6)




class MultivariateStudentsTPDF{
private:

    Eigen::Matrix<double,Eigen::Dynamic,1> mean ;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_inv;

    double det, logdet;
    double normalising_term, log_normalising_term;
	double degf; // degrees of freedom

public:




    MultivariateStudentsTPDF(){
		degf = 3;
    }
	
	
	/*
	void set_degf(double val){
		degf = val;
	}
	*/

    Eigen::Matrix<double,Eigen::Dynamic,1> get_mean() const{
        return mean;
    }

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> get_covar() const{
        return covar;
    }

    void setCovar(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> covar_, const double degf, const bool use_factorization = true){
        covar = covar_;
        const double p = double(covar.rows());
        if (use_factorization == false){
            covar_inv = covar.inverse();
            det = covar.determinant();
			logdet = std::log(det);
        }
        else {
            Eigen::FullPivHouseholderQR< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > qr(covar);
            covar_inv = qr.inverse();
            logdet = qr.logAbsDeterminant();
            det = std::exp(logdet);
        }
        const double lognumerator = std::lgamma((degf+p)/2.0);
        const double logdenominator = std::lgamma(degf/2.0) + ((p/2.0) * std::log(degf)) + ((p/2.0) * std::log(M_PI)) + (logdet/2.0);
        log_normalising_term = lognumerator - logdenominator;
        normalising_term = std::log(log_normalising_term);
		//log_normalising_term = 0.5 * (-std::log(std::pow(2.0 * M_PI, covar.rows())) - logdet);
		//normalising_term = std::exp(log_normalising_term);
    }

    void setMean(Eigen::Matrix<double,Eigen::Dynamic,1> mean_){
        mean = mean_;
    }

    double  getPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){
        /*
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
        const double exp_term =  std::exp((-0.5 * (val_cent.transpose()* covar_inv * val_cent))(0,0));
        const double prob = normalising_term * exp_term;
        */
        return std::exp(this->getLogPdf(val)); //prob;
    }


    double  getLogPdf(Eigen::Matrix<double,Eigen::Dynamic,1> val){
		const double p = double(covar.rows());
        Eigen::Matrix<double,Eigen::Dynamic,1>  val_cent = val - mean;
		const double mat_term =  -((degf+p)/2.0) * std::log(1 + ((1.0/degf)*(val_cent.transpose()* covar_inv * val_cent)(0,0)));
		//const double lognumerator = std::lgamma((degf+p)/2.0);
		//const double logdenominator = std::lgamma(degf/2.0) + ((p/2.0) * std::log(degf)) + ((p/2.0) * std::log(M_PI)) + (logdet/2.0);
		const double log_prob = log_normalising_term + mat_term;
        return log_prob;
    }



};


#endif // MULTIVARIATE_T_DIST_PDF_H_INCLUDED
