//
//  mulitvariate_normal_rand.h
//  ode_mcmc_code
//
//  Created by jmacdon on 27/12/2016.
//  Copyright © 2016 James MacDonald. All rights reserved.
//

#ifndef mulitvariate_normal_rand_h
#define mulitvariate_normal_rand_h

#include "MersenneTwister.h"
#include <eigen3/Eigen/Dense>

namespace Eigen {
	
	
	namespace internal {
		struct normal_dist_functor
		{
			MTRand::MTRand_shared_ptr rand_gen;
			
			template<typename Index>
			inline  double operator() (Index, Index = 0) const { return rand_gen->randNorm(); }
		};
		
		
		
	} // end namespace internal
	
class multivariate_normal_rand
{
	
	
	
private:
	
	Matrix<double,Dynamic,Dynamic> _covar;
	Matrix<double,Dynamic,Dynamic> _transform;
	Matrix< double, Dynamic, 1> _mean;
	bool _use_cholesky;
	SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic> > _eigenSolver; // drawback: this creates a useless eigenSolver when using Cholesky decomposition, but it yields access to eigenvalues and vectors
	MTRand::MTRand_shared_ptr rand_gen;
	internal::normal_dist_functor randN; // Gaussian functor
	
public:
	
	
	
	multivariate_normal_rand(const Matrix<double,Dynamic,1>& mean,
							 const Matrix<double,Dynamic,Dynamic>& covar,
							 const bool use_cholesky=false,
							 const MTRand::uint32 &seed=111111)
							 :_use_cholesky(use_cholesky), rand_gen(new MTRand(seed))
	{
		
		randN.rand_gen = this->rand_gen;
		setMean(mean);
		setCovar(covar);
	}
	
	multivariate_normal_rand(const multivariate_normal_rand& old):
	_covar(old._covar), _transform(old._transform), _mean(old._mean), _use_cholesky(old._use_cholesky),
	_eigenSolver(old._eigenSolver), rand_gen(new MTRand(*old.rand_gen)), randN()  {
		randN.rand_gen = this->rand_gen;
	}
	
	multivariate_normal_rand(const Matrix<double,Dynamic,1>& mean,
							 const Matrix<double,Dynamic,Dynamic>& covar,
							 const bool use_cholesky,
							 const MTRand::MTRand_shared_ptr rand_gen_)
	:_use_cholesky(use_cholesky)
	{
		set_rand_gen(rand_gen_);
		setMean(mean);
		setCovar(covar);
	}
	
	void set_rand_gen(MTRand::MTRand_shared_ptr new_rand_gen){
		rand_gen = new_rand_gen;
		randN.rand_gen = this->rand_gen;
	}
	
	
	void setMean(const Matrix<double,Dynamic,1>& mean) {
		_mean = mean;
	}
	void setCovar(const Matrix<double,Dynamic,Dynamic>& covar)
	{
		_covar = covar;
		
		// Assuming that we'll be using this repeatedly,
		// compute the transformation matrix that will
		// be applied to unit-variance independent normals
		
		if (_use_cholesky)
		{
	  Eigen::LLT<Eigen::Matrix<double,Dynamic,Dynamic> > cholSolver(_covar);
	  // We can only use the cholesky decomposition if
	  // the covariance matrix is symmetric, pos-definite.
	  // But a covariance matrix might be pos-semi-definite.
	  // In that case, we'll go to an EigenSolver
	  if (cholSolver.info()==Eigen::Success)
	  {
		  // Use cholesky solver
		  _transform = cholSolver.matrixL();
	  }
	  else
	  {
		  throw std::runtime_error("Failed computing the Cholesky decomposition. Use solver instead");
	  }
		}
		else
		{
	  _eigenSolver = SelfAdjointEigenSolver<Matrix<double,Dynamic,Dynamic> >(_covar);
	  _transform = _eigenSolver.eigenvectors()*_eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
		}
	}
	
	/// Draw nn samples from the gaussian and return them
	/// as columns in a Dynamic by nn matrix
	Matrix<double,Dynamic,-1> samples(int nn)
	{
		return (_transform * Matrix<double,Dynamic,-1>::NullaryExpr(_covar.rows(),nn,randN)).colwise() + _mean;
	}
 
	
};
} // namespace Eigen




#endif /* mulitvariate_normal_rand_h */
