//
// Copyright (c) 2015 CNRS
// Authors: Florent Lamiraux
//
// This file is part of hpp-walkgen
// hpp-walkgen is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-walkgen is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-walkgen  If not, see
// <http://www.gnu.org/licenses/>.

#include <iostream>
#include <sstream>
#include <roboptim/trajectory/cubic-b-spline.hh>
#include <hpp/walkgen/bspline-based.hh>

namespace hpp {
  namespace walkgen {
    value_type SplineBased::gravity = 9.81;
    SplineBasedPtr_t SplineBased::create (const value_type& height)
    {
      SplineBased* ptr = new SplineBased (height);
      SplineBasedPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    SplineBased::SplineBased (const value_type& height) : height_ (height)
    {
    }

    void SplineBased::init (const SplineBasedWkPtr_t& self)
    {
      weakPtr_ = self;
    }

    void SplineBased::timeSequence (const Times_t& times)
    {
      tau_ = times;
      boundaryConditions_.clear ();
    }

    void SplineBased::stepSequence (const Steps_t& steps)
    {
      steps_ = steps;
      boundaryConditions_.clear ();
      std::size_t p = steps.size ();
      m_ = 2*p + 6;
      if (p < 2) {
	throw std::runtime_error
	  ("Step sequence should contain at least 2 steps");
      }
      zmpRefInit_ = .5 * (steps [0].position + steps [1].position);
      zmpRefEnd_ = .5 * (steps [p-2].position + steps [p-1].position);
    }

    void SplineBased::add (const BoundaryCondition& boundaryCondition)
    {
      boundaryConditions_.push_back (boundaryCondition);
    }

    void SplineBased::buildPolynomialVector () const
    {
      value_type un_sur_omega_2 = sqrt (gravity / height_);
      const polynomials3vectors_t& bases = comTrajectory_->basisPolynomials ();
      assert ((size_type)bases.size () == m_ - 4);
      const std::vector< value_type >& knots = comTrajectory_->knotVector ();
      Z_.clear ();
      for (std::size_t iBasis = 0; iBasis < bases.size (); ++iBasis) {
	const polynomials3vector_t& basis = bases [iBasis];
	std::vector <PiecewisePoly3> v (4);
	for (std::size_t iPoly = std::max (0, (int)(3 - iBasis));
	     iPoly < 4 && (size_type)(iBasis + iPoly) < m_ - 4; ++iPoly) {
	  const Polynomial3& poly = basis [iPoly];
	  v [iPoly] = PiecewisePoly3
	    (knots [iBasis + iPoly], knots [iBasis + iPoly + 1],
	     poly - un_sur_omega_2 * poly.derivative <2> ());
	}
	Z_.push_back (v);
      }
    }

    static value_type integral (value_type lower,value_type upper,
				const PiecewisePoly3& P1,
				const PiecewisePoly3& P2)
    {
      assert (lower == P1.lower);
      assert (upper == P1.upper);
      assert (lower == P2.lower);
      assert (upper == P2.upper);
      value_type length = P1.upper - P1.lower;
      value_type res = 41./840 * P1 [0] * P2 [0] +
	9./35 * P1 [1] * P2 [1] +
	9./280 * P1 [2] * P2 [2] +
	34./105 * P1 [3] * P2 [3] +
	9./280 * P1 [4] * P2 [4] +
	9./35 * P1 [5] * P2 [5] +
	41./840 * P1 [6] * P2 [6];
      res *= length;
      return res;      
    }

    void SplineBased::defineProblem () const
    {
      std::size_t p = steps_.size ();
      if (p < 2) {
	std::ostringstream oss;
	oss << "size of step vector (" << p << ") should be at least 2.";
	throw std::runtime_error (oss.str ());
      }
      if (tau_.size () != 2*p) {
	std::ostringstream oss;
	oss << "size of time vector (" << tau_.size ()
	    << ") does not fit " << "size of step vector (" << steps_.size ()
	    << ").";
	throw std::runtime_error (oss.str ());
      }
      // build knot vector
      std::vector <value_type> knots (m_);
      double tau0 = tau_ [0];
      knots [0] = tau0 - 3; knots [1] = tau0 - 2;  knots [2] = tau0 - 2;
      for (std::size_t i=3; (size_type)i < m_-3; ++i) {
	knots [i] = tau_ [i-3];
      }
      knots [m_-3] = tau_ [2*p-1] + 1; knots [m_-2] = tau_ [2*p-1] + 2;
      knots [m_-1] = tau_ [2*p-1] + 3;
      // create spline
      vector_t parameters (2*(m_-4));
      comTrajectory_ = CubicBSplinePtr_t (new CubicBSpline
					  (2, knots, parameters));
      // Resize quadratic program
      H0_.resize (m_-4, m_-4);
      b0_.resize (m_-4);
      b1_.resize (m_-4);

      // Fill H0_
      buildPolynomialVector ();
      H0_.setZero ();
      for (size_type i = 0; i < m_ - 4; ++i) {
	for (size_type k = 0; k < m_ - 4 && i + k < m_ - 4; ++k) {
	  for (size_type j=std::max(0, (size_type)(3-i-k));
	       j < 4 - k && i + k + j < m_ - 4; ++j) {
	    H0_ (i, i+k) += integral (knots [i+k+j], knots [i+k+j+1],
				      Z_ [i][k+j], Z_ [i+k][j]);
	  }
	  H0_ (i+k, i) = H0_ (i, i+k);
	}
      }
      // Fill zmpref
      zmpRef0_.clear ();
      zmpRef1_.clear ();
      zmpRef0_.push_back (PiecewisePoly3 (tau_ [0], tau_ [1], zmpRefInit_ [0],
					  steps_ [1][0]));
      zmpRef1_.push_back (PiecewisePoly3 (tau_ [0], tau_ [1], zmpRefInit_ [1],
					  steps_ [1][1]));
      zmpRef0_.push_back (PiecewisePoly3 (tau_ [1], tau_ [2], steps_ [1][0],
					  steps_ [1][0]));
      zmpRef1_.push_back (PiecewisePoly3 (tau_ [1], tau_ [2], steps_ [1][1],
					  steps_ [1][1]));
      for (size_type i=1; i < (size_type)p-1; ++i) {
	zmpRef0_.push_back (PiecewisePoly3 (tau_ [2*i], tau_ [2*i+1],
					    steps_ [i][0], steps_ [i+1][0]));
	zmpRef1_.push_back (PiecewisePoly3 (tau_ [2*i], tau_ [2*i+1],
					    steps_ [i][1], steps_ [i+1][1]));
	zmpRef0_.push_back (PiecewisePoly3 (tau_ [2*i+1], tau_ [2*i+2],
					    steps_ [i+1][0], steps_ [i+1][0]));
	zmpRef1_.push_back (PiecewisePoly3 (tau_ [2*i+1], tau_ [2*i+2],
					    steps_ [i+1][1], steps_ [i+1][1]));
      }
      zmpRef0_.push_back (PiecewisePoly3 (tau_ [2*p-2], tau_ [2*p-1],
					  steps_ [p-1][0], zmpRefEnd_ [0]));
      zmpRef1_.push_back (PiecewisePoly3 (tau_ [2*p-2], tau_ [2*p-1],
					  steps_ [p-1][1], zmpRefEnd_ [1]));
      // Fill b0_ and b1_
      b0_.setZero (); b1_.setZero ();
      for (size_type i = 0; i < m_ - 4; ++i) {
	for (std::size_t j = std::max (0, 3-i); j < 4 && j < 2*p+2-i; ++j) {
	  b0_ [i]  += integral (knots [i+j], knots [i+j+1],
				Z_ [i][j], zmpRef0_ [i+j-3]);
	  b1_ [i]  += integral (knots [i+j], knots [i+j+1],
				Z_ [i][j], zmpRef1_ [i+j-3]);
	}
      }
    }
    CubicBSplinePtr_t SplineBased::solve () const
    {
      defineProblem ();
      vector_t param0 = H0_.colPivHouseholderQr ().solve (b0_);
      vector_t param1 = H0_.colPivHouseholderQr ().solve (b1_);
      vector_t parameters (2*(m_-4));
      for (size_type i=0; i < m_-4; ++i) {
	parameters [2*i] = param0 [i];
	parameters [2*i+1] = param1 [i];
      }
      comTrajectory_->setParameters (parameters);
      return comTrajectory_;
    }
  } // namespace walkgen
} // namespace hpp
