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
#include <hpp/util/debug.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/walkgen/bspline-based.hh>
#include <hpp/walkgen/foot-trajectory.hh>
#include <hpp/constraints/macros.hh>

namespace hpp {
  namespace walkgen {
    value_type SplineBased::gravity = 9.81;
    size_type SplineBased::l = 3;

    SplineBasedPtr_t SplineBased::create (const value_type& height)
    {
      SplineBased* ptr = new SplineBased (height);
      SplineBasedPtr_t shPtr (ptr);
      ptr->init (shPtr);
      return shPtr;
    }

    SplineBased::SplineBased (const value_type& height) :
      height_ (height), defaultStepHeight_ (0),
      leftFoot_ (createFootDevice ()),
      rightFoot_ (createFootDevice ())
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

    void SplineBased::footPrintSequence (const FootPrints_t& footPrints)
    {
      footPrints_ = footPrints;
      boundaryConditions_.clear ();
      std::size_t p = footPrints.size ();
      m_ = (2*p - 3)*l + 7;
      if (p < 3) {
	throw std::runtime_error
	  ("FootPrint sequence should contain at least 3 foot prints");
      }
      if (tau_.size () != 2*p-2) {
	std::ostringstream oss;
	oss << "size of time vector (" << tau_.size ()
	    << ") does not fit " << "size of step vector ("
	    << footPrints_.size ()
	    << ").";
	throw std::runtime_error (oss.str ());
      }
      // build knot vector
      vector_t knots (m_);
      double tau0 = tau_ [0];
      knots [0] = tau0 - 3.; knots [1] = tau0 - 2.;  knots [2] = tau0 - 1.;
      for (size_type k=0; k < (size_type)(2*p-3); ++k) {
	for (size_type j=0; j < l; ++j) {
	  knots [3+k*l+j] = ((value_type)(l-j))/(value_type)l * tau_ [k] +
	    ((value_type)j)/(value_type)l * tau_ [k+1];
	}
      }
      knots [m_-4] = tau_ [2*p-3];
      knots [m_-3] = knots [m_-4] + 1.; knots [m_-2] = knots [m_-4] + 2.;
      knots [m_-1] = knots [m_-4] + 3.;
      // create spline
      vector_t parameters (2*(m_-4));
      comTrajectory_ = CubicBSplinePtr_t (new CubicBSpline
					  (2, knots, parameters));
      // Resize quadratic program
      H0_.resize (m_-4, m_-4);
      b0_.resize (m_-4);
      b1_.resize (m_-4);
      A0_.resize (0, m_-4);
      c0_.resize (0);
      c1_.resize (0);
      // Initialize initial and final values of zmpRef.
      zmpRefInit_ = .5 * (footPrints [0].position + footPrints [1].position);
      zmpRefEnd_ = .5 * (footPrints [p-2].position + footPrints [p-1].position);
      // Set step height to default value
      stepHeights_ = StepHeights_t (p-2, defaultStepHeight_);
    }

    void SplineBased::add (const BoundaryCondition& boundaryCondition)
    {
      throw std::runtime_error ("Not implemented yet");
      boundaryConditions_.push_back (boundaryCondition);
    }

    void SplineBased::setInitialComState (const vector2_t& position,
					  const vector2_t& velocity)
    {
      // Fill right hand side of constraint
      c0_.conservativeResize (c0_.rows () + 2);
      c1_.conservativeResize (c1_.rows () + 2);
      c0_.bottomRows <2> () [0] = position [0];
      c0_.bottomRows <2> () [1] = velocity [0];
      c1_.bottomRows <2> () [0] = position [1];
      c1_.bottomRows <2> () [1] = velocity [1];
      // Fill matrix of constraint
      value_type t3 = tau_ [0];
      A0_.conservativeResize (A0_.rows () + 2, m_-4);
      A0_.bottomRows <2> ().setZero ();
      const polynomials3vectors_t& bases = comTrajectory_->basisPolynomials ();
      for (size_type j=0; j < 4; ++j) {
	A0_.bottomRows <2> () (0, j) = bases [j][3-j] (t3);
	A0_.bottomRows <2> () (1, j) = bases [j][3-j].derivative (t3, 1);
      }
    }

    void SplineBased::setEndComState (const vector2_t& position,
				      const vector2_t& velocity)
    {
      // Fill right hand side of constraint
      c0_.conservativeResize (c0_.rows () + 2);
      c1_.conservativeResize (c1_.rows () + 2);
      c0_.bottomRows <2> () [0] = position [0];
      c0_.bottomRows <2> () [1] = velocity [0];
      c1_.bottomRows <2> () [0] = position [1];
      c1_.bottomRows <2> () [1] = velocity [1];
      // Fill matrix of constraint
      value_type tm_4 = tau_ [tau_.size () - 1];
      A0_.conservativeResize (A0_.rows () + 2, m_-4);
      A0_.bottomRows <2> ().setZero ();
      const polynomials3vectors_t& bases = comTrajectory_->basisPolynomials ();
      for (size_type j=0; j < 3; ++j) {
	A0_.bottomRightCorner (2, 3) (0, j) = bases [m_-7+j][3-j] (tm_4);
	A0_.bottomRightCorner (2, 3) (1, j) =
	  bases [m_-7+j][3-j].derivative (tm_4, 1);
      }
    }

    void SplineBased::buildPolynomialVector () const
    {
      value_type un_sur_omega_2 = height_ / gravity;
      const polynomials3vectors_t& bases = comTrajectory_->basisPolynomials ();
      assert ((size_type)bases.size () == m_ - 4);
      const vector_t& knots = comTrajectory_->knotVector ();
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

    value_type SplineBased::integral (value_type lower,value_type upper,
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
      std::size_t p = footPrints_.size ();
      const vector_t& knots = comTrajectory_->knotVector ();
      // Fill H0_
      buildPolynomialVector ();
      H0_.setZero ();
      for (size_type i = 0; i < m_ - 4; ++i) {
	for (size_type k = 0; k < m_ - 4 && i + k < m_ - 4; ++k) {
	  for (size_type j=std::max((size_type)0, (size_type)(3-i-k));
	       j < 4 - k && i + k + j < m_ - 4; ++j) {
	    H0_ (i, i+k) += integral (knots [i+k+j], knots [i+k+j+1],
				      Z_ [i][k+j], Z_ [i+k][j]);
	  }
	  H0_ (i+k, i) = H0_ (i, i+k);
	}
      }
      hppDout (info, "H0_ = " << std::endl << H0_);
      // Fill zmpref
      zmpRef0_.clear ();
      zmpRef1_.clear ();
      for (size_type j=0; j < l; ++j) {
	zmpRef0_.push_back
	  (PiecewisePoly3 (knots [3+j], knots [4+j],
			   zmpRefInit_ [0]*(value_type)(l-j)/(value_type)l +
			   footPrints_ [1][0]*(value_type)j/(value_type)l,
			   zmpRefInit_ [0]*(value_type)(l-j-1)/(value_type)l +
			   footPrints_ [1][0]*(value_type)(j+1)/(value_type)l));
	zmpRef1_.push_back
	  (PiecewisePoly3 (knots [3+j], knots [4+j],
			   zmpRefInit_ [1]*(value_type)(l-j)/(value_type)l +
			   footPrints_ [1][1]*(value_type)j/(value_type)l,
			   zmpRefInit_ [1]*(value_type)(l-j-1)/(value_type)l +
			   footPrints_ [1][1]*(value_type)(j+1)/(value_type)l));
      }
      for (size_type j=0; j < l; ++j) {
	zmpRef0_.push_back (PiecewisePoly3 (knots [3+l+j], knots [4+l+j],
					    footPrints_ [1][0],
					    footPrints_ [1][0]));
	zmpRef1_.push_back (PiecewisePoly3 (knots [3+l+j], knots [4+l+j],
					    footPrints_ [1][1],
					    footPrints_ [1][1]));
      }
      for (size_type i=1; i < (size_type)p-2; ++i) {
	for (size_type j=0; j < l; ++j) {
	  zmpRef0_.push_back
	    (PiecewisePoly3
	     (knots [3+(2*i)*l+j], knots [4+(2*i)*l+j],
	      footPrints_ [i][0]*(value_type)(l-j)/(value_type)l +
	      footPrints_ [i+1][0]*(value_type)j/(value_type)l,
	      footPrints_ [i][0]*(value_type)(l-j-1)/(value_type)l +
	      footPrints_ [i+1][0]*(value_type)(j+1)/(value_type)l));
	  zmpRef1_.push_back
	    (PiecewisePoly3
	     (knots [3+(2*i)*l+j], knots [4+(2*i)*l+j],
	      footPrints_ [i][1]*(value_type)(l-j)/(value_type)l +
	      footPrints_ [i+1][1]*(value_type)j/(value_type)l,
	      footPrints_ [i][1]*(value_type)(l-j-1)/(value_type)l +
	      footPrints_ [i+1][1]*(value_type)(j+1)/(value_type)l));
	}
	for (size_type j=0; j < l; ++j) {
	  zmpRef0_.push_back (PiecewisePoly3
			      (knots [3+(1+2*i)*l+j], knots [4+(1+2*i)*l+j],
			       footPrints_ [i+1][0], footPrints_ [i+1][0]));
	  zmpRef1_.push_back (PiecewisePoly3
			      (knots [3+(1+2*i)*l+j], knots [4+(1+2*i)*l+j],
			       footPrints_ [i+1][1], footPrints_ [i+1][1]));
	}
      }
      for (size_type j=0; j < l; ++j) {
	zmpRef0_.push_back
	  (PiecewisePoly3
	   (knots [3+(2*p-4)*l+j], knots [4+(2*p-4)*l+j],
	    footPrints_ [p-2][0]*(value_type)(l-j)/(value_type)l +
	    zmpRefEnd_ [0]*(value_type)j/(value_type)l,
	    footPrints_ [p-2][0]*(value_type)(l-j-1)/(value_type)l +
	    zmpRefEnd_ [0]*(value_type)(j+1)/(value_type)l));
	zmpRef1_.push_back
	  (PiecewisePoly3
	   (knots [3+(2*p-4)*l+j], knots [4+(2*p-4)*l+j],
	    footPrints_ [p-2][1]*(value_type)(l-j)/(value_type)l +
	    zmpRefEnd_ [1]*(value_type)j/(value_type)l,
	    footPrints_ [p-2][1]*(value_type)(l-j-1)/(value_type)l +
	    zmpRefEnd_ [1]*(value_type)(j+1)/(value_type)l));
      }
      // Fill b0_ and b1_
      b0_.setZero (); b1_.setZero ();
      for (size_type i = 0; i < m_ - 4; ++i) {
	for (size_type j = std::max ((size_type)0, 3-i); j < 4 &&
	       j < m_-4-i; ++j) {
	  b0_ [i]  += integral (knots [i+j], knots [i+j+1],
				Z_ [i][j], zmpRef0_ [i+j-3]);
	  b1_ [i]  += integral (knots [i+j], knots [i+j+1],
				Z_ [i][j], zmpRef1_ [i+j-3]);
	}
      }
      hppDout (info, "b0_=" << std::endl << b0_.transpose ());
      hppDout (info, "b1_=" << std::endl << b1_.transpose ());
    }

    CubicBSplinePtr_t SplineBased::solve () const
    {
      defineProblem ();
      vector_t X0, X1;
      if (A0_.rows () > 0) {
	// constraints are defined
	Eigen::JacobiSVD <matrix_t> svd (A0_, Eigen::ComputeThinU |
					 Eigen::ComputeFullV);
        HPP_DEBUG_SVDCHECK (svd);
	hppDout (info, "A0 singular values: "
		 << svd.singularValues ().transpose ());
	vector_t X_00 = svd.solve (c0_);
	vector_t X_10 = svd.solve (c1_);
	const matrix_t& V = svd.matrixV ();
	size_type rank = svd.rank ();
	matrix_t V0THi (V.rightCols (m_-4-rank).transpose ()*H0_);
	Eigen::JacobiSVD <matrix_t> svd1 (V0THi*V.rightCols (m_-4-rank),
					  Eigen::ComputeThinU |
					  Eigen::ComputeThinV);
        HPP_DEBUG_SVDCHECK (svd1);
	hppDout (info, "V0^T H0 V0 singular values: "
		 << svd1.singularValues ().transpose ());
	vector_t rhs0 (V.rightCols (m_-4-rank).transpose ()*b0_ - V0THi*X_00);
	vector_t rhs1 (V.rightCols (m_-4-rank).transpose ()*b1_ - V0THi*X_10);
	vector_t u0 (svd1.solve (rhs0));
	vector_t u1 (svd1.solve (rhs1));
	X0 = X_00 + V.rightCols (m_-4-rank) * u0;
	X1 = X_10 + V.rightCols (m_-4-rank) * u1;
	// Check optimality
	vector_t grad0 (H0_*X0 - b0_);
	vector_t grad1 (H0_*X1 - b1_);
	hppDout (info, "grad0 (X0)^T V0="
		 << grad0.transpose ()*V.rightCols (m_-4-rank));
	hppDout (info, "grad1 (X1)^T V0="
		 << grad1.transpose ()*V.rightCols (m_-4-rank));
      } else {
	Eigen::JacobiSVD <matrix_t> svd2 (H0_, Eigen::ComputeThinU |
					 Eigen::ComputeThinV);
        HPP_DEBUG_SVDCHECK (svd2);
	hppDout (info, "H0 singular values: "
		 << svd2.singularValues ().transpose ());
	X0 = svd2.solve (b0_);
	X1 = svd2.solve (b1_);
	vector_t grad0 (H0_*X0 - b0_);
	vector_t grad1 (H0_*X1 - b1_);
	hppDout (info, "grad0 (X0)^T=" << grad0.transpose ());
	hppDout (info, "grad1 (X1)^T=" << grad1.transpose ());
      }
      hppDout (info, "X0=" << X0.transpose ());
      hppDout (info, "X1=" << X1.transpose ());
      vector_t parameters (2*(m_-4));
      for (size_type i=0; i < m_-4; ++i) {
	parameters [2*i] = X0 [i];
	parameters [2*i+1] = X1 [i];
      }
      comTrajectory_->setParameters (parameters);
      hppDout (info, "cost = " << cost ());
      // Compute trajectory of the feet
      computeFootTrajectory ();

      return comTrajectory_;
    }

    value_type SplineBased::cost (const vector_t& controlPoints)
    {
      comTrajectory_->setParameters (controlPoints);
      return cost ();
    }

    value_type SplineBased::cost () const
    {
      // Compute cost at minumum
      value_type cost = 0;
      value_type un_sur_omega_2 = height_ / gravity;
      for (std::size_t i=0; i < zmpRef0_.size (); ++i) {
	value_type lower = zmpRef0_ [i].lower;
	value_type upper = zmpRef0_ [i].upper;
	vector7_t values;
	for (std::size_t j = 0; j < 7; ++j) {
	  value_type t = (value_type)(6-j)/6. * lower +
	    (value_type)j/6. * upper;
	  vector_t zmp = (*comTrajectory_) (t) -
	    un_sur_omega_2 * comTrajectory_->derivative (t, 2);
	  values [j] = zmp [0] - zmpRef0_ [i][j];
	}
	PiecewisePoly3 px1 (lower, upper, values);
	PiecewisePoly3 px2 (lower, upper, values);
	cost += integral (lower, upper, px1, px2);

	for (std::size_t j = 0; j < 7; ++j) {
	  value_type t = (value_type)(6-j)/6. * lower +
	    (value_type)j/6. * upper;
	  vector_t zmp = (*comTrajectory_) (t) -
	    un_sur_omega_2 * comTrajectory_->derivative (t, 2);
	  values [j] = zmp [1] - zmpRef1_ [i][j];
	}
	PiecewisePoly3 py1 (lower, upper, values);
	PiecewisePoly3 py2 (lower, upper, values);
	cost += integral (lower, upper, py1, py2);
      }
      return cost;
    }

    void SplineBased::computeFootTrajectory () const
    {
      leftFootTraj_ = PathVector::create (leftFoot_->configSize (),
					  leftFoot_->numberDof ());
      rightFootTraj_ = PathVector::create (rightFoot_->configSize (),
					   rightFoot_->numberDof ());
      // Determine which foot is moving first
      vector2_t delta (footPrints_ [1].position - footPrints_ [0].position);
      vector2_t u (footPrints_ [0].orientation);
      PathVectorPtr_t firstFoot, secondFoot;
      if (u [0]*delta [1] - u [1]*delta [0] > 0) {
	//second foot print is at the left of first foot print
	firstFoot = rightFootTraj_;
	secondFoot = leftFootTraj_;
      } else {
	firstFoot = leftFootTraj_;
	secondFoot = rightFootTraj_;
      }
      // First both feet remain static
      firstFoot->appendPath (SupportFoot::create (footPrints_ [0], 0,
						  tau_ [1] - tau_ [0]));
      secondFoot->appendPath (SupportFoot::create (footPrints_ [1], 0,
						   tau_ [1] - tau_ [0]));
      // Loop
      std::size_t i=0, j=0;
      while (true) {
        // Previous trajectory is a double support trajectory
	++j;
	// Single support
	firstFoot->appendPath (Step::create
			       (footPrints_ [i], footPrints_ [i+2],
				0., stepHeights_ [i],
				tau_ [j+1] - tau_ [j]));
	secondFoot->appendPath (SupportFoot::create
				(footPrints_ [i+1], 0,
				 tau_ [j+1] - tau_ [j]));
	
	++j;
	// Double support
	firstFoot->appendPath (SupportFoot::create
			       (footPrints_ [i+2], 0,
				tau_ [j+1] - tau_ [j]));
	secondFoot->appendPath (SupportFoot::create
				(footPrints_ [i+1], 0,
				 tau_ [j+1] - tau_ [j]));
	++i;
        /// Swap the foot
        firstFoot.swap (secondFoot);

        /// Check if this step was the last one
	if (i+2 == footPrints_.size ()) break;
      }
      // Last, both feet remain static
      // firstFoot->appendPath (SupportFoot::create (footPrints_ [i], 0,
						  // tau_ [j+1] - tau_ [j]));
      // secondFoot->appendPath (SupportFoot::create (footPrints_ [i+1], 0,
						   // tau_ [j+1] - tau_ [j]));
    }

  } // namespace walkgen
} // namespace hpp
