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

#ifndef HPP_WALKGEN_SPLINE_BASED_HH
# define HPP_WALKGEN_SPLINE_BASED_HH

# include<Eigen/StdVector>

# include <roboptim/trajectory/cubic-b-spline.hh>
# include <hpp/walkgen/config.hh>
# include <hpp/walkgen/fwd.hh>

namespace hpp {
  namespace walkgen {

    using roboptim::trajectory::Polynomial3;
    using roboptim::trajectory::CubicBSpline;
    using roboptim::trajectory::CubicBSplinePtr_t;
    typedef roboptim::trajectory::CubicBSpline::polynomials3vectors_t
    polynomials3vectors_t;
    typedef roboptim::trajectory::CubicBSpline::polynomials3vector_t
    polynomials3vector_t;

    typedef Eigen::Matrix <value_type, 7, 1> vector7_t;

    /// Polynomial function of degree 3 restricted to an interval
    ///
    /// Represented by the values of the function taken at 7 equally spaced
    /// parameter values containing the boundaries of the interval of
    /// definition.
    class PiecewisePoly3
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      /// Constructor
      ///
      /// \param tmin, tmax, interval of definition
      /// \param polynomial polynomial of degree 3 defining the spline
      ///        basis function.
      PiecewisePoly3 (const value_type& tmin, const value_type& tmax,
		      const Polynomial3& polynomial) : lower (tmin),
						       upper (tmax)
      {
	value_type delta = (upper-lower)/6;
	values_ [0] = polynomial (lower);
	values_ [6] = polynomial (upper);
	value_type t = lower + delta;
	for (unsigned int i=1; i<6; ++i) {
	  values_ [i] = polynomial (t);
	  t += delta;
	}
      }
      /// Constructor for affine function
      ///
      /// \param tmin, tmax, interval of definition
      /// \param v0, v1 values at tmin and tmax respectively.
      PiecewisePoly3 (const value_type& tmin, const value_type& tmax,
		      const value_type& v0, const value_type& v1) :
	lower (tmin), upper (tmax)
      {
	value_type delta_p = (upper-lower)/6;
	value_type delta_v = (v1-v0)/6;
	values_ [0] = v0;
	values_ [6] = v1;
	value_type t = lower + delta_p;
	value_type v = v0 + delta_v;
	for (unsigned int i=1; i<6; ++i) {
	  values_ [i] = v;
	  t += delta_p;
	  v += delta_v;
	}
      }
      /// Empty constructor
      /// set all fields to Nan
      PiecewisePoly3 () : lower (sqrt (-1.)), upper (sqrt (-1.))
      {
	values_.fill (sqrt (-1.));
      }
      /// Get value at equally spaced parameter
      ///
      /// \param index index of the parameter between 0 and 6.
      const value_type& operator[] (const size_type& index) const
      {
	return values_ [index];
      }
      value_type lower;
      value_type upper;
    private:
      vector7_t values_;
    }; // class PiecewisePoly3


    /// Walkg motion generator for humanoid legged robot
    ///
    /// This class computes the reference trajectory of the center of mass
    /// of a humanoid robot, given as input a list of time-stamped steps.
    ///
    /// The trajectory of the center of mass is obtained by optimization for the
    /// linear so called table cart model.
    class HPP_WALKGEN_DLLAPI SplineBased
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      static value_type gravity;
      struct Step {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Step (const value_type& abs, const value_type& ord) :
	  position (abs, ord)
	{
	}
	vector2_t position;
	const value_type& operator[] (size_type index) const
	{
	  return position [index];
	}
      }; // struct Step
      struct BoundaryCondition
      {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	BoundaryCondition (const value_type& time, const vector2_t& pos,
			   const vector2_t& vel) : t (time),
						   position (pos),
						   velocity (vel)
	{
	}
	value_type t;
	vector2_t position;
	vector2_t velocity;
      }; // struct BoundaryCondition
      typedef std::vector <Step, Eigen::aligned_allocator <Step> > Steps_t;
      typedef std::vector <BoundaryCondition,
			   Eigen::aligned_allocator <BoundaryCondition> >
	BoundaryConditions_t;
      /// Create instance and return shared pointer
      ///
      /// \param height height of the center of mass
      static SplineBasedPtr_t create (const value_type& height);

      /// set time sequence
      ///
      /// \param times sequence of times corresponding to support phases.
      void timeSequence (const Times_t& times);

      /// get time sequence
      ///
      /// \return sequence of times corresponding to support phases.
      const Times_t& timeSequence () const
      {
	return tau_;
      }

      /// set sequence of steps
      ///
      void stepSequence (const Steps_t& steps);

      /// get sequence of steps
      const Steps_t& stepSequence () const
      {
	return steps_;
      }

      /// Specify zmpref boundary conditions
      ///
      /// \param init, end initial and end values of zmpref trajectory
      /// \note If not specified, the middle of the first (respectively last)
      ///       steps.
      /// \note setting step sequence discards zmp boundary conditions.
      void zmpRefBoundaryConditions (const vector2_t& init,
				     const vector2_t& end)
      {
	zmpRefInit_ = init;
	zmpRefEnd_ = end;
      }
      /// Add a boundary condition on center of mass trajectory
      ///
      /// \param t time when position of center of mass is constrained,
      /// \param position required position of the center of mass.
      /// \note setting time or step sequence discards boundary conditions.
      void add (const BoundaryCondition& boundaryCondition);

      /// Set initial position and velocity of center of mass
      /// \param position position of the center of mass at start \f$t=\tau_0\f$
      /// \param velocity velocity of the center of mass at start \f$t=\tau_0\f$
      void setInitialComState (const vector2_t& position,
			       const vector2_t& velocity);
      /// Set final position and velocity of center of mass
      /// \param position position of the center of mass at end
      ///        \f$t=\tau_{2p-3}\f$
      /// \param velocity velocity of the center of mass at start
      ///        \f$t=\tau_{2p-3}\f$
      void setEndComState (const vector2_t& position,
			   const vector2_t& velocity);
      
      /// Solve quadratic program and return resulting cubic B spline
      ///
      CubicBSplinePtr_t solve () const;
      /// Getter to representation of zmp ref abscissa
      const std::vector <PiecewisePoly3>& zmpRefx () const
      {
	return zmpRef0_;
      }
      /// Getter to representation of zmp ref abscissa
      const std::vector <PiecewisePoly3>& zmpRefy () const
      {
	return zmpRef1_;
      }
    protected:
      /// Constructor
      ///
      /// \param height height of the center of mass
      SplineBased (const value_type& height);
      /// initialization
      /// Store weak pointer to instance
      void init (const SplineBasedWkPtr_t& self);
    private:
      void defineProblem () const;
      void buildPolynomialVector () const;
      value_type height_;
      SplineBasedWkPtr_t weakPtr_;
      mutable size_type m_;
      Times_t tau_;
      Steps_t steps_;
      BoundaryConditions_t boundaryConditions_;
      vector2_t zmpRefInit_;
      vector2_t zmpRefEnd_;
      mutable CubicBSplinePtr_t comTrajectory_;
      mutable matrix_t H0_;
      mutable vector_t b0_, b1_;
      mutable matrix_t A0_;
      mutable vector_t c0_, c1_;
      mutable std::vector < std::vector <PiecewisePoly3> > Z_;
      mutable std::vector <PiecewisePoly3> zmpRef0_;
      mutable std::vector <PiecewisePoly3> zmpRef1_;
    }; // class SplineBased
  } // namespace walkgen
} // namespace hpp
#endif // HPP_WALKGEN_SPLINE_BASED_HH
