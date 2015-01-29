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

#ifndef HPP_WALKGEN_FWD_HH
# define HPP_WALKGEN_FWD_HH

# include <vector>
# include <Eigen/Core>
# include <hpp/util/pointer.hh>
# include <roboptim/trajectory/fwd.hh>

namespace roboptim {
  namespace trajectory {
    HPP_PREDEF_CLASS (CubicBSpline);
    typedef boost::shared_ptr <CubicBSpline> CubicBSplinePtr_t;
  } // namespace trajectory
} // namespace roboptim

namespace hpp {
  namespace walkgen {
    HPP_PREDEF_CLASS (SplineBased);
    typedef boost::shared_ptr <SplineBased> SplineBasedPtr_t;

    typedef double value_type;
    typedef Eigen::Matrix <value_type, Eigen::Dynamic, 1> vector_t;
    typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
    typedef matrix_t::Index size_type;
    typedef std::vector <value_type> Times_t;
    typedef Eigen::Matrix <value_type, 2, 1> vector2_t;

  } // namespace walkgen
} // namespace hpp

#endif // HPP_WALKGEN_FWD_HH
