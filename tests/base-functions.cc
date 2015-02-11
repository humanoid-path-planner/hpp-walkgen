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

#include <boost/format.hpp>
#include <hpp/walkgen/bspline-based.hh>

#define BOOST_TEST_MODULE base_functions
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE (test_hpp_walkgen)

using hpp::walkgen::PiecewisePoly3;
using hpp::walkgen::SplineBased;
using hpp::walkgen::SplineBasedPtr_t;
using hpp::walkgen::Times_t;
using roboptim::trajectory::CubicBSpline;
using roboptim::trajectory::CubicBSplinePtr_t;
using hpp::walkgen::value_type;
using hpp::walkgen::size_type;
using hpp::walkgen::vector2_t;
using hpp::walkgen::vector_t;


BOOST_AUTO_TEST_CASE (base_functions)
{
    // Define six steps
  SplineBased::Steps_t steps;
  typedef SplineBased::Step Step;
  typedef SplineBased::ZmpTrajectory_t ZmpTrajectory_t;
  typedef SplineBased::ZmpTrajectories_t ZmpTrajectories_t;

  steps.push_back (Step (0, -.1));
  steps.push_back (Step (0, .1));
  steps.push_back (Step (0.2, -.1));
  steps.push_back (Step (0.4, .1));
  steps.push_back (Step (0.6, -.1));
  steps.push_back (Step (0.6, .1));

  // Define times
  Times_t times;
  value_type start = 1., end = 3.4;
  times.push_back (start);
  times.push_back (1.3);
  times.push_back (1.7);
  times.push_back (1.8);
  times.push_back (2.2);
  times.push_back (2.3);
  times.push_back (2.7);
  times.push_back (2.8);
  times.push_back (3.2);
  times.push_back (end);

  value_type height = .8;
  SplineBasedPtr_t pg (SplineBased::create (height));
  pg->timeSequence (times);
  pg->stepSequence (steps);
  pg->solve ();

  size_type i = 0;
  for (ZmpTrajectories_t::const_iterator it = pg->zmpBasisFunctions ().begin ();
       it != pg->zmpBasisFunctions ().end (); ++it) {
    const ZmpTrajectory_t& b = *it;
    std::cout << (boost::format
		  ("set term wxt persist title 'Z%1%' %1%") % i)
		  << std::endl
	      << "set xlabel 't'" << std::endl
	      << "set ylabel 'zmp'" << std::endl;
    std::cout << (boost::format
		  ("plot '-' using 1:2 title 'Z%1%' with points\n") % i);
    ++i;
    std::cout << start << "\t" << 0. << std::endl;
    for (ZmpTrajectory_t::const_iterator it1 = b.begin (); it1 != b.end ();
	 ++it1) {
      const PiecewisePoly3& P = *it1;
      value_type t;
      for (size_type i=0; i<7; ++i) {
	t = ((6-i)*P.lower + i*P.upper)/6.;
	if (!std::isnan (P [i])) {
	  std::cout << t << "\t" << P [i] << std::endl;
	}
      }
    }
    std::cout << end << "\t" << 0. << std::endl;
    std::cout << "e" << std::endl;
  }
}
}
