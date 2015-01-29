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

#include <hpp/walkgen/bspline-based.hh>

#define BOOST_TEST_MODULE four_steps
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

BOOST_AUTO_TEST_CASE (four_steps)
{
  // Define four steps
  SplineBased::Steps_t steps;
  typedef SplineBased::Step Step;
  steps.push_back (Step (0, -.1));
  steps.push_back (Step (0, .1));
  steps.push_back (Step (0.2, -.1));
  steps.push_back (Step (0.2, .1));

  // Define times
  Times_t times;
  times.push_back (1.);
  times.push_back (1.3);
  times.push_back (1.7);
  times.push_back (1.8);
  times.push_back (2.2);
  times.push_back (2.3);
  times.push_back (2.7);
  times.push_back (2.9);

  value_type height = .8;
  SplineBasedPtr_t pg (SplineBased::create (height));
  pg->timeSequence (times);
  pg->stepSequence (steps);

  CubicBSplinePtr_t comTrajectory = pg->solve ();

  std::cout << "set term wxt persist title 'walk trajectory'" << std::endl
	    << "set xlabel 'x'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' title 'zmp_ref' with circles\n";
  const std::vector <PiecewisePoly3>& zmpRefx = pg->zmpRefx ();
  const std::vector <PiecewisePoly3>& zmpRefy = pg->zmpRefy ();
  BOOST_CHECK (zmpRefx.size () == zmpRefy.size ());
  for (std::size_t i=0; i<zmpRefx.size (); ++i) {
    const PiecewisePoly3& px = zmpRefx [i];
    const PiecewisePoly3& py = zmpRefy [i];
    for (std::size_t j=0; j<7; ++j) {
      std::cout << px [j] << "\t" << py [j] << std::endl;
    }
  }
}
}
