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

#define BOOST_TEST_MODULE six_steps
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

BOOST_AUTO_TEST_CASE (six_steps)
{
  // Define twenty steps
  SplineBased::Steps_t steps;
  typedef SplineBased::Step Step;
  steps.push_back (Step (0, -.1));
  steps.push_back (Step (0, .1));
  steps.push_back (Step (0.2, -.1));
  steps.push_back (Step (0.4, .1));
  steps.push_back (Step (0.6, -.1));
  steps.push_back (Step (0.8, .1));
  steps.push_back (Step (1.0, -.1));
  steps.push_back (Step (1.2, .1));
  steps.push_back (Step (1.4, -.1));
  steps.push_back (Step (1.6, .1));
  steps.push_back (Step (1.8, -.1));
  steps.push_back (Step (2.0, .1));
  steps.push_back (Step (2.2, -.1));
  steps.push_back (Step (2.4, .1));
  steps.push_back (Step (2.6, -.1));
  steps.push_back (Step (2.8, .1));
  steps.push_back (Step (3.0, -.1));
  steps.push_back (Step (3.2, .1));
  steps.push_back (Step (3.4, -.1));
  steps.push_back (Step (3.4, .1));

  // Define times
  Times_t times;
  times.push_back (1.);
  times.push_back (1.3);
  times.push_back (1.7);
  times.push_back (1.8);
  times.push_back (2.2);
  times.push_back (2.3);
  times.push_back (2.7);
  times.push_back (2.8);
  times.push_back (3.2);
  times.push_back (3.3);
  times.push_back (3.7);
  times.push_back (3.8);
  times.push_back (4.2);
  times.push_back (4.3);
  times.push_back (4.7);
  times.push_back (4.8);
  times.push_back (5.2);
  times.push_back (5.3);
  times.push_back (5.7);
  times.push_back (5.8);
  times.push_back (6.2);
  times.push_back (6.3);
  times.push_back (6.7);
  times.push_back (6.8);
  times.push_back (7.2);
  times.push_back (7.3);
  times.push_back (7.7);
  times.push_back (7.8);
  times.push_back (8.2);
  times.push_back (8.3);
  times.push_back (8.7);
  times.push_back (8.8);
  times.push_back (9.2);
  times.push_back (9.3);
  times.push_back (9.7);
  times.push_back (9.8);
  times.push_back (10.2);
  times.push_back (10.4);

  value_type height = .8;
  SplineBasedPtr_t pg (SplineBased::create (height));
  pg->timeSequence (times);
  pg->stepSequence (steps);

  // define boundary conditions
  vector2_t position; position.setZero ();
  vector2_t velocity; velocity.setZero ();
  pg->setInitialComState (position, velocity);
  position [0] = 3.4;
  pg->setEndComState (position, velocity);

  CubicBSplinePtr_t comTrajectory = pg->solve ();

  std::cout << "set term wxt persist title 'zmp ref' 0 font ',5'" << std::endl
	    << "set xlabel 'x'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 2:3 title 'zmp ref' with points\n";
  std::vector <PiecewisePoly3> zmpRefx = pg->zmpRefx ();
  std::vector <PiecewisePoly3> zmpRefy = pg->zmpRefy ();
  BOOST_CHECK (zmpRefx.size () == zmpRefy.size ());
  for (std::size_t i=0; i<zmpRefx.size (); ++i) {
    const PiecewisePoly3& px = zmpRefx [i];
    const PiecewisePoly3& py = zmpRefy [i];
    value_type t = px.lower;
    value_type delta = (px.upper - px.lower)/6.;
    for (std::size_t j=0; j<7; ++j) {
      std::cout << t << "\t" << px [j] << "\t" << py [j] << std::endl;
      t += delta;
    }
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'zmp ref x' 1 font ',5'" << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'x'" << std::endl;
  std::cout << "plot '-' using 1:2 title 'zmp ref x' with points\n";
  zmpRefx = pg->zmpRefx ();
  zmpRefy = pg->zmpRefy ();
  BOOST_CHECK (zmpRefx.size () == zmpRefy.size ());
  for (std::size_t i=0; i<zmpRefx.size (); ++i) {
    const PiecewisePoly3& px = zmpRefx [i];
    const PiecewisePoly3& py = zmpRefy [i];
    value_type t = px.lower;
    value_type delta = (px.upper - px.lower)/6.;
    for (std::size_t j=0; j<7; ++j) {
      std::cout << t << "\t" << px [j] << "\t" << py [j] << std::endl;
      t += delta;
    }
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'zmp ref y' 2 font ',5'" << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 1:3 title 'zmp ref y' with points\n";
  zmpRefx = pg->zmpRefx ();
  zmpRefy = pg->zmpRefy ();
  BOOST_CHECK (zmpRefx.size () == zmpRefy.size ());
  for (std::size_t i=0; i<zmpRefx.size (); ++i) {
    const PiecewisePoly3& px = zmpRefx [i];
    const PiecewisePoly3& py = zmpRefy [i];
    value_type t = px.lower;
    value_type delta = (px.upper - px.lower)/6.;
    for (std::size_t j=0; j<7; ++j) {
      std::cout << t << "\t" << px [j] << "\t" << py [j] << std::endl;
      t += delta;
    }
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'com trajectory' 3 font ',5'"
	    << std::endl
	    << "set xlabel 'x'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 2:3 title 'com' with lines\n";
  // com
  value_type dt = 0.01;
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    std::cout << t << "\t" << com [0] << "\t" << com [1] << std::endl;
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'com trajectory x' 4 font ',5'"
	    << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'x'" << std::endl;
  std::cout << "plot '-' using 1:2 title 'com x' with lines\n";
  // com x
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    std::cout << t << "\t" << com [0] << "\t" << com [1] << std::endl;
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'com trajectory y' 5 font ',5'"
	    << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 1:3 title 'com y' with lines\n";
  // com y
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    std::cout << t << "\t" << com [0] << "\t" << com [1] << std::endl;
  }
  std::cout << "e" << std::endl;
  // zmp
  std::cout << "set term wxt persist title 'zmp' 6 font ',5'"
	    << std::endl
	    << "set xlabel 'x'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 2:3 title 'zmp' with lines\n";
  value_type un_sur_omega_2 = sqrt (height / SplineBased::gravity);
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    vector_t zmp = com - un_sur_omega_2 * comTrajectory->derivative (t, 2);
    std::cout << t << "\t" << zmp [0] << "\t" << zmp [1] << std::endl;
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'zmp x' 7 font ',5'"
	    << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'x'" << std::endl;
  std::cout << "plot '-' using 1:2 title 'zmp x' with lines\n";
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    vector_t zmp = com - un_sur_omega_2 * comTrajectory->derivative (t, 2);
    std::cout << t << "\t" << zmp [0] << "\t" << zmp [1] << std::endl;
  }
  std::cout << "e" << std::endl;
  std::cout << "set term wxt persist title 'zmp y' 8 font ',5'"
	    << std::endl
	    << "set xlabel 't'" << std::endl
	    << "set ylabel 'y'" << std::endl;
  std::cout << "plot '-' using 1:3 title 'zmp y' with lines\n";
  for (double t = comTrajectory->timeRange ().first;
       t <= comTrajectory->timeRange ().second; t+=dt) {
    vector_t com = (*comTrajectory) (t);
    vector_t zmp = com - un_sur_omega_2 * comTrajectory->derivative (t, 2);
    std::cout << t << "\t" << zmp [0] << "\t" << zmp [1] << std::endl;
  }
  std::cout << "e" << std::endl;
}
}
