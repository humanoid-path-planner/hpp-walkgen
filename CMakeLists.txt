#
# Copyright (c) 2015, 2020 CNRS
# Authors: Florent Lamiraux, Guilhem Saurel
#
#
# This file is part of hpp-walkgen
# hpp-walkgen is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# hpp-walkgen is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# hpp-walkgen  If not, see
# <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

SET(PROJECT_NAME hpp-walkgen)
SET(PROJECT_DESCRIPTION "Walk generator for hpp framework")

SET(PROJECT_USE_CMAKE_EXPORT TRUE)
SET(PROJECT_USE_KEYWORD_LINK_LIBRARIES TRUE)
SET(CXX_DISABLE_WERROR TRUE)

INCLUDE(cmake/hpp.cmake)

COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

ADD_PROJECT_DEPENDENCY("hpp-core" REQUIRED)
ADD_REQUIRED_DEPENDENCY("roboptim-trajectory >= 1")

IF(${ROBOPTIM_TRAJECTORY_VERSION} VERSION_GREATER "3.1")
  ADD_DEFINITIONS(-DROBOPTIM_TRAJECTORY_32)
  PKG_CONFIG_APPEND_CFLAGS("-DROBOPTIM_TRAJECTORY_32")
ENDIF(${ROBOPTIM_TRAJECTORY_VERSION} VERSION_GREATER "3.1")

SET(${PROJECT_NAME}_HEADERS
  include/hpp/walkgen/bspline-based.hh
  include/hpp/walkgen/foot-print.hh
  include/hpp/walkgen/foot-trajectory.hh
  include/hpp/walkgen/fwd.hh
)

SET(${PROJECT_NAME}_SOURCES
  src/bspline-based.cc
  src/foot-trajectory.cc
  )

ADD_LIBRARY(${PROJECT_NAME} SHARED ${${PROJECT_NAME}_SOURCES} ${${PROJECT_NAME}_HEADERS})
TARGET_LINK_LIBRARIES (${PROJECT_NAME} PUBLIC hpp-core::hpp-core)
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} roboptim-trajectory)

INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)

ADD_SUBDIRECTORY(tests)

PKG_CONFIG_APPEND_LIBS(${PROJECT_NAME})

IF(INSTALL_DOCUMENTATION)
  INSTALL(FILES doc/walkgen.pdf DESTINATION
    share/doc/${PROJECT_NAME}/doxygen-html/figures)
ENDIF(INSTALL_DOCUMENTATION)

INSTALL(FILES package.xml DESTINATION share/${PROJECT_NAME})
