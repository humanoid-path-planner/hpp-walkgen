#
# Copyright (c) 2015 CNRS
# Authors: Florent Lamiraux
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CXX_DISABLE_WERROR TRUE)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)

SET(PROJECT_NAME hpp-walkgen)
SET(PROJECT_DESCRIPTION
  "Walk generator for hpp framework"
  )
SET(PROJECT_URL "")

SETUP_PROJECT()

# Activate hpp-util logging if requested
SET (HPP_DEBUG FALSE CACHE BOOL "trigger hpp-util debug output")
IF (HPP_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHPP_DEBUG")
ENDIF()

SET(${PROJECT_NAME}_HEADERS
  include/hpp/walkgen/bspline-based.hh
  include/hpp/walkgen/foot-print.hh
  include/hpp/walkgen/foot-trajectory.hh
  include/hpp/walkgen/fwd.hh
)

ADD_REQUIRED_DEPENDENCY(hpp-util >= 3)
ADD_REQUIRED_DEPENDENCY(hpp-core >= 3)
ADD_REQUIRED_DEPENDENCY(roboptim-trajectory >= 1)

PKG_CONFIG_APPEND_LIBS("hpp-walkgen")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

INSTALL(FILES doc/walkgen.pdf
  DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/doc/${PROJECT_NAME}/doxygen-html/figures
  )

SETUP_PROJECT_FINALIZE()
