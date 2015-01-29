# Copyright (c) 2014, LAAS-CNRS
# Author: Florent Lamiraux
#
# This file is part of hpp-walkgen.
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
# hpp-walkgen. If not, see <http:#www.gnu.org/licenses/>.

#
# This file computes the integral of a polynomial between 0 and 1 as
# a linear combination of the evaluation of the polynomial at regularly
# spaced points between 0 and 1.

from sympy import Matrix, Symbol, Integer
from fractions import Fraction

def computeCoefficients (n):
    m_data = []
    v_data = []

    for i in range (n+1):
        m_data.append ([(Integer(j)/Integer(n))**i for j in range (n+1)])
        v_data.append (1/Integer (i+1))

    M = Matrix (m_data)
    v = Matrix (v_data)
    beta = M**-1*v
    return beta


for n in range (2, 7):
    beta = computeCoefficients (n)

    # Integrale sur [0,1] de x**k = 1/(k+1)

    for k in range (7):
        v_data = [(Integer (j)/Integer (n))**k for j in range (n+1)]
        power = Matrix (v_data)
        x = beta.transpose ()*power
        I = 1/Integer (k+1)
        print ("n = %i, k = %i,\tIntegrale = %s\t linear combination = %s,\t difference = %s"%(n, k, I, x, x[0] - I))

