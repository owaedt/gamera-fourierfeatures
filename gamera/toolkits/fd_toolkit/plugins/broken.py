#
# Copyright (C) 2013      Christian Brandt, Christoph Dalitz
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from gamera.plugin import * 

#sets length of fourier descirptors
# this has to be a multiple of 4 and 3
FDLENGTH = 60
if not (FDLENGTH % 3 == 0 and FDLENGTH % 4 == 0):
	raise RuntimeError("FDLENGTH not OK")

class fdbroken(PluginFunction):
	category = "Features"
	self_type = ImageType([ONEBIT])
	return_type = FloatVector(length=FDLENGTH)
	feature_function = True


#fdbroken_a_phase_s1
class fdbroken_a_phase_s1(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (19) in
[Dalitz2013]_. 
The coefficients c_0 and c_1 are used for scale and phase normalisation. 

The complex values of the coefficients A(0), A(N/2-1) A(1), A(N/2-2), ... 
are returned as alternating real and imaginary parts."""
	pass

#fdbroken_a_phase
class fdbroken_a_phase(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (19) in
[Dalitz2013]_. The coefficients c_r and c_s are used for scale and phase
normalisation. c_r and c_s are the two coefficients with the largest and
second largest value (r >= 0, 0 < s < N/2).

The complex values of the coefficients A(0), A(N/2-1) A(1), A(N/2-2), ...
are returned as alternating real and imaginary parts."""
	pass	


#fdbroken_a
class fdbroken_a(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq.(19) in
[Dalitz2013]_. The coefficient c_0 is used for scale normalisation.

The absolute values of 
the coefficients A(0), A(N-1) A(1), A(N-2), ... are returned.

Reference:

  .. [Dalitz2013] Dalitz, C., Brandt, C., Goebbels, S., Kolanus, D.:
     "Fourier Descriptors for Broken Shapes." EURASIP Journal on
     Advances in Signal Processing 2013:161, 2013

.. note:: In this study, *fdbroken_a* was the best performing Fourier
          descriptor on the NEUMES data set of broken shapes.
"""
	pass
	
#fdbroken_b
class fdbroken_b(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (20) in
[Dalitz2013]_.

The values of the coefficients B(1), B(2), ... B(N) are returned."""
	pass

#fdbroken_c_phase_s1
class fdbroken_c_phase_s1(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (22) in
[Dalitz2013]_. The coefficients c_0 = 1/N  sum(r(t)) is used for scale
normalisation. The coefficient c_1 is used for phase normalization.

The complex values of the coefficients C(0), ... C(N/4-1) are returned as
alternating real and imaginary parts."""
	pass

#fdbroken_c_phase
class fdbroken_c_phase(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (22) in
[Dalitz2013]_. The coefficients c_r and c_s are used for scale and phase
normalisation. c_r and c_s are the two coefficients with the largest and
second largest absolute value (r >= 0, 0 < s < N/2).

The complex values of the coefficients C(0), 
... C(N/4-1) are returned as alternating real and imaginary parts."""
	pass

#fdbroken_c
class fdbroken_c(fdbroken):
	"""Fourier descriptor for broken shapes according to Eq. (22) in
[Dalitz2013]_. The coefficient c_0 is used for scale normalisation.

The absolute values of the coefficients C(0), ... C(N/2-1) are returned."""
	pass





class FourierBroken(PluginModule):
    import os, sys, gamera
    gamera_root = os.path.dirname(gamera.__file__)
    cppfile = os.path.join(gamera_root,"src/geostructs/kdtree.cpp")
    if not os.path.exists(cppfile):
        gamera_root = os.path.abspath(os.path.join(sys.exec_prefix,"gamera"))
        cppfile = os.path.join(gamera_root,"src/geostructs/kdtree.cpp")
    cpp_sources=[cppfile, 
			os.path.join(gamera_root,"src/geostructs/delaunaytree.cpp")
			]

    cpp_headers=["dft.hpp", "broken.hpp"]
    cpp_namespaces=["FdToolkit"]
    category = "FourierBroken"
    functions = [fdbroken_a, fdbroken_a_phase, 
			fdbroken_a_phase_s1,	fdbroken_b, 
			fdbroken_c, fdbroken_c_phase, 
			fdbroken_c_phase_s1,
			]
    extra_compile_args = ["-g", "-DFDLENGTH=%d" % FDLENGTH]
    author = "Christian Brandt and Christoph Dalitz"
    url = "http://gamera.sf.net/"

module = FourierBroken()
