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

# sets length of fourier descriptors
# this has to be a multiple of 4 and 3
FDLENGTH = 60
if not (FDLENGTH % 3 == 0 and FDLENGTH % 4 == 0):
	raise RuntimeError("FDLENGTH not OK")


class fdsingle(PluginFunction):
	"""FDS base class"""
	category = "Features"
	self_type = ImageType([ONEBIT])
	return_type = FloatVector(length=FDLENGTH)
	feature_function = True

#fdsingle_complex_position_r1
class fdsingle_complex_position_r1(fdsingle):
	"""Complex position Fourier descriptor for unbroken shapes according to 
Eq. (6) in [Dalitz2013]_. The coefficient c_1 is used for scale normalisation.

The absolute values of the coefficients l(1), l(N-1), l(2), l(N-2), ... are 
returned."""
	pass

#fdsingle_complex position
class fdsingle_complex_position(fdsingle):
	"""Complex position Fourier descriptor for unbroken shapes according to 
Eq. (6) in [Dalitz2013]_. The coefficient c_r is used for scale normalisation.
c_r is the coefficient with the largest absolute value.

The absolute values
of the coefficients l(1), l(N-1), l(2), l(N-2), ... are returned.

.. note:: In the study [Dalitz2013]_, this was the best performing Fourier
          descriptor on the MPEG7 data set of connected shapes.
"""
	pass

#complex_psotion_phase_r1
class fdsingle_complex_position_phase_r1(fdsingle):
	"""Complex position Fourier descriptor for unbroken shapes according to 
Eq. (6) in [Dalitz2013]_. The coefficients c_1 and c_r are used for scale and
phase normalisation. c_r is the coefficient with the largest absolute value.

The complex values of the coefficients l(1), l(N-1) l(2), l(N-2), ... are 
returned as alternating real and imaginary parts."""
	pass

#complex_psotion_phase
class fdsingle_complex_position_phase(fdsingle):
	"""Complex position Fourier descriptor for unbroken shapes according to 
Eq. (6) in [Dalitz2013]_. The coefficients c_r and c_s are used for scale and
phase normalisation. c_r and c_s are the two coefficients with the largest and 
second largest absolute value. (r >= 0, 0 < s < N/2).

The complex values of the coefficients l(1), l(N-1) l(2), l(N-2), ... are
returned as alternating real and imaginary parts."""
	pass

#granlund
class fdsingle_granlund(fdsingle):
	"""Granlund's Fourier descriptor for unbroken shapes according to 
Eq. (7) in [Dalitz2013]_. The coefficient c_1 is used for normalisation.
c_r is the coefficient with the largest absolute value.

The complex values of the coefficients d(2), d(3), ... are returned as
alternating real and imaginary parts."""
	pass

#elliptic
class fdsingle_elliptic(fdsingle):
	"""Elliptic Fourier descriptor for unbroken shapes according to Eq. 
(12) in [Dalitz2013]_. The coefficient I_1 is used for normalisation.

The values of the coefficients I(1), J(1), k1(1), I(2), J(2), K(2), ... are
returned."""
	pass

#real_position
class fdsingle_real_position(fdsingle):
	"""Real position Fourier descriptor for unbroken shapes according to 
Eq. (13) in [Dalitz2013]_. The coefficient c_1 is used for normalisation.

The values of the coefficients with indexes starting at 1 are returned."""
	pass

#centroid_distance
class fdsingle_centroid_distance(fdsingle):
	"""Centroid distance Fourier descriptor for unbroken shapes according
to Eq. (16) in [Dalitz2013]_. The coefficient c_0 is used for scale
normalisation.

The absolute values of the coefficients R(1), R(2), ... are returned."""
	pass


#centroid_distance_phase
class fdsingle_centroid_distance_phase(fdsingle):
	"""Centroid distance Fourier descriptor for unbroken shapes according
to Eq. (16) in [Dalitz2013]_. The coefficients c_0 and c_s are used for scale
and phase normalisation. c_s is the coefficient with the largest absolute
value (s >= 1).

The complex values of the coefficients R(1), R(2), ... are returned 
as alternating real and imaginary parts."""
	pass


class FourierSingle(PluginModule):
    cpp_headers=["dft.hpp", "single.hpp"]
    cpp_namespaces=["FdToolkit"]
    category = "FourierSingle"
    functions = [fdsingle_complex_position,
			fdsingle_complex_position_r1,
			fdsingle_complex_position_phase, fdsingle_granlund, 
			fdsingle_complex_position_phase_r1,
                 fdsingle_elliptic, 
				 fdsingle_centroid_distance,
				 fdsingle_centroid_distance_phase,
				 fdsingle_real_position]

    extra_compile_args = ["-g", "-DFDLENGTH=%d" % FDLENGTH]
    author = "Christian Brandt and Christoph Dalitz"
    url = "http://gamera.sf.net/"

module = FourierSingle()
