/*
 * Copyright (C) 2013      Christian Brandt and Christoph Dalitz
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef _COMPLEX_VECTOR_HPP
#define _COMPLEX_VECTOR_HPP

#include "gamera.hpp"
#include <complex>
#include <cmath>

using namespace std;


namespace Gamera {
namespace FdToolkit {

	typedef std::complex<double> Complex;
	typedef std::vector<Complex> ComplexVector;

	/* ******************************************************************** */
	FloatVector* complexToRealImagFloat(ComplexVector* v) {
		size_t inSize = v->size();
		FloatVector* result = new FloatVector(inSize*2);

		for(size_t k = 0; k < inSize; k++) {
			(*result)[2*k] = (*v)[k].real();
			(*result)[2*k+1] = (*v)[k].imag();
		}

		return result;
	}



	/* ******************************************************************** */
	FloatVector* abs(ComplexVector* in) {
		FloatVector* res = new FloatVector(in->size());

		for(size_t k = 0; k < in->size(); k++) {
			(*res)[k] = std::abs((*in)[k]);
		}

		return res;
	}

    ComplexVector* floatVectorToComplexVector(FloatVector*v) {
		size_t inSize = v->size();
		ComplexVector* result = new ComplexVector(inSize);

		for(size_t k = 0; k < inSize; k++) {
			(*result)[k] = Complex((*v)[k], 0);
		}

		return result;
    }

	void getCrCsMax(ComplexVector* c_k, size_t* r, Complex* c_r, 
                    size_t* s, Complex* c_s, size_t start=1, size_t end=0) {
		size_t max=0, secondMax = 0;
		double absMax = 0,
			   absSecondMax = 0;
        if (end == 0) end = c_k->size();

	    for(size_t i = start; i < end; i++) {

			double curAbs = std::abs((*c_k)[i]);
			if(curAbs > absMax) { 
				absSecondMax = absMax;
				secondMax = max;
				max = i;
				absMax = curAbs;
			} else if (curAbs > absSecondMax) {
				absSecondMax = curAbs;
				secondMax = i;
			}

	    }

		*r = max;
		*s = secondMax;
		*c_r = c_k->at(max);
		*c_s = c_k->at(secondMax);
    }

    //get two non-zero factors
	void getCrCs(ComplexVector* c_k, size_t* r, Complex* c_r, 
		size_t* s, Complex* c_s) {
	    size_t N = c_k->size();
	    *r = 1;
	    *c_r = (*c_k)[*r];
	    for(size_t i = 1; i < c_k->size()/2; i++) {
		    *s = N-i;
		    *c_s = (*c_k)[*s];
	     	if(std::abs(*c_s) > 1e-5) {
			    break;
		    }

		    *s = i+1;
		    *c_s = (*c_k)[*s];
	     	if(std::abs(*c_s) > 1e-5) {
			    break;
		    }
	    }
    }


}
}
#endif
