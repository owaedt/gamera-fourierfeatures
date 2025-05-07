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

#ifndef _DFT_HPP
#define _DFT_HPP

#include "gamera.hpp"
#include <complex>
#include <cmath>
#include "complexvector.hpp"

using namespace std;


namespace Gamera {
namespace FdToolkit {
    typedef std::complex<double> Complex;
    typedef vector<Complex> ComplexVector;

    /* ******************************************************************** */
    // note that cut offset is the count of value to take from the left + from the right
    ComplexVector* cutComplexDft(ComplexVector* in, int numCoeff) {
        int dftSize = (signed)in->size();
        if(numCoeff % 2 == 0) {
            throw std::runtime_error("even number of coefficients in "
                    "cutComplexDft is not allowed");
        }
        
        ComplexVector *c_k = new ComplexVector(numCoeff);

        int numCoeffHalf = numCoeff/2;
        if(dftSize < numCoeff) {
            numCoeffHalf = dftSize/2;
        }
        int targetIdx = 0;

        for(int k = 0; k <= numCoeffHalf; k++) {
            Complex sum(0.0,0.0);
            Complex prod(1.0,0.0);
            Complex expfac =
                std::exp(Complex(0.0, (-2 * M_PI * k) / dftSize));
            for(int t = 0; t < dftSize; t++) {
                //sum += (*in)[t] * std::exp(
                //        Complex(0.0, -2 * M_PI * k * t / dftSize));
                sum += (*in)[t] * prod;
                prod *= expfac;
            }
            sum /= dftSize;
            //(*c_k).at(targetIdx) = sum;
            (*c_k)[targetIdx] = sum;
            targetIdx++;
        }

        if(dftSize < numCoeff) {
            targetIdx = numCoeff-numCoeffHalf;
        }
        
        for(int k = dftSize - numCoeffHalf; k < dftSize; k++) {
            Complex sum(0.0,0.0);
            Complex prod(1.0,0.0);
            Complex expfac =
                std::exp(Complex(0.0, (-2 * M_PI * k) / dftSize));
            for(int t = 0; t < dftSize; t++) {
                //sum += (*in)[t] * std::exp(
                //        Complex(0.0, -2 * M_PI * k * t / dftSize));
                sum += (*in)[t] * prod;
                prod *= expfac;
            }
            
            sum /= dftSize;
            (*c_k)[targetIdx] = sum; //.at(targetIdx) = sum;
            targetIdx++;
        }

        return c_k;
    }

#define DEBUG_K_FOR_CUT 0

    int kForCutIndex(int dftSize, int numCoeff, int cutIndex) {
#if DEBUG_K_FOR_CUT
      std::cout << "dftSize = " << dftSize << ", numCoeff = "
              << numCoeff << ", cutIndex = " << cutIndex;
#endif

        int k;
        if(dftSize < numCoeff)
            dftSize = numCoeff;

        if(cutIndex <= numCoeff/2)
            k = cutIndex;
        else 
            k = dftSize-numCoeff+cutIndex;

#if DEBUG_K_FOR_CUT
        std::cout << ", k = " << k << std::endl;
#endif

        return k;
    }

#define SINCOS_FILL_MIDDLE 0 // set to 1 for filling transform in middle

    /* ******************************************************************** */
    FloatVector* cutRealCosT(FloatVector* in, int numCoeff) {
        int dftSize = (signed)in->size();
#if SINCOS_FILL_MIDDLE
        int numCoeffToFill = dftSize -numCoeff;
        if(numCoeffToFill <= 0) {
            numCoeff = dftSize;
         }
        FloatVector* res = new FloatVector(numCoeff);
#else
        FloatVector* res = new FloatVector(numCoeff, 0.0);
         if(numCoeff > dftSize/2) {
            numCoeff = dftSize/2;
        }
#endif

        for(int k = 0; k < numCoeff; k++) {
            double sum = 0.0;
            //assumes trailing values as 0 if in->size() < cutOffset
            for(int t = 0; t < dftSize; t++) {
                sum += (*in)[t] * cos(2 * M_PI * k * t / dftSize);
             } 
            sum /= dftSize;
            (*res)[k] = sum;
         }

#if SINCOS_FILL_MIDDLE
        if(numCoeffToFill < 0) {
            res->insert(res->begin() + res->size() / 2+1, 
                    -numCoeffToFill+1, 0.0);
//            std::cout << numCoeffToFill << ":" <<
         }
#else
#endif

        return res;
    } 


    /* ******************************************************************** */
    FloatVector* cutRealSinT(FloatVector* in, int numCoeff) {
        int dftSize = (signed)in->size();
#if SINCOS_FILL_MIDDLE
        int numCoeffToFill = dftSize -numCoeff;
        if(numCoeffToFill <= 0) {
            numCoeff = dftSize;
        }
        FloatVector* res = new FloatVector(numCoeff);
#else
        FloatVector* res = new FloatVector(numCoeff, 0.0);
        if(numCoeff > dftSize/2) {
            numCoeff = dftSize/2;
        }
#endif

        for(int k = 0; k < numCoeff; k++) {
            double sum = 0.0;
            //assumes trailing values as 0 if in->size() < cutOffset
            for(int t = 0; t < dftSize; t++) {
                sum += (*in)[t] * sin(2 * M_PI * k * t / dftSize);
            }
            sum /= dftSize;
            (*res)[k] = sum;
        }

#if SINCOS_FILL_MIDDLE
        if(numCoeffToFill < 0) {
            res->insert(res->begin() + res->size() / 2 + 1, 
                    -numCoeffToFill+1, 0.0);
        }
#else
#endif


        return res;
    } 

}
}
#endif
