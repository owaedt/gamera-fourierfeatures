/*
 * Copyright (C) 2013       Christian Brandt and Christoph Dalitz
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

#ifndef _SINGLE_HPP
#define _SINGLE_HPP

#include "dft.hpp"
#include "plugins/contour.hpp"
#include <numeric>
namespace Gamera {
namespace FdToolkit {
    const double zeroEpsilon = 1e-6;

    /* ******************************************************************** */
    ComplexVector* combinePointsComplexPlaneSequence(PointVector* p) {
        ComplexVector* res = new ComplexVector(p->size());
        for(size_t k = 0; k < p->size(); k++) {
            (*res)[k] = Complex((*p)[k].x(), (*p)[k].y());
        }

        return res;
    }




    /* ******************************************************************** */
    FloatVector* fourierDescriptorComplexPositionR1(ComplexVector* in, int dftCount) {
        ComplexVector* ft = cutComplexDft(in, dftCount+1);
        ft->erase(ft->begin()); //remove c_0
        FloatVector* res = new FloatVector(dftCount);
        double normBy = std::abs(ft->front());
        for(int k = 0; k < dftCount/2; k++) {
          res->at(2*k) = std::abs(ft->at(k)) / normBy;
          res->at(2*k+1) = std::abs(ft->at(dftCount-1-k)) / normBy;
        }
        delete ft;
        return res;
    }


    /* ******************************************************************** */
    // wrapper for contour Points
    FloatVector* fourierDescriptorComplexPositionR1(PointVector* in, int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        FloatVector* res = fourierDescriptorComplexPositionR1(c, dftCount);
        delete c;
        return res;
    } 

    /* ******************************************************************** */
    FloatVector* fourierDescriptorComplexPosition(ComplexVector* in, int dftCount) {
        ComplexVector* ft = cutComplexDft(in, dftCount+1);
        ft->erase(ft->begin()); //remove c_0
        size_t r,s;
        Complex c_r, c_s;
        getCrCsMax(ft, &r, &c_r, &s, &c_s);
        double normBy = std::abs(c_r);
        FloatVector* res = new FloatVector(dftCount);
        for(int k = 0; k < dftCount/2; k++) {
          res->at(2*k) = std::abs(ft->at(k)) / normBy;
          res->at(2*k+1) = std::abs(ft->at(dftCount-1-k)) / normBy;
        }
        delete ft;
        return res;
    }


    /* ******************************************************************** */
    // wrapper for contour Points
    FloatVector* fourierDescriptorComplexPosition(PointVector* in, int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        FloatVector* res = fourierDescriptorComplexPosition(c, dftCount);
        delete c;
        return res;
    } 



    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorComplexPositionPhaseR1(ComplexVector* in, 
            int dftCount) {
        ComplexVector* c_k = cutComplexDft(in, dftCount+1);
        ComplexVector* l_k = new ComplexVector(c_k->size()-1);
        ComplexVector* l_k_sorted = new ComplexVector(c_k->size()-1);

        Complex c_r, c_s;
        size_t r,s;
        getCrCs(c_k, &r, &c_r, &s, &c_s);
        int cutR = kForCutIndex(in->size(), dftCount+1 , r);
        int cutS = kForCutIndex(in->size(), dftCount+1 , s);
        double    argr = std::arg(c_r);
        double    args = std::arg(c_s);
        double absr = std::abs(c_r);
        for(size_t k = 1; k < c_k->size(); k++) {
        //for(size_t k = 1; k < c_k->size(); k++) {
            int cutK = kForCutIndex(in->size(), dftCount+1, k);
            Complex e = std::exp(Complex(0.0, 
                            (cutS - cutR) * std::arg((*c_k)[k]) + 
                            (cutK - cutS) * argr + 
                            (cutR - cutK) * args
                            )
                        );
            (*l_k)[k-1] = (std::abs((*c_k)[k]) / absr) * e;
/*                  std::pow(c_r * std::abs(c_s) / (c_s * std::abs(c_r)), 
                         int((k - r) / (s - r)));*/
                //std::pow((*c_k)[k] / c_r, (k - r) / (s - r));
        }
        for(int k = 0; k < dftCount/2; k++) {
          l_k_sorted->at(2*k) = l_k->at(k);
          l_k_sorted->at(2*k+1) = l_k->at(dftCount-1-k);
        }         

        delete c_k;
        delete l_k;
        return l_k_sorted;
    }


    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorComplexPositionPhaseR1(ComplexVector* in, 
            int dftCount) {
        ComplexVector* fd = complexFourierDescriptorComplexPositionPhaseR1(in, dftCount/2);
        FloatVector* res = complexToRealImagFloat(fd);
        delete fd;
        return res;
    }


    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorComplexPositionPhaseR1(PointVector* in, 
            int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        ComplexVector* res = complexFourierDescriptorComplexPositionPhaseR1(c, dftCount);
        delete c;
        return res;
    }


    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorComplexPositionPhaseR1(PointVector* in, 
            int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        FloatVector* res = floatFourierDescriptorComplexPositionPhaseR1(c, dftCount);
        delete c;
        return res;
    }






    /* ******************************************************************** */

    ComplexVector* complexFourierDescriptorComplexPositionPhase(ComplexVector* in, 
            int dftCount) {
        ComplexVector* c_k = cutComplexDft(in, dftCount+1);
        ComplexVector* l_k = new ComplexVector(c_k->size()-1);
        ComplexVector* l_k_sorted = new ComplexVector(c_k->size()-1);

        Complex c_r, c_s;
        size_t r,s;
        getCrCsMax(c_k, &r, &c_r, &s, &c_s);
        int cutR = kForCutIndex(in->size(), dftCount+1 , r);
        int cutS = kForCutIndex(in->size(), dftCount+1 , s);
        //printf("|c_1|=%.4f, |c_2|=%.4f, |c_3|=%.4f, r=%d, s=%d\n",
        //         std::abs(c_k->at(1)), std::abs(c_k->at(2)), std::abs(c_k->at(3)), cutR, cutS);
        double    argr = std::arg(c_r);
        double    args = std::arg(c_s);
        double absr = std::abs(c_r);
        for(size_t k = 1; k < c_k->size(); k++) {
        //for(size_t k = 1; k < c_k->size(); k++) {
            int cutK = kForCutIndex(in->size(), dftCount+1, k);
            Complex e = std::exp(Complex(0.0, 
                            (cutS - cutR) * std::arg((*c_k)[k]) + 
                            (cutK - cutS) * argr + 
                            (cutR - cutK) * args
                            )
                        );
            (*l_k)[k-1] = (std::abs((*c_k)[k]) / absr) * e;
/*                  std::pow(c_r * std::abs(c_s) / (c_s * std::abs(c_r)), 
                         int((k - r) / (s - r)));*/
                //std::pow((*c_k)[k] / c_r, (k - r) / (s - r));
        }
        for(int k = 0; k < dftCount/2; k++) {
          l_k_sorted->at(2*k) = l_k->at(k);
          l_k_sorted->at(2*k+1) = l_k->at(dftCount-1-k);
        }         

        delete c_k;
        delete l_k;
        return l_k_sorted;
    }


    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorComplexPositionPhase(ComplexVector* in, 
            int dftCount) {
        ComplexVector* fd = complexFourierDescriptorComplexPositionPhase(in, dftCount/2);
        FloatVector* res = complexToRealImagFloat(fd);
        delete fd;
        return res;
    }


    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorComplexPositionPhase(PointVector* in, 
            int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        ComplexVector* res = complexFourierDescriptorComplexPositionPhase(c, dftCount);
        delete c;
        return res;
    }


    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorComplexPositionPhase(PointVector* in, 
            int dftCount) {
        ComplexVector* c = combinePointsComplexPlaneSequence(in);
        FloatVector* res = floatFourierDescriptorComplexPositionPhase(c, dftCount);
        delete c;
        return res;
    }




    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorGranlund(PointVector* p, 
            int dftCount) {

        ComplexVector* in = combinePointsComplexPlaneSequence(p);
        ComplexVector* c_k = cutComplexDft(in, 2*(dftCount+1)+1);
        ComplexVector* d_k = new ComplexVector(dftCount/2);

        Complex norm = c_k->at(1)*c_k->at(1);
        for(size_t k = 2; k <= (size_t)dftCount/2+1; k++) {
          //(*d_k).at(k-2) = ((*c_k).at(k+1) * 
          //                 (*c_k).at(c_k->size() - k + 1) )
          //              / ((*c_k).at(1) * (*c_k).at(1));
          d_k->at(k-2) = c_k->at(k+1)*c_k->at(c_k->size()-k+1) / norm;
        }

        delete c_k;
        delete in;
        FloatVector* res = complexToRealImagFloat(d_k);
        delete d_k;
        return res;
    }

    /* ******************************************************************** */
    template <typename T>
    inline int sgn(T x) {
    if(x > 0)
        return 1;
    else if(x == 0)
        return 0;
    else
        return -1;
    }

    template <typename T>
    inline T square(T x) {
        return x*x;
    }

    void splitPointCoordinates(PointVector* p, 
            FloatVector** x, FloatVector** y) {

        *x = new FloatVector(p->size());
        *y = new FloatVector(p->size());

        for(size_t k = 0; k < p->size(); k++) {
            (**x)[k] = (*p)[k].x();
            (**y)[k] = (*p)[k].y();
        }
    }

    FloatVector* floatFourierDescriptorElliptic(PointVector* p, int numCoeff) {
        FloatVector *a_xk, *a_yk, *b_xk, *b_yk;
        ComplexVector* c_xk, *c_yk;
        FloatVector *x, *y;
        int i = 1;

        splitPointCoordinates(p, &x, &y);

        a_xk = cutRealCosT(x, numCoeff/3+2);
        a_yk = cutRealCosT(y, numCoeff/3+2);
        b_xk = cutRealSinT(x, numCoeff/3+2);
        b_yk = cutRealSinT(y, numCoeff/3+2);
        
        //TODO: find a better way of calculating c_xk and x_yk
        ComplexVector* x_cplx = floatVectorToComplexVector(x);
        ComplexVector* y_cplx = floatVectorToComplexVector(y);
        c_xk = cutComplexDft(x_cplx, numCoeff*2+1);
        c_yk = cutComplexDft(y_cplx, numCoeff*2+1);

        delete x_cplx;
        delete y_cplx;
        delete x;
        delete y;

        FloatVector *res = new FloatVector(numCoeff);

        double axi_ayi, bxi_byi;
        double cyi_abs = std::abs((*c_yk)[i]);
        double cxi_abs = std::abs((*c_xk)[i]);
        axi_ayi = (*a_xk)[i] * (*a_yk)[i];
        bxi_byi = (*b_xk)[i] * (*b_yk)[i];

        for(size_t k = 1; k <= (size_t)numCoeff/3; k++) {
            double Ik=0, Jk=0, Kk=0;
            double axk_ayk, bxk_byk;
            axk_ayk = (*a_xk)[k] * (*a_yk)[k];
            bxk_byk = (*b_xk)[k] * (*b_yk)[k];
            Ik = (*a_xk)[k] * (*a_xk)[k] + (*b_xk)[k] * (*b_xk)[k]
               + (*a_yk)[k] * (*a_yk)[k] + (*b_yk)[k] * (*b_yk)[k];
            Jk = (*a_xk)[k] * (*b_yk)[k] - (*b_xk)[k] * (*a_yk)[k];
            
            {
            k++;
            double cxk_abs = std::abs((*c_xk)[k]);
            double cyk_abs = std::abs((*c_yk)[k]);
            axk_ayk = (*a_xk).at(k) * (*a_yk).at(k);
            bxk_byk = (*b_xk).at(k) * (*b_yk).at(k);
            //axk_ayk = (*a_xk)[k] * (*a_yk)[k];
            //bxk_byk = (*b_xk)[k] * (*b_yk)[k];
            Kk = sgn((axk_ayk + bxk_byk) * (cyi_abs*cyi_abs - cxi_abs*cxi_abs)
                   + (axi_ayi + bxi_byi) * (cxk_abs*cxk_abs - cyk_abs*cyk_abs)
                   ) * 
                (cxi_abs*cxi_abs * cxk_abs * cxk_abs +
                 cyi_abs*cyi_abs * cyk_abs * cyk_abs + 
                 2 * (axi_ayi + bxi_byi) * (axk_ayk + bxk_byk));
            k--;
            }

            (*res)[3*(k-1)] = Ik;
            (*res)[3*(k-1) + 1] = Jk;
            (*res)[3*(k-1) + 2] = Kk;
        }


        double normalizeI1 = (*res)[0];
        if (normalizeI1 < 0.00001) {
          normalizeI1 = 1.0;
          (*res)[0] = 1.0;
        }
        for(size_t k = 0; k < (size_t)numCoeff/3; k++) {
            (*res)[3*k] /= normalizeI1; 
            (*res)[3*k + 1] /= normalizeI1; 
            (*res)[3*k + 2] /= square(normalizeI1); 
        }

        delete a_xk;
        delete a_yk;
        delete b_xk;
        delete b_yk;
        delete c_xk;
        delete c_yk;
        return res;
    }


#define DEBUG_SHRIDHAR 0
    FloatVector* floatFourierDescriptorRealPosition(PointVector* p, 
            int numCoeff) {
        FloatVector *a_xk, *a_yk, *b_xk, *b_yk;
        FloatVector *x, *y;
        FloatVector* res = new FloatVector(numCoeff, 0.0);

        splitPointCoordinates(p, &x, &y);

        a_xk = cutRealCosT(x, 2*numCoeff+1);
        a_yk = cutRealSinT(y, 2*numCoeff+1);
        b_xk = cutRealSinT(x, 2*numCoeff+1);
        b_yk = cutRealCosT(y, 2*numCoeff+1);

        delete x;
        delete y;

        size_t i = 1;
        double norm = square((*a_xk)[i]) + square((*a_yk)[i]) + 
                           square((*b_xk)[i]) + square((*b_yk)[i]);
        for(i = 2; i < (size_t)numCoeff; i++) {
            double curNorm = square((*a_xk)[i]) + square((*a_yk)[i]) + 
                           square((*b_xk)[i]) + square((*b_yk)[i]);
            if(curNorm > norm)
                norm = curNorm;
        }

        norm = sqrt(norm);
        if (norm == 0.0) {
          (*res)[0] = 1.0;
        }
        else {
          for(size_t k = 1; k < (size_t)numCoeff; k++) {
            (*res)[k-1] = sqrt(square((*a_xk)[k]) + square((*b_xk)[k]) + 
                               square((*a_yk)[k]) + square((*b_yk)[k])
                               ) / norm;
          }
        }

        delete a_xk;
        delete a_yk;
        delete b_xk;
        delete b_yk;
        return res;
    }


    ComplexVector* complexFourierDescriptorCentroidDistance(PointVector* p, 
            int numCoeff) {
        size_t pointCount = p->size();
        ComplexVector* r = new ComplexVector(pointCount);

        double meanX=0, meanY=0;
        for(size_t t = 0; t < pointCount; t++) {
            meanX += (*p)[t].x();
            meanY += (*p)[t].y();
        }

        meanX /= pointCount;
        meanY /= pointCount;


        for(size_t t = 0; t < pointCount; t++) {
            (*r)[t] = sqrt(square((*p)[t].x() - meanX) + 
                           square((*p)[t].y() - meanY));
        } 
        ComplexVector* c_k = cutComplexDft(r, 2*numCoeff+1);

        ComplexVector* d_k = new ComplexVector(numCoeff);
        double c0abs = std::abs((*c_k)[0]);
        
        Complex cs = (*c_k)[1];
        double csabs = std::abs(cs);
        size_t s = 1;
        while(csabs < zeroEpsilon && s < c_k->size()) {
            s++;
            cs = (*c_k)[s];
            csabs = std::abs(cs);
        }

        double csarg = std::arg(cs);
        for(size_t k = 0; k < (size_t)numCoeff; k++) {
            (*d_k)[k] = std::exp(Complex(0.0, 
                        s * std::arg((*c_k)[k]) - k * csarg)) 
                * std::abs((*c_k)[k]) / c0abs;
        } 

        delete r;
        delete c_k;
        return d_k;
    }


    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorCentroidDistance(PointVector* in, 
            int dftCount) {
        ComplexVector* r = complexFourierDescriptorCentroidDistance(in, dftCount);
        FloatVector* res = abs(r);
        delete r;
        return res;
    }





    ComplexVector* complexFourierDescriptorCentroidDistancePhase(
            PointVector* p, int numCoeff) {
        size_t pointCount = p->size();
        ComplexVector* r = new ComplexVector(pointCount);

        double meanX=0, meanY=0;
        for(size_t t = 0; t < pointCount; t++) {
            meanX += (*p)[t].x();
            meanY += (*p)[t].y();
        }

        meanX /= pointCount;
        meanY /= pointCount;


        for(size_t t = 0; t < pointCount; t++) {
            (*r)[t] = sqrt(square((*p)[t].x() - meanX) + 
                           square((*p)[t].y() - meanY));
        } 
        ComplexVector* c_k = cutComplexDft(r, 2*numCoeff+1);

        ComplexVector* d_k = new ComplexVector(numCoeff);
        
        Complex cs;
        size_t s = 1;
    
        double csAbsMax = 0.0;
        for(size_t k = 1; k < (size_t)numCoeff; k++) {
            double curAbs = std::abs(c_k->at(k));
            if(curAbs > csAbsMax) {
                csAbsMax = curAbs;
                s = k;
                cs = c_k->at(k);
            }
        }
        //printf("|c_1|=%.4f, |c_2|=%.4f, |c_3|=%.4f, s=%d\n",
        //         std::abs(c_k->at(1)), std::abs(c_k->at(2)), std::abs(c_k->at(3)), s);

        double c0abs = std::abs(c_k->at(0));
        double csarg = std::arg(cs);
        for(size_t k = 0; k < (size_t)numCoeff; k++) {
            (*d_k)[k] = std::exp(Complex(0.0, 
                        s * std::arg((*c_k)[k]) - k * csarg)) 
                * std::abs((*c_k)[k]) / c0abs;
        } 

        delete r;
        delete c_k;
        return d_k;
    }



    /* ******************************************************************** */
    FloatVector* floatFourierDescriptorCentroidDistancePhase(PointVector* in, 
            int dftCount) {
        ComplexVector* r = complexFourierDescriptorCentroidDistancePhase(in, dftCount / 2);
        FloatVector* res = complexToRealImagFloat(r);
        delete r;
        return res;
    }




}

    /* ******************************************************************** */
    template <class T, FloatVector* (*descriptor)(PointVector*, int)>
    inline void imageFourierDescriptor(T &m, feature_t* buf, int dftCount) {
        PointVector* p = contour_pavlidis(m);
        if (p->size() > 1) {
          FloatVector* res = descriptor(p, dftCount);
          for(size_t k = 0; k < res->size(); k++) {
            buf[k] = res->at(k);
          }
          delete res;
        }
        else if (p->size() == 0){
          for(int k = 0; k < dftCount; k++) {
            buf[k] = 0.0;
          }
        }
        else if (p->size() == 1){
          buf[0] = 1.0;
          for(int k = 1; k < dftCount; k++) {
            buf[k] = 0.0;
          }
        }
        delete p;
    }

    template <class T>
    void fdsingle_complex_position_r1(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::fourierDescriptorComplexPositionR1>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdsingle_complex_position(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::fourierDescriptorComplexPosition>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdsingle_complex_position_phase_r1(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorComplexPositionPhaseR1>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdsingle_complex_position_phase(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorComplexPositionPhase>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdsingle_granlund(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorGranlund>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdsingle_elliptic(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorElliptic>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdsingle_real_position(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorRealPosition>(m, buf, FDLENGTH);
    }


    template <class T>
    void fdsingle_centroid_distance(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorCentroidDistance>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdsingle_centroid_distance_phase(T &m, feature_t* buf) {
        imageFourierDescriptor<T, FdToolkit::floatFourierDescriptorCentroidDistancePhase>(m, buf, FDLENGTH);
    }



}


#endif
