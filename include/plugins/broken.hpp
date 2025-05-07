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

#ifndef _BROKEN_HPP
#define _BROKEN_HPP

#include "plugins/geometry.hpp"
#include "plugins/segmentation.hpp"
#include "dft.hpp"
#include "plugins/contour.hpp"
#include "geostructs/kdtree.hpp"

namespace Gamera {
namespace FdToolkit {
    typedef std::vector<FloatPoint> FloatPointVector;
    /* ************ */

    /* ******************************************************************** */
    /** write linear interpolation with distance 1 between Points a and b 
     * into *res 
     * */
    void interpolatePoints(FloatPointVector* res, 
            Point a, Point b) {
        FloatPoint q(a.x(), a.y());
        FloatPoint p(b.x(), b.y());

        int dist = (int)q.distance(p);
        FloatPoint dq = (p - q) / FloatPoint(dist,dist);
        
        for(int n = 1; n < dist; n++) {
            q = q + dq;
            res->push_back(q);
        }
        res->push_back(p);
    }


    /* ******************************************************************** */
    /**
     * do linear interpolation with distance 1 of polygon given by points
     * */
    FloatPointVector* interpolatePolygonPoints(PointVector*  points) {
        size_t pointCount = points->size();
        FloatPointVector* res = new FloatPointVector;

        for(size_t i = 0; i < pointCount; i++) { 
            interpolatePoints(res, (*points)[(i-1 + pointCount) % pointCount], 
                    (*points)[i]);
         }

        return res;
    } 

    /* ******************************************************************** */
    /** calculate distance of every point in hullPoints to the nearest point 
     * in contourPoints using a KdTree
     * */
    FloatVector* minimumContourHullDistances(FloatPointVector* hullPoints, 
            PointVector* contourPoints) {
        FloatVector* res = new FloatVector(hullPoints->size());
        Gamera::Kdtree::KdNodeVector nodes;
        for(size_t i = 0; i < contourPoints->size(); i++) {
            Gamera::Kdtree::CoordPoint p;
            p.push_back((*contourPoints)[i].x());
            p.push_back((*contourPoints)[i].y());
            nodes.push_back(KdNode(p));
         }

        Gamera::Kdtree::KdTree tree(&nodes);

        for(size_t i = 0; i < hullPoints->size(); i++) {
            Gamera::Kdtree::KdNodeVector neighbors;
            Gamera::Kdtree::CoordPoint point;
            double x = (*hullPoints)[i].x();
            double y = (*hullPoints)[i].y();
            point.push_back(x); 
            point.push_back(y);
        
            tree.k_nearest_neighbors(point, 1, &neighbors);

            double dx = neighbors[0].point[0] - x;
            double dy = neighbors[0].point[1] - y;
            double dist = sqrt(dx*dx + dy*dy);
            
            if(dist < 1.0) {
                dist = 0.0;
             }

            (*res)[i] = dist;
         }

        return res;    
    }



    /* ******************************************************************** 
     * descriptions of fourier descriptors are given in their python
     * documentation
     * ******************************************************************** */


    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorBrokenAPhaseS1(FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances, int numCoeff) {
        // calculate cplx_dat(t) = r(t) - jd(t)
        ComplexVector* cplx_dat = 
            new ComplexVector(interpolatedHullPoints->size());

        double meanX=0, meanY = 0;
        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            meanX += (*interpolatedHullPoints)[i].x();
            meanY += (*interpolatedHullPoints)[i].y();
        }    

        meanX /= interpolatedHullPoints->size();
        meanY /= interpolatedHullPoints->size();

        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            double x_r = (*interpolatedHullPoints)[i].x() - meanX;
            double y_r = (*interpolatedHullPoints)[i].y() - meanY;
            double r = sqrt(x_r*x_r + y_r*y_r);
            (*cplx_dat)[i] = Complex(r, (*distances)[i]);
        }

        ComplexVector* c_k = cutComplexDft(cplx_dat, numCoeff+1);
        //for (size_t k=0; k<c_k->size(); k++)
        //  printf("%.6f + i%.6f\n", c_k->at(k).real(), c_k->at(k).imag());

        delete cplx_dat;

        //normalization
        double c0abs = std::abs((*c_k)[0]);
        
        size_t s = 1;
        Complex cs = (*c_k)[s];
        double csabs = std::abs(cs);
        const double zeroEpsilon = 1e-6;
        while(csabs < zeroEpsilon && s < c_k->size()-1) {
            s++;
            cs = (*c_k)[s];
            csabs = std::abs(cs);
        }
        s = kForCutIndex(interpolatedHullPoints->size(), numCoeff+1, s);
        //printf("s=%d\n", s);

        ComplexVector* d_k = new ComplexVector(numCoeff+1);
        ComplexVector* d_k_sorted = new ComplexVector(numCoeff);
        double csarg = std::arg(cs);
        for(size_t k = 0; k < (size_t)numCoeff+1; k++) {
          int cutK = kForCutIndex(interpolatedHullPoints->size(), numCoeff+1, k);
          (*d_k)[k] = std::exp(Complex(0.0, 
                s * std::arg((*c_k)[k]) - cutK * csarg)) 
                * std::abs((*c_k)[k]) / c0abs;
        }
        for(int k = 0; k < numCoeff/2; k++) {
          d_k_sorted->at(2*k) = d_k->at(k);
          d_k_sorted->at(2*k+1) = d_k->at(numCoeff-k);
        }        

        delete c_k;
        delete d_k;
        return d_k_sorted;
    }  


    FloatVector* floatFourierDescriptorBrokenAPhaseS1(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenAPhaseS1(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount/2);

        FloatVector* res = complexToRealImagFloat(c_k);
        delete c_k;
        return res;
    }



    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorBrokenAPhase(FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances, int numCoeff) {
        // calculate cplx_dat(t) = r(t) - jd(t)
        ComplexVector* cplx_dat = 
            new ComplexVector(interpolatedHullPoints->size());

        double meanX=0, meanY = 0;
        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            meanX += (*interpolatedHullPoints)[i].x();
            meanY += (*interpolatedHullPoints)[i].y();
        }    

        meanX /= interpolatedHullPoints->size();
        meanY /= interpolatedHullPoints->size();

        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            double x_r = (*interpolatedHullPoints)[i].x() - meanX;
            double y_r = (*interpolatedHullPoints)[i].y() - meanY;
            double r = sqrt(x_r*x_r + y_r*y_r);
            (*cplx_dat)[i] = Complex(r, (*distances)[i]);
        }

        ComplexVector* c_k = cutComplexDft(cplx_dat, numCoeff+1);

        delete cplx_dat;

        Complex cs, cr;
        size_t s,r;
        //getCrCsMax(c_k, &r, &cr, &s, &cs, 0, numCoeff+1);
        // restricting on the first coefficients is better
        getCrCsMax(c_k, &r, &cr, &s, &cs, 0, numCoeff/2);

        //printf("|c_0|=%.4f, |c_1|=%.4f, |c_2|=%.4f, ", std::abs(c_k->at(0)), std::abs(c_k->at(1)), std::abs(c_k->at(2)));
        //printf("r=%d, ", r);
        s = kForCutIndex(interpolatedHullPoints->size(), numCoeff+1, s);
        //printf("s=%d, arg(c_1)=%.4f, arg(c_s)=%.4f\n", s, std::arg(c_k->at(1)), std::arg(cs));
        
        ComplexVector* d_k = new ComplexVector(numCoeff+1);
        double csarg = std::arg(cs);
        double crabs = std::abs(cr);
        ComplexVector* d_k_sorted = new ComplexVector(numCoeff);
        for(size_t k = 0; k < (size_t)numCoeff+1; k++) {
          int cutK = kForCutIndex(interpolatedHullPoints->size(), numCoeff+1, k);
          (*d_k)[k] = std::exp(Complex(0.0, 
                s * std::arg((*c_k)[k]) - cutK * csarg)) 
                * std::abs((*c_k)[k]) / crabs;
        }
        for(int k = 0; k < numCoeff/2; k++) {
          d_k_sorted->at(2*k) = d_k->at(k);
          d_k_sorted->at(2*k+1) = d_k->at(numCoeff-k);
        }        

        delete c_k;
        delete d_k;
        return d_k_sorted;
    }


    FloatVector* floatFourierDescriptorBrokenAPhase(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenAPhase(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount/2);

        FloatVector* res = complexToRealImagFloat(c_k);
        delete c_k;
        return res;
    }

    FloatVector* floatFourierDescriptorBrokenA(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenAPhase(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount);

        FloatVector* res = abs(c_k);
        delete c_k;
        return res;
    }


    FloatVector* floatFourierDescriptorBrokenB(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {
        // calculate cplx_dat(t) = r(t) - jd(t)
        ComplexVector* x = 
            new ComplexVector(interpolatedHullPoints->size());
        ComplexVector* y = 
            new ComplexVector(interpolatedHullPoints->size());
        ComplexVector* d = 
            new ComplexVector(interpolatedHullPoints->size());

        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            (*x)[i] = (*interpolatedHullPoints)[i].x();
            (*y)[i] = (*interpolatedHullPoints)[i].y();
            (*d)[i] = (*distances)[i];
        }

        ComplexVector* c_x = cutComplexDft(x, 2*(dftCount+1)+1);
        ComplexVector* c_y = cutComplexDft(y, 2*(dftCount+1)+1);
        ComplexVector* c_d = cutComplexDft(d, 2*(dftCount+1)+1);
        
        delete x;
        delete y;
        delete d;

        double cx1abs = std::abs((*c_x)[1]);
        double cy1abs = std::abs((*c_y)[1]);
        double norm = sqrt(cx1abs*cx1abs + cy1abs * cy1abs);

        FloatVector* res = new FloatVector(dftCount);
        for(int k = 1; k <= dftCount; k++) {
            double cxk, cyk, cdk;
            cxk = std::abs((*c_x)[k]);
            cyk = std::abs((*c_y)[k]);
            cdk = std::abs((*c_d)[k]);
            (*res)[k-1] = sqrt(cxk*cxk + cyk*cyk + cdk*cdk) / norm;
        }

        delete c_x;
        delete c_y;
        delete c_d;

        return res;
    }


    /* ******************************************************************** */
    ComplexVector* complexFourierDescriptorBrokenCPhaseS1(FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances, int numCoeff) {
        // calculate cplx_dat(t) = r(t) - jd(t)
        ComplexVector* cplx_dat = 
            new ComplexVector(interpolatedHullPoints->size());

        double meanX=0, meanY = 0;
        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            meanX += (*interpolatedHullPoints)[i].x();
            meanY += (*interpolatedHullPoints)[i].y();
        }  

        meanX /= interpolatedHullPoints->size();
        meanY /= interpolatedHullPoints->size();

        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            double x_r = (*interpolatedHullPoints)[i].x() - meanX;
            double y_r = (*interpolatedHullPoints)[i].y() - meanY;
            double r = sqrt(x_r*x_r + y_r*y_r);
            (*cplx_dat)[i] = Complex(r -  (*distances)[i], 0);
        }

        ComplexVector* c_k = cutComplexDft(cplx_dat, (2*numCoeff)+1);



        //normalization
        double c0abs = 0;
        for(size_t k = 0; k < cplx_dat->size(); k++) {
          c0abs += std::max<double>(0, cplx_dat->at(k).real());
        }
        c0abs /= cplx_dat->size();
        
        delete cplx_dat;
        
        size_t s = 1;
        Complex cs = (*c_k)[s];
        double csabs = std::abs(cs);
        const double zeroEpsilon = 1e-6;
        while(csabs < zeroEpsilon && s < c_k->size()-1) {
            s++;
            cs = (*c_k)[s];
            csabs = std::abs(cs);
        }

        ComplexVector* d_k = new ComplexVector(numCoeff);
        double csarg = std::arg(cs);
        for(size_t k = 0; k < (size_t)numCoeff; k++) {
            int cutK = k; //kForCutIndex(interpolatedHullPoints->size(), numCoeff+1 , k);
            (*d_k)[k] = std::exp(Complex(0.0, 
                        s * std::arg((*c_k)[k]) - cutK * csarg)) 
                * std::abs((*c_k)[k]) / c0abs;
        }  
        delete c_k;
        return d_k;
    }


    FloatVector* floatFourierDescriptorBrokenCPhaseS1(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenCPhaseS1(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount/2);

        FloatVector* res = complexToRealImagFloat(c_k);
        delete c_k;
        return res;
    }





    ComplexVector* complexFourierDescriptorBrokenCPhase(FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances, int numCoeff) {
        // calculate cplx_dat(t) = r(t) - jd(t)
        ComplexVector* cplx_dat = 
            new ComplexVector(interpolatedHullPoints->size());

        double meanX=0, meanY = 0;
        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            meanX += (*interpolatedHullPoints)[i].x();
            meanY += (*interpolatedHullPoints)[i].y();
        }  

        meanX /= interpolatedHullPoints->size();
        meanY /= interpolatedHullPoints->size();

        for(size_t i = 0; i < interpolatedHullPoints->size(); i++) {
            double x_r = (*interpolatedHullPoints)[i].x() - meanX;
            double y_r = (*interpolatedHullPoints)[i].y() - meanY;
            double r = sqrt(x_r*x_r + y_r*y_r);
            (*cplx_dat)[i] = Complex(r -  (*distances)[i], 0);
        }

        ComplexVector* c_k = cutComplexDft(cplx_dat, (2*numCoeff)+1);



        //normalization
        /**************************
        double c0abs = 0;
        for(size_t k = 0; k < numCoeff; k++) {
            c0abs += std::max<double>(0, c_k->at(k).real());
        }
        c0abs /= numCoeff;
        
        delete cplx_dat;
        
        size_t s = 1;
        Complex cs = (*c_k)[s];
        double csabs = std::abs(cs);
        const double zeroEpsilon = 1e-6;
        for(size_t k  = 2; k < numCoeff; k++) {
            double a = std::abs(c_k->at(k));
            if(a > csabs) {
                cs = c_k->at(k);
                csabs = a;
                s = k;
            }
        }
        **********************************/
        Complex cs,cr;
        size_t r,s;
        getCrCsMax(c_k, &r, &cr, &s, &cs, 0, numCoeff);
        double c0abs = std::abs(cr);
        //printf("|c_0|=%.4f, |c_1|=%.4f, |c_2|=%.4f, ", std::abs(c_k->at(0)), std::abs(c_k->at(1)), std::abs(c_k->at(2)));
        //printf("r=%d, ", r);
        //printf("s=%d, arg(c_1)=%.4f, arg(c_s)=%.4f\n", s, std::arg(c_k->at(1)), std::arg(cs));

        ComplexVector* d_k = new ComplexVector(numCoeff);
        double csarg = std::arg(cs);
        for(size_t k = 0; k < (size_t)numCoeff; k++) {
            int cutK = k;
            (*d_k)[k] = std::exp(Complex(0.0, 
                        s * std::arg((*c_k)[k]) - cutK * csarg)) 
                * std::abs((*c_k)[k]) / c0abs;
        }  
        delete c_k;
        return d_k;
    }


    FloatVector* floatFourierDescriptorBrokenCPhase(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenCPhase(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount/2);

        FloatVector* res = complexToRealImagFloat(c_k);
        delete c_k;
        return res;
    }

    FloatVector* floatFourierDescriptorBrokenC(
            FloatPointVector* interpolatedHullPoints, 
            PointVector* contourPoints, FloatVector* distances,
            int dftCount) {

        ComplexVector* c_k = complexFourierDescriptorBrokenCPhase(
                interpolatedHullPoints, 
                contourPoints, distances, dftCount);

        FloatVector* res = abs(c_k);
        delete c_k;
        return res;
    }





}


    /* ******************************************************************** */
    template <class T, FloatVector* (*descriptor)(FdToolkit::FloatPointVector*, PointVector*, FloatVector*, int) >
    inline void imageFourierDescriptorContourHullDistance(T &m, feature_t* buf, 
            int dftCount) {

		typename ImageFactory<T>::view_type *tmp = simple_image_copy(m);
        ImageList* ccs = cc_analysis(*tmp);
        PointVector* p = new PointVector;

        for(ImageList::iterator cc_it = ccs->begin(); cc_it != ccs->end(); 
                 cc_it++) {
            Cc* cc = static_cast<Cc*>(*cc_it);
            Point orig = cc->origin();
            PointVector* cc_p = contour_pavlidis(*cc);
            
            for(PointVector::iterator p_it = cc_p->begin(); 
                     p_it != cc_p->end(); p_it++) {
                
                p->push_back(*p_it + orig);
            }
            delete *cc_it;
            delete cc_p;
        }
        delete ccs;
		delete tmp->data();
		delete tmp;

		if(p->size() == 0) {
          for (int k = 0; k < dftCount; k++) 
            buf[k] = 0;
          delete p;
          return;
		}
		else if(p->size() == 1) {
          buf[0] = 1.0;
          for (int k = 1; k < dftCount; k++) 
            buf[k] = 0;
          delete p;
          return;
		}


        PointVector* hullPoints = convex_hull_from_points(p);
        FdToolkit::FloatPointVector* interpolatedHullPoints = 
            FdToolkit::interpolatePolygonPoints(hullPoints);
        delete hullPoints;

#if DEBUG_BROKEN_DISTANCE        
        std::cout << "|CH_i| = " << interpolatedHullPoints->size() << std::endl;
        std::cout << "|Cont| = " << p->size() << std::endl;
#endif

        FloatVector* distances = FdToolkit::minimumContourHullDistances(
                interpolatedHullPoints, p);

#if DEBUG_BROKEN_DISTANCE
        std::cout << "dist = ["; 
        for(int i = 0; i < distances->size(); i++)
            std::cout 
                << interpolatedHullPoints->at(i).x() << ", " 
                << interpolatedHullPoints->at(i).y() << ", " 
                << distances->at(i) <<  ";";
        std::cout << "]\n";

//        std::cout << p->size() << std::endl;
#endif

        FloatVector* res = descriptor(interpolatedHullPoints, p, distances,
                dftCount);

        for(size_t k = 0; k < res->size(); k++) {
            buf[k] = (*res)[k];
        }
        delete interpolatedHullPoints;
        delete distances;
        delete p;
        delete res;
    }

    template <class T>
    void fdbroken_a_phase_s1(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenAPhaseS1>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdbroken_a_phase(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenAPhase>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdbroken_a(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenA>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdbroken_b(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenB>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdbroken_c_phase_s1(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenCPhaseS1>(m, buf, FDLENGTH);
    }
    
    template <class T>
    void fdbroken_c_phase(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenCPhase>(m, buf, FDLENGTH);
    }

    template <class T>
    void fdbroken_c(T &m, feature_t* buf) {
        imageFourierDescriptorContourHullDistance<T, FdToolkit::floatFourierDescriptorBrokenC>(m, buf, FDLENGTH);
    }
    



}


#endif
