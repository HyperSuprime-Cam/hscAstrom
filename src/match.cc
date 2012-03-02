#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include "hsc/meas/astrom/match.h"
#include "lsst/utils/ieee.h"

#include "boost/scoped_array.hpp"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::astrom;
using namespace lsst::afw::table;

// Algorithm is based on V.Tabur 2007, PASA, 24, 189-198
// "Fast Algorithms for Matching CCD Images to a Stellar Catalogue"

// Compair Source based on its PsfFlux
// Ordering is bright to faint
bool cmpSrc(SourceRecord const& a, SourceRecord const& b) {
    float aFlux = a.getPsfFlux();
    float bFlux = b.getPsfFlux();
    if (lsst::utils::isnan(aFlux)) {
        aFlux = 0.0;
    }
    if (lsst::utils::isnan(bFlux)) { 
        bFlux = 0.0;
    }
    return aFlux > bFlux;
}

bool cmpPair(SourcePair const &a, SourcePair const &b) {
    return a.distance > b.distance;
}


SourceCatalog hsc::meas::astrom::selectPoint(SourceCatalog const &a,
                                             unsigned int num,
                                             unsigned int start) {
    // copy and sort array of pointers on psfFlux
    SourceCatalog b(a);
    b.sort(b.getTable()->getPsfFluxKey());
    if (start == 0) {
	if (num < b.size()) {
	    // chop off the end
            b = SourceCatalog(b.getTable(), b.begin(), b.begin() + num, true);
	}
        return b;
    }
    if (start + num < b.size()) {
        return SourceCatalog(b.getTable(), b.begin() + start, b.begin() + start +  num);
    }
    return SourceCatalog(b.getTable(), b.begin() + start, b.end());
}

std::vector<SourcePair> searchPair(std::vector<SourcePair> &a, SourcePair &p, double e) {
    std::vector<SourcePair> v;

    for (size_t i = 0; i < a.size(); i++) {
	double dd = fabs(a[i].distance - p.distance);
	double dpa = fabs(a[i].pa - p.pa);
	if (dd < e && dpa < 0.03) {
//	if (dd < e) {
	    a[i].deltaD = dd;
	    v.push_back(a[i]);
	}
    }

    return v;
}

std::vector<SourcePair>::iterator searchPair3(std::vector<SourcePair> &a,
					      const SourcePair &p,
					      const SourcePair &q,
					      const double &e,
					      const double &dpa,
					      const double &e_dpa = 0.03) {
    std::vector<SourcePair>::iterator idx = a.end();
    double dd_min = 1.E+10;
    //double dpa_min = e_dpa;

    for (std::vector<SourcePair>::iterator i = a.begin(); i < a.end(); i++) {
	double dd = fabs(i->distance - p.distance);
#if 1
	if (dd < e &&
	    fabs(p.pa - i->pa - dpa) < e_dpa &&
	    dd < dd_min &&
	    (i->first == q.first)) {
	    dd_min = dd;
	    idx = i;
	}
#else
	if (dd < e &&
	    fabs(p.pa - i->pa - dpa) < dpa_min) {
	    dpa_min = fabs(p.pa - i->pa - dpa);
	    idx = i;
	}
#endif
    }

    return idx;
}

void transform(int order, double *coeff, double x, double y,
	       double *xn, double *yn) {
    int ncoeff = (order + 1) * (order + 2) / 2;
    *xn = 0.0;
    *yn = 0.0;
    int n = 0;
    for (int i = 0; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    *xn += coeff[n] * pow(x, j) * pow(y, k);
	    *yn += coeff[n+ncoeff] * pow(x, j) * pow(y, k);
	    n++;
	}
    }
}

double *polyfit(int order,
		SourceCatalog const &img,
		SourceCatalog const &cat) {
    int ncoeff = (order + 1) * (order + 2) / 2;
    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 0; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    int *flag = new int[img.size()];
    for (size_t k = 0; k < img.size(); k++) {
	flag[k] = 1;
    }

    double *a_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    double *coeff = new double[ncoeff*2];

    for (int loop = 0; loop < 1; loop++) {
	for (int i = 0; i < ncoeff; i++) {
	    for (int j = 0; j < ncoeff; j++) {
		a_data[i*ncoeff+j] = 0.0;
		for (size_t k = 0; k < img.size(); k++) {
		    if (flag[k] == 1) {
			a_data[i*ncoeff+j] += pow(img[k].getX(), xorder[i]) * 
			    pow(img[k].getY(), yorder[i]) * 
			    pow(img[k].getX(), xorder[j]) * 
			    pow(img[k].getY(), yorder[j]);
		    }
		}
	    }
	    b_data[i] = c_data[i] = 0.0;
	    for (unsigned int k = 0; k < img.size(); k++) {
		if (flag[k] == 1) {
		    b_data[i] += pow(img[k].getX(), xorder[i]) * 
			pow(img[k].getY(), yorder[i]) * 
			cat[k].getX();
		    c_data[i] += pow(img[k].getX(), xorder[i]) * 
			pow(img[k].getY(), yorder[i]) * 
			cat[k].getY();
		}
	    }
	}

	gsl_matrix_view a = gsl_matrix_view_array(a_data, ncoeff, ncoeff);
	gsl_vector_view b = gsl_vector_view_array(b_data, ncoeff);
	gsl_vector_view c = gsl_vector_view_array(c_data, ncoeff);

	gsl_vector *x = gsl_vector_alloc(ncoeff);
	gsl_vector *y = gsl_vector_alloc(ncoeff);

	int s;

	gsl_permutation *p = gsl_permutation_alloc(ncoeff);

	gsl_linalg_LU_decomp(&a.matrix, p, &s);
	gsl_linalg_LU_solve(&a.matrix, p, &b.vector, x);
	gsl_linalg_LU_solve(&a.matrix, p, &c.vector, y);

	for (int i = 0; i < ncoeff; i++) {
	    coeff[i] = x->data[i];
	    coeff[i+ncoeff] = y->data[i];
	}

	double S, Sx, Sy, Sxx, Syy;
	S = Sx = Sy = Sxx = Syy = 0.0;
	for (size_t k = 0; k < img.size(); k++) {
	    if (flag[k] == 1) {
		double x0 = img[k].getX();
		double y0 = img[k].getY();
		double x1, y1;
		transform(order, coeff, x0, y0, &x1, &y1);
		S   += 1.;
		Sx  += (x1 - cat[k].getX());
		Sxx += (x1 - cat[k].getX()) * (x1 - cat[k].getX());
		Sy  += (y1 - cat[k].getY());
		Syy += (y1 - cat[k].getY()) * (y1 - cat[k].getY());
	    }
	}
	double x_sig = sqrt((Sxx - Sx * Sx / S) / S);
	double y_sig = sqrt((Syy - Sy * Sy / S) / S);
	//std::cout << x_sig << " " << y_sig << std::endl;
    
	for (size_t k = 0; k < img.size(); k++) {
	    double x0 = img[k].getX();
	    double y0 = img[k].getY();
	    double x1, y1;
	    transform(order, coeff, x0, y0, &x1, &y1);
	    if (fabs(x1-cat[k].getX()) > 2. * x_sig ||
		fabs(y1-cat[k].getY()) > 2. * y_sig) {
		flag[k] = 0;
	    }
	}

	gsl_permutation_free(p);
	gsl_vector_free(x);
	gsl_vector_free(y);
    }

    delete [] flag;
    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] b_data;
    delete [] c_data;

    return coeff;
}

SourceCatalog::const_iterator searchNearestPoint(SourceCatalog const &cat, double& x, double& y, double& e) {
    SourceCatalog::const_iterator index = cat.end();
    double d_min, d;

    d_min = e;

    for (SourceCatalog::const_iterator i = cat.begin(); i < cat.end(); i++) {
	d = sqrt(pow(((*i).getX()-x), 2) + pow(((*i).getY()-y), 2));
	if (d < d_min) {
	    index = i;
	    d_min = d;
	    break;
	}
    }

    return index;
}

SourceMatchVector FinalVerify(double *coeff,
                              SourceCatalog const &src,
                              SourceCatalog const &cat,
                              double matchingAllowanceInPixel,
                              bool verbose) {
    SourceCatalog srcMat;
    SourceCatalog catMat;
    SourceMatchVector matPair;

    double x0, y0, x1, y1;
    int num = 0, num_prev = -1;
    for (size_t i = 0; i < src.size(); i++) {
	x0 = src[i].getX();
	y0 = src[i].getY();
	transform(1, coeff, x0, y0, &x1, &y1);
	double e = matchingAllowanceInPixel;
	SourceCatalog::const_iterator p = searchNearestPoint(cat, x1, y1, e);
	if (p != cat.end()) {
	    num++;
	    srcMat.push_back(src[i]);
	    catMat.push_back(*p);
	}
    }

    //std::cout << num << std::endl;
    int order = 1;
    if (num > 5) {
	coeff = polyfit(order, srcMat, catMat);

	for (int j = 0; j < 100; j++) {
	    srcMat.clear();
	    catMat.clear();
	    matPair.clear();
	    num = 0;
	    for (size_t i = 0; i < src.size(); i++) {
		x0 = src[i].getX();
		y0 = src[i].getY();
		transform(order, coeff, x0, y0, &x1, &y1);
		double e = matchingAllowanceInPixel;
		SourceCatalog::const_iterator p = searchNearestPoint(cat, x1, y1, e);
		if (p != cat.end()) {
		    num++;
		    srcMat.push_back(src[i]);
		    catMat.push_back(*p);
		    matPair.push_back(SourceMatch(p, src.get(i), 0.0));
		}
	    }
	    //std::cout << "# of objects matched: " << num << " " << num_prev << std::endl;
	    if (num == num_prev) break;
	    //if (num > 50) order = 3;
	    coeff = polyfit(order, srcMat, catMat);
	    num_prev = num;
	}
	if (verbose) {
	  //std::cout << "# of objects matched: " << num << std::endl;
	  //for (int i = 0; i < 10; i++) {
	  //printf("%2d %12.5e %12.5e\n", i, coeff[i], coeff[10+i]);
	  //}
	  //printf("\n");
	}
    } else {
	for (unsigned int i = 0; i < srcMat.size(); i++) {
	    matPair.push_back(SourceMatch(catMat.get(i), srcMat.get(i), 0.0));
	}
    }

    return matPair;
}

std::vector<SourceMatch>
hsc::meas::astrom::match(SourceCatalog const &src,
			 SourceCatalog const &cat,
			 int numBrightStars,
			 unsigned int minNumMatchedPair,
			 double matchingAllowanceInPixel,
			 double offsetAllowedInPixel,
			 bool verbose) {
    // Select brightest Nsub stars from list of objects
    // Process both detected from image and external catalog
    int Nsub = numBrightStars;
    SourceCatalog const& srcSub = selectPoint(src, Nsub);
    SourceCatalog const& catSub = selectPoint(cat, srcSub.size()+25);
    //std::cout << srcSub.size() << " " << catSub.size() << std::endl;

    unsigned int catSize = catSub.size();
    unsigned int srcSize = srcSub.size();
    
    /*
    std::ofstream of("zzz");
    for (unsigned int i = 0; i < srcSub.size(); i++) {
	of << srcSub[i].getX() << " " << srcSub[i].getY()<< " " << srcSub[i].getPsfFlux() << " " << catSub[i].getX() << " " << catSub[i].getY() << " " << catSub[i].getPsfFlux()<< std::endl;
    }
    of.close();
    */
    // Construct a list of Pair of objects in catalog
    std::vector<SourcePair> catPair;
    for (size_t i = 0; i < catSize-1; i++) {
	for (size_t j = i+1; j < catSize; j++) {
	    catPair.push_back(SourcePair(catSub.get(i), catSub.get(j)));
	}
    }

    // Sort catPair on distance
    std::sort(catPair.begin(), catPair.end(), cmpPair);

    // Construct a list of Pair of objects in source
    std::vector<SourcePair> srcPair;
    for (size_t i = 0; i < srcSize-1; i++) {
	for (size_t j = i+1; j < srcSize; j++) {
	    srcPair.push_back(SourcePair(srcSub.get(i), srcSub.get(j)));
	}
    }

    // Sort srcPair on distance
    std::sort(srcPair.begin(), srcPair.end(), cmpPair);

    std::vector<SourceMatch> matPair;

    size_t m = 5;		// Number of objects to define the shape
    double e = matchingAllowanceInPixel; // Error allowed for matching
    for (size_t ii = 0; ii < srcPair.size(); ii++) {
	SourcePair p = srcPair[ii];

	    std::vector<SourcePair> q = searchPair(catPair, p, e);

	    // If candidate pairs are found
	    if (q.size() != 0) {

		std::vector<SourcePair> srcMatPair;
		std::vector<SourcePair> catMatPair;

		// Go through candidate pairs
		for (size_t l = 0; l < q.size(); l++) {

		    double dpa = p.pa - q[l].pa;

		    srcMatPair.clear();
		    catMatPair.clear();

		    srcMatPair.push_back(p);
		    catMatPair.push_back(q[l]);

		    //std::cout << "p dist: " << p.distance << " pa: " << p.pa << std::endl;
		    //std::cout << "q dist: " << q[l].distance << " pa: " << q[l].pa << std::endl;

		    for (size_t k = 0; k < srcSub.size(); k++) {
			if (p.first == srcSub.get(k) || p.second == srcSub.get(k)) continue;

			SourcePair pp(p.first, srcSub.get(k));

			std::vector<SourcePair>::iterator r = searchPair3(catPair, pp, q[l], e, dpa);
			if (r != catPair.end()) {
			    srcMatPair.push_back(pp);
			    catMatPair.push_back(*r);
			    //std::cout << "  p dist: " << pp.distance << " pa: " << pp.pa << std::endl;
			    //std::cout << "  r dist: " << (*r).distance << " pa: " << (*r).pa << std::endl;
			    if (srcMatPair.size() == m-1) {
				break;
			    }
			}
		    }

		    bool goodMatch = false;
		    if (srcMatPair.size() == m - 1) {
			for (size_t k = 0; k < srcMatPair.size(); k++) {
			    //std::cout << "  " << srcMatPair[k].first.getX() << " " << srcMatPair[k].first.getY() << " " << srcMatPair[k].second.getX() << " " << srcMatPair[k].second.getY() << std::endl;
			}
			for (size_t k = 0; k < catMatPair.size(); k++) {
			    //std::cout << "  " << catMatPair[k].first.getX() << " " << catMatPair[k].first.getY() << " " << catMatPair[k].second.getX() << " " << catMatPair[k].second.getY() << std::endl;
			}

			goodMatch = true;
			for (size_t k = 1; k < catMatPair.size(); k++) {
			    if (catMatPair[0].first != catMatPair[k].first) {
				goodMatch = false;
			    }
			}

			//std::cout << "  flag: " << goodMatch << std::endl;
		    }

		    if (goodMatch && srcMatPair.size() == m-1) {

			SourceCatalog srcMat(src.getTable());
			SourceCatalog catMat(cat.getTable());

			srcMat.push_back(*srcMatPair[0].first);
			catMat.push_back(*catMatPair[0].first);
			for (size_t k = 0; k < srcMatPair.size(); k++) {
			    srcMat.push_back(*srcMatPair[k].second);
			    catMat.push_back(*catMatPair[k].second);
			}

			double *coeff = polyfit(1, srcMat, catMat);

			//std::cout << coeff[0] << " " << coeff[1] << " " << coeff[2] << std::endl;
			//std::cout << coeff[3] << " " << coeff[4] << " " << coeff[5] << std::endl;
			//std::cout << fabs(coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1.) << std::endl;
			//std::cout << std::endl;

			if (fabs(coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1.) > 0.008 ||
			    fabs(coeff[0]) > offsetAllowedInPixel || fabs(coeff[3]) > offsetAllowedInPixel) {
			    delete [] coeff;
			    continue;
			} else {
			    for (size_t k = 0; k < srcMat.size(); k++) {
				//std::cout << srcMat[k].getX() << " " << srcMat[k].getY() << " " << catMat[k].getX() << " " << catMat[k].getY() << std::endl;
			    }
			    for (size_t k = 1; k < srcMat.size(); k++) {
				//std::cout << "line(" << srcMat[0].getX() << "," << srcMat[0].getY() << "," << srcMat[k].getX() << "," << srcMat[k].getY() << ") # line=0 0" << std::endl;
			    }
			    for (size_t k = 1; k < srcMat.size(); k++) {
				//std::cout << "line(" << catMat[0].getX() << "," << catMat[0].getY() << "," << catMat[k].getX() << "," << catMat[k].getY() << ") # line=0 0 color=red" << std::endl;
			    }
			    matPair = FinalVerify(coeff, src, cat, matchingAllowanceInPixel, verbose);
			    //std::cout << matPair.size() << std::endl;
			    if (matPair.size() <= minNumMatchedPair) {
				delete [] coeff;
				continue;
			    } else {
				//std::cout << "Finish" << std::endl;
				goto FIN;
			    }
			}

		    } // if
		} // for l
		//std::cout << std::endl;
	    }	// if
    }	// for ii

 FIN:

    return matPair;
}
