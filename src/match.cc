#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include "hsc/meas/astrom/match.h"
#include "lsst/utils/ieee.h"
#include "lsst/afw/image/Wcs.h"

#include "boost/scoped_array.hpp"
#include "boost/shared_array.hpp"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::astrom;
using namespace lsst::afw::table;

namespace {

// Algorithm is based on V.Tabur 2007, PASA, 24, 189-198
// "Fast Algorithms for Matching CCD Images to a Stellar Catalogue"

bool cmpPair(ProxyPair const &a, ProxyPair const &b) {
    return a.distance > b.distance;
}

// Compair Source based on its PsfFlux
// Ordering is bright to faint
struct CompareProxyFlux {

    bool operator()(RecordProxy const & a, RecordProxy const & b) const {
        double aFlux = a.record->get(key);
        double bFlux = b.record->get(key);
        if (lsst::utils::isnan(aFlux)) {
            aFlux = 0.0;
        }
        if (lsst::utils::isnan(bFlux)) { 
            bFlux = 0.0;
        }
        return aFlux > bFlux;
    }

    Key<double> key;
};

ProxyVector makeVector(SourceCatalog const & sources) {
    ProxyVector r; r.reserve(sources.size());
    for (SourceCatalog::const_iterator i = sources.begin(); i != sources.end(); ++i) {
        r.push_back(RecordProxy(i, i->getCentroid()));
    }
    return r;
}

ProxyVector makeVector(SimpleCatalog const & refs, lsst::afw::image::Wcs const & wcs) {
    ProxyVector r; r.reserve(refs.size());
    for (SimpleCatalog::const_iterator i = refs.begin(); i != refs.end(); ++i) {
        r.push_back(RecordProxy(i, wcs.skyToPixel(i->getCoord())));
    }
    return r;    
}

ProxyVector selectPoint(
    ProxyVector const &a, Key<double> const & key, std::size_t num, std::size_t start=0
) {
    // copy and sort array of pointers on apFlux
    CompareProxyFlux cmp = {key};
    ProxyVector b(a);
    std::sort(b.begin(), b.end(), cmp);
    std::size_t end = std::min(start + num, b.size());
    return ProxyVector(b.begin() + start, b.begin() + end);
}

std::vector<ProxyPair> searchPair(std::vector<ProxyPair> &a, ProxyPair &p, double &e, double &e_dpa) {
    std::vector<ProxyPair> v;

    for (size_t i = 0; i < a.size(); i++) {
	double dd = fabs(a[i].distance - p.distance);
	double dpa = fabs(a[i].pa - p.pa);
	if (dd < e && dpa < e_dpa) {
//	if (dd < e) {
	    a[i].deltaD = dd;
	    v.push_back(a[i]);
	}
    }

    return v;
}

std::vector<ProxyPair>::iterator searchPair3(std::vector<ProxyPair> &a,
                                             const ProxyPair &p,
                                             const ProxyPair &q,
                                             const double &e,
                                             const double &dpa,
                                             const double &e_dpa = 0.02) {
    std::vector<ProxyPair>::iterator idx = a.end();
    double dd_min = 1.E+10;
    //double dpa_min = e_dpa;

    for (std::vector<ProxyPair>::iterator i = a.begin(); i < a.end(); i++) {
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

void transform(
    int order, boost::shared_array<double> const & coeff, double x, double y, double *xn, double *yn
) {
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

boost::shared_array<double> polyfit(int order, ProxyVector const &img, ProxyVector const &cat) {
    int ncoeff = (order + 1) * (order + 2) / 2;
    boost::scoped_array<int> xorder(new int[ncoeff]);
    boost::scoped_array<int> yorder(new int[ncoeff]);

    int n = 0;
    for (int i = 0; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    boost::scoped_array<int> flag(new int[img.size()]);
    for (size_t k = 0; k < img.size(); k++) {
	flag[k] = 1;
    }

    boost::scoped_array<double> a_data(new double[ncoeff*ncoeff]);
    boost::scoped_array<double> b_data(new double[ncoeff]);
    boost::scoped_array<double> c_data(new double[ncoeff]);

    boost::shared_array<double> coeff(new double[ncoeff*2]);

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

	gsl_matrix_view a = gsl_matrix_view_array(a_data.get(), ncoeff, ncoeff);
	gsl_vector_view b = gsl_vector_view_array(b_data.get(), ncoeff);
	gsl_vector_view c = gsl_vector_view_array(c_data.get(), ncoeff);

	boost::shared_ptr<gsl_vector> x(gsl_vector_alloc(ncoeff), gsl_vector_free);
	boost::shared_ptr<gsl_vector> y(gsl_vector_alloc(ncoeff), gsl_vector_free);

	int s;

	boost::shared_ptr<gsl_permutation> p(gsl_permutation_alloc(ncoeff), gsl_permutation_free);

	gsl_linalg_LU_decomp(&a.matrix, p.get(), &s);
	gsl_linalg_LU_solve(&a.matrix, p.get(), &b.vector, x.get());
	gsl_linalg_LU_solve(&a.matrix, p.get(), &c.vector, y.get());

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

    }

    return coeff;
}

ProxyVector::const_iterator searchNearestPoint(ProxyVector const &cat, double& x, double& y, double& e) {
    ProxyVector::const_iterator index = cat.end();
    double d_min, d;

    d_min = e;

    for (ProxyVector::const_iterator i = cat.begin(); i < cat.end(); i++) {
	d = sqrt(pow((i->getX()-x), 2) + pow((i->getY()-y), 2));
	if (d < d_min) {
	    index = i;
	    d_min = d;
	    break;
	}
    }

    return index;
}

ReferenceMatchVector FinalVerify(
    boost::shared_array<double> coeff,
    ProxyVector const & cat,
    ProxyVector const & src,
    double matchingAllowanceInPixel,
    bool verbose
) {
    ProxyVector srcMat;
    ProxyVector catMat;
    ReferenceMatchVector matPair;

    double x0, y0, x1, y1;
    int num = 0, num_prev = -1;
    for (size_t i = 0; i < src.size(); i++) {
        x0 = src[i].getX();
        y0 = src[i].getY();
        transform(1, coeff, x0, y0, &x1, &y1);
        double e = matchingAllowanceInPixel;
        ProxyVector::const_iterator p = searchNearestPoint(cat, x1, y1, e);
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
                ProxyVector::const_iterator p = searchNearestPoint(cat, x1, y1, e);
                if (p != cat.end()) {
                    num++;
                    srcMat.push_back(src[i]);
                    catMat.push_back(*p);
                    matPair.push_back(
                        ReferenceMatch(*p, boost::static_pointer_cast<SourceRecord>(src[i].record), 0.0)
                    );
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
            matPair.push_back(
                ReferenceMatch(catMat[i], boost::static_pointer_cast<SourceRecord>(srcMat[i].record), 0.0)
            );
        }
    }

    return matPair;
}

} // anonymous namespace


ProxyVector hsc::meas::astrom::makeProxies(SourceCatalog const & sources) {
    ProxyVector r;
    r.reserve(sources.size());
    for (SourceCatalog::const_iterator i = sources.begin(); i != sources.end(); ++i) {
        r.push_back(RecordProxy(i, i->getCentroid()));
    }
    return r;
}

ProxyVector hsc::meas::astrom::makeProxies(SimpleCatalog const & refs, lsst::afw::image::Wcs const & wcs) {
    ProxyVector r;
    r.reserve(refs.size());
    for (SimpleCatalog::const_iterator i = refs.begin(); i != refs.end(); ++i) {
        r.push_back(RecordProxy(i, wcs.skyToPixel(i->getCoord())));
    }
    return r;
}


ReferenceMatchVector
hsc::meas::astrom::match(
    SimpleCatalog const &cat,
    SourceCatalog const &src,
    lsst::afw::image::Wcs const & wcs,
    int numBrightStars,
    unsigned int minNumMatchedPair,
    double matchingAllowanceInPixel,
    int catOffset,
    double offsetAllowedInPixel,
    double rotationAllowedInRad,
    double angleDiffFrom90,
    bool verbose
) {
    // Select brightest Nsub stars from list of objects
    // Process both detected from image and external catalog
    int Nsub = numBrightStars;
    ProxyVector proxyCat = makeProxies(cat, wcs);
    ProxyVector proxySrc = makeProxies(src);
    ProxyVector srcSub = selectPoint(proxySrc, src.getTable()->getApFluxKey(), Nsub);
    ProxyVector catSub = selectPoint(proxyCat, cat.getSchema().find<double>("flux").key, srcSub.size()+25, catOffset);
    if (verbose)
	std::cout << "Catalog sizes: " << srcSub.size() << " " << catSub.size() << std::endl;

    unsigned int catSize = catSub.size();
    unsigned int srcSize = srcSub.size();
/*    
    CompareProxyFlux cmp = {cat.getSchema().find<double>("flux").key};
    std::sort(proxyCat.begin(), proxyCat.end(), cmp);
    cmp = {src.getTable()->getPsfFluxKey()};
    std::sort(proxySrc.begin(), proxySrc.end(), cmp);
*/
    /*
    std::ofstream of("zzz");
    for (unsigned int i = 0; i < srcSub.size(); i++) {
	of << srcSub[i].getX() << " " << srcSub[i].getY()<< " " << catSub[i].getX() << " " << catSub[i].getY() << std::endl;
    }
    of.close();
    */
    // Construct a list of Pair of objects in catalog
    std::vector<ProxyPair> catPair;
    for (size_t i = 0; i < catSize-1; i++) {
        for (size_t j = i+1; j < catSize; j++) {
            catPair.push_back(ProxyPair(catSub[i], catSub[j]));
        }
    }

    // Sort catPair on distance
    std::sort(catPair.begin(), catPair.end(), cmpPair);

    // Construct a list of Pair of objects in source
    std::vector<ProxyPair> srcPair;
    for (size_t i = 0; i < srcSize-1; i++) {
        for (size_t j = i+1; j < srcSize; j++) {
            srcPair.push_back(ProxyPair(srcSub[i], srcSub[j]));
        }
    }

    // Sort srcPair on distance
    std::sort(srcPair.begin(), srcPair.end(), cmpPair);

    ReferenceMatchVector matPair;
    ReferenceMatchVector matPairSave;
    std::vector<ReferenceMatchVector> matPairCand;

    size_t m = 6;		// Number of objects to define the shape
    double e = matchingAllowanceInPixel; // Error allowed for matching
    double e_dpa = rotationAllowedInRad;
    for (size_t ii = 0; ii < srcPair.size(); ii++) {
	ProxyPair p = srcPair[ii];

	std::vector<ProxyPair> q = searchPair(catPair, p, e, e_dpa);

	// If candidate pairs are found
	if (q.size() != 0) {

	    std::vector<ProxyPair> srcMatPair;
	    std::vector<ProxyPair> catMatPair;

	    // Go through candidate pairs
	    for (size_t l = 0; l < q.size(); l++) {

		double dpa = p.pa - q[l].pa;

		srcMatPair.clear();
		catMatPair.clear();

		srcMatPair.push_back(p);
		catMatPair.push_back(q[l]);

		if (verbose) {
		    std::cout << "p dist: " << p.distance << " pa: " << p.pa << std::endl;
		    std::cout << "q dist: " << q[l].distance << " pa: " << q[l].pa << std::endl;
		}

		for (size_t k = 0; k < srcSub.size(); k++) {
		    if (p.first == srcSub[k] || p.second == srcSub[k]) continue;

		    ProxyPair pp(p.first, srcSub[k]);
                
		    std::vector<ProxyPair>::iterator r = searchPair3(catPair, pp, q[l], e, dpa, e_dpa);
		    if (r != catPair.end()) {
			srcMatPair.push_back(pp);
			catMatPair.push_back(*r);
			if (verbose) {
			    std::cout << "  p dist: " << pp.distance << " pa: " << pp.pa << std::endl;
			    std::cout << "  r dist: " << (*r).distance << " pa: " << (*r).pa << std::endl;
			}
			if (srcMatPair.size() == m-1) {
			    break;
			}
		    }
		}

		bool goodMatch = false;
		if (srcMatPair.size() == m - 1) {
		    goodMatch = true;
		    for (size_t k = 1; k < catMatPair.size(); k++) {
			if (catMatPair[0].first != catMatPair[k].first) {
			    goodMatch = false;
			}
		    }
		}

		if (goodMatch && srcMatPair.size() == m-1) {

		    ProxyVector srcMat;
		    ProxyVector catMat;

		    srcMat.push_back(srcMatPair[0].first);
		    catMat.push_back(catMatPair[0].first);
		    for (size_t k = 0; k < srcMatPair.size(); k++) {
			srcMat.push_back(srcMatPair[k].second);
			catMat.push_back(catMatPair[k].second);
		    }

		    boost::shared_array<double> coeff = polyfit(1, srcMat, catMat);

		    if (verbose) {
			for (size_t k = 0; k < srcMat.size(); k++) {
			    std::cout << "circle(" << srcMat[k].getX() << "," << srcMat[k].getY() << ",10) # color=green" << std::endl;
			    std::cout << "circle(" << catMat[k].getX() << "," << catMat[k].getY() << ",10) # color=red" << std::endl;
			    std::cout << "line(" << srcMat[0].getX() << "," << srcMat[0].getY() << "," << srcMat[k].getX() << "," << srcMat[k].getY() << ") # line=0 0 color=green" << std::endl;
			    std::cout << "line(" << catMat[0].getX() << "," << catMat[0].getY() << "," << catMat[k].getX() << "," << catMat[k].getY() << ") # line=0 0 color=red" << std::endl;
			}
		    }

		    double a = coeff[1];
		    double b = coeff[2];
		    double c = coeff[4];
		    double d = coeff[5];
		    double theta = acos((a*b+c*d)/(sqrt(a*a+c*c)*sqrt(b*b+d*d))) / M_PI * 180.;
		    if (verbose) {
                        std::cout << "Linear fit from match:" << std::endl;
			std::cout << coeff[0] << " " << coeff[1] << " " << coeff[2] << std::endl;
			std::cout << coeff[3] << " " << coeff[4] << " " << coeff[5] << std::endl;
			std::cout << coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1. << std::endl;
			std::cout << theta << std::endl;
		    }
		    if (((fabs(coeff[1] * coeff[5] - coeff[2] * coeff[4] - 1.) > 0.02 || fabs(theta - 90.) > 0.25) &&
			 fabs(theta - 90.) > angleDiffFrom90) ||
			fabs(coeff[0]) > offsetAllowedInPixel || fabs(coeff[3]) > offsetAllowedInPixel) {
			if (verbose)
			    std::cout << "Bad; continuing" << std::endl;
			continue;
		    } else {

			double x0, y0, x1, y1;
			int num = 0;
			srcMat.clear();
			catMat.clear();
			for (size_t i = 0; i < srcSub.size(); i++) {
			    x0 = srcSub[i].getX();
			    y0 = srcSub[i].getY();
			    transform(1, coeff, x0, y0, &x1, &y1);
			    ProxyVector::const_iterator p = searchNearestPoint(catSub, x1, y1, e);
			    if (p != catSub.end()) {
				num++;
				srcMat.push_back(srcSub[i]);
				catMat.push_back(*p);
                                if (verbose) {
                                    std::cout << "Match: " << x0 << "," << y0 << " --> " << x1 << "," << y1 <<
                                        " <==> " << p->getX() << "," << p->getY() << std::endl;
                                }
			    }
			}
                        if (num <= 5) {
                            // Can get matrix = 0,0,0,0; everything matches a single catalog object
                            if (verbose) {
                                std::cout << "Insufficient initial matches; continuing" << std::endl;
                            }
                            continue;
                        }
			coeff = polyfit(1, srcMat, catMat);
                        if (verbose) {
                            std::cout << "Coefficients from initial matching:" << std::endl;
                            for (size_t i = 0; i < 6; ++i) {
                                std::cout << coeff[i] << " ";
                            }
                            std::cout << std::endl;
                        }

			matPair = FinalVerify(coeff, proxyCat, proxySrc, matchingAllowanceInPixel, verbose);
			/*
			for (int k = 0; k < matPair.size(); k++) {
			    double flux_cat = matPair[k].first->get(matPair[k].first->getSchema().find<double>("flux").key);
			    double flux_src = matPair[k].second->get(matPair[k].second->getTable()->getPsfFluxKey());
			    double delta_mag = 2.5 * log10(flux_src/flux_cat);
			    std::cout << delta_mag  << std::endl;
			}
			*/
			if (verbose)
			    std::cout << "Number of matches: " << matPair.size() << " vs " <<
                                minNumMatchedPair << std::endl;
			if (matPair.size() <= minNumMatchedPair) {
			    if (verbose)
				std::cout << "Insufficient final matches; continuing" << std::endl;
			    if (matPair.size() > matPairSave.size()) {
			        matPairSave = matPair;
			    }
			    continue;
			} else {
			    if (verbose)
				std::cout << "Finish" << std::endl;
                            matPairCand.push_back(matPair);
                            if (matPairCand.size() == 3)
                                goto END;
			}
		    }
                
		} // if
	    } // for l
	} // if
    } // for ii

 END:
    if (matPairCand.size() == 0) {
        return matPairSave;
    } else {
        size_t nmatch = matPairCand[0].size();
        ReferenceMatchVector matPairRet = matPairCand[0];
        for (size_t i = 1; i < matPairCand.size(); i++) {
            if (matPairCand[i].size() > nmatch) {
                nmatch = matPairCand[i].size();
                matPairRet = matPairCand[i];
            }
        }
        return matPairRet;
    }
}
