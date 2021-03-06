#include "fitsio.h"

#include <boost/make_shared.hpp>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>

#include "hsc/meas/astrom/match.h"
#include "hsc/meas/astrom/sipfit.h"
#include "lsst/afw/table/Source.h"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace lsst::afw::table;
namespace afwGeom = lsst::afw::geom;

#include <gsl/gsl_linalg.h>

namespace hsc { namespace meas { namespace astrom { namespace {

double calXi(double a, double d, double A, double D);
double calXi_a(double a, double d, double A, double D);
double calXi_d(double a, double d, double A, double D);
double calXi_A(double a, double d, double A, double D);
double calXi_D(double a, double d, double A, double D);
double calEta(double a, double d, double A, double D);
double calEta_a(double a, double d, double A, double D);
double calEta_d(double a, double d, double A, double D);
double calEta_A(double a, double d, double A, double D);
double calEta_D(double a, double d, double A, double D);

boost::shared_array<double> sipfit(
    int order,
    ProxyVector const &cat,
    ProxyVector const &img,
    SourceTable const & imgTable
) {
    int ncoeff = (order + 1) * (order + 2) / 2 - 1;
    boost::scoped_array<int> xorder(new int[ncoeff]);
    boost::scoped_array<int> yorder(new int[ncoeff]);

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    boost::scoped_array<double> a_data(new double[ncoeff*ncoeff]);
    boost::scoped_array<double> b_data(new double[ncoeff]);
    boost::scoped_array<double> c_data(new double[ncoeff]);
    Key< Covariance<Point<float> > > covKey = imgTable.getCentroidErrKey();
    for (int i = 0; i < ncoeff; i++) {
	for (int j = 0; j < ncoeff; j++) {
	    a_data[i*ncoeff+j] = 0.0;
	    for (unsigned int k = 0; k < img.size(); k++) {
		double w = std::sqrt(img[k].record->get(covKey(0,0)));
		if (w <= 0.0) w = 1.0;
		a_data[i*ncoeff+j] += pow(img[k].getX(), xorder[i]) * 
		                      pow(img[k].getY(), yorder[i]) * 
		                      pow(img[k].getX(), xorder[j]) * 
		                      pow(img[k].getY(), yorder[j]) * w;
	    }
	}
	b_data[i] = c_data[i] = 0.0;
	// Subtract img[k]->getX()
        //          img[k]->getY() to
	// account for Ap, Bp definition of TAN-SIP.
	//     u = U + F(U)
        //     v = V + G(V)
	for (unsigned int k = 0; k < img.size(); k++) {
	    double w = std::sqrt(img[k].record->get(covKey(0,0)));
	    if (w <= 0.0) w = 1.0;
	    b_data[i] += pow(img[k].getX(), xorder[i]) * 
		         pow(img[k].getY(), yorder[i]) * 
		         (cat[k].getX()-img[k].getX()) * w;
	    c_data[i] += pow(img[k].getX(), xorder[i]) * 
		         pow(img[k].getY(), yorder[i]) * 
		         (cat[k].getY()-img[k].getY()) * w;
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
    /*
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, x.get());
    gsl_linalg_cholesky_solve(&a.matrix, &c.vector, y.get());
    */
    boost::shared_array<double> coeff(new double[ncoeff*2]);
    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = x->data[i];
	coeff[i+ncoeff] = y->data[i];
    }

    return coeff;
}

} // anonymous namespace

lsst::afw::image::Wcs::Ptr
fitTANSIP(int order,
          ReferenceMatchVector const &matPair,
          lsst::afw::coord::Coord const &crvalo,
          lsst::afw::geom::PointD const &crpix,
          bool verbose) {
    int npair = matPair.size();
    assert(npair);
    CONST_PTR(SourceTable) srcTable = matPair.front().second->getTable();
    std::vector<int> flag;
    boost::scoped_array<double> x(new double[npair]);
    boost::scoped_array<double> y(new double[npair]);
    boost::scoped_array<double> u(new double[npair]);
    boost::scoped_array<double> v(new double[npair]);

    double ra, dec;

    lsst::afw::geom::PointD crval = crvalo.getPosition(lsst::afw::geom::radians);

    int ncoeff = (order+1)*(order+2)/2 - 1;
    boost::shared_array<double> coeff;
    int ndim = ncoeff * 2 + 2;

    boost::scoped_array<int> xorder(new int[ncoeff]);
    boost::scoped_array<int> yorder(new int[ncoeff]);

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    boost::scoped_array<double> a_data(new double[ndim*ndim]);
    boost::scoped_array<double> b_data(new double[ndim]);

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
	ra = matPair[i].first->getRa().asRadians();
	dec = matPair[i].first->getDec().asRadians();
	double xi    = calXi  (ra, dec, crval[0], crval[1]);
	double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	double eta   = calEta  (ra, dec, crval[0], crval[1]);
	double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	double u = matPair[i].second->getX() - crpix[0];
	double v = matPair[i].second->getY() - crpix[1];
	int i0 = ncoeff * 2 * iexp;
	for (int k = 0; k < ncoeff; k++) {
	    for (int l = 0; l < ncoeff; l++) {
		a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
		                            pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
		a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
		                                          pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
	    }
	    b_data[i0+k] += xi * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	}
	int j0 = ncoeff * 2 * nexp;
	for (int k = 0; k < ncoeff; k++) {
	    a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
	    a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
	    a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
	    a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
	    a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
	}
	a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w1;
	a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w1;
	a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w1;
	a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w1;
	
	b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w1;
	b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w1;
    }

    gsl_matrix_view a = gsl_matrix_view_array(a_data.get(), ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data.get(), ndim);

    boost::shared_ptr<gsl_vector> c(gsl_vector_alloc(ndim), gsl_vector_free);

    //int s;
    /*
      boost::shared_ptr<gsl_permutation> p(gsl_permutation_alloc(ndim), gsl_permutation_free);

      gsl_linalg_LU_decomp(&a.matrix, p.get(), &s);
      gsl_linalg_LU_solve(&a.matrix, p.get(), &b.vector, c.get());
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c.get());

    if (verbose) {
	for (int i = 0; i < ncoeff; i++) {
	    printf("%2d %12.5e %12.5e\n", i, c->data[i], c->data[ncoeff+i]);
	}
	printf("\n");
	printf("   %12.5e %12.5e\n", c->data[ncoeff*2], c->data[ncoeff*2+1]);
	printf("\n");
    }

    crval[0] += c->data[ncoeff*2];
    crval[1] += c->data[ncoeff*2+1];

    Eigen::Matrix2d cd; cd << c->data[0], c->data[1], c->data[ncoeff], c->data[ncoeff+1];
    double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
    
    Eigen::MatrixXd sipA = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipB = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 2; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipA(i,j) = ( cd(1,1)*c->data[n+j] - cd(0,1)*c->data[ncoeff+n+j]) / D;
	    sipB(i,j) = (-cd(1,0)*c->data[n+j] + cd(0,0)*c->data[ncoeff+n+j]) / D;
	}
    }
    if (verbose) {
	std::cout << "sipA" << std::endl;
	std::cout << sipA << std::endl << std::endl;
	std::cout << "sipB" << std::endl;
	std::cout << sipB << std::endl << std::endl;
    }

    ProxyVector cat;
    ProxyVector img;
    for (int i = 0; i < npair; i++) {
        ra  = matPair[i].first->getRa().asRadians();
        dec = matPair[i].first->getDec().asRadians();
        x[i] = calXi(ra, dec, crval[0], crval[1]);
        y[i] = calEta(ra, dec, crval[0], crval[1]);
        double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
        RecordProxy s(
            matPair[i].first,
            lsst::afw::geom::Point2D((cd(1,1)*x[i]-cd(0,1)*y[i])/D,
                                     (-cd(1,0)*x[i]+cd(0,0)*y[i])/D)
        );
        cat.push_back(s);
        
        u[i] = matPair[i].second->getX() - crpix[0];
        v[i] = matPair[i].second->getY() - crpix[1];
        RecordProxy s2(matPair[i].second, lsst::afw::geom::Point2D(u[i], v[i]));
        img.push_back(s2);
    }
    coeff = sipfit(order, img, cat, *srcTable);
    if (verbose) {
	for (int i = 0; i < ncoeff-2; i++) {
	    printf("%2d %12.5e %12.5e\n", i, coeff[i], coeff[ncoeff-2+i]);
	}
	printf("\n");
    }

    Eigen::MatrixXd sipAp = Eigen::MatrixXd::Zero(order+1,order+1);
    Eigen::MatrixXd sipBp = Eigen::MatrixXd::Zero(order+1,order+1);
    for (int k = 1; k <= order; k++) {
	for (int i = k; i >= 0; i--) {
	    int j = k - i;
	    int n = k*(k+1)/2 - 1;
	    sipAp(i,j) = coeff[n+j];
	    sipBp(i,j) = coeff[ncoeff+n+j];
	}
    }

    if (verbose) {
	std::cout << "sipAp" << std::endl;
	std::cout << sipAp << std::endl << std::endl;
	std::cout << "sipBp" << std::endl;
	std::cout << sipBp << std::endl << std::endl;
    }

    crval[0] *= R2D;
    crval[1] *= R2D;
    cd *= R2D;
    return boost::make_shared<lsst::afw::image::TanWcs>(crval, crpix, cd, sipA, sipB, sipAp, sipBp);
}

lsst::afw::image::Wcs::Ptr
fitTAN(ReferenceMatchVector const &matPair,
       bool verbose) {
    int npair = matPair.size();
    ProxyVector img;
    ProxyVector cat;
    boost::scoped_array<double> x(new double[npair]);
    boost::scoped_array<double> y(new double[npair]);
    boost::scoped_array<double> u(new double[npair]);
    boost::scoped_array<double> v(new double[npair]);

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    for (int i = 0; i < npair; i++) {
	lsst::afw::geom::Point3D v = matPair[i].first->getCoord().getVector();
        cx += v[0];
	cy += v[1];
	cz += v[2];
        Sx += matPair[i].second->getX();
        Sy += matPair[i].second->getY();
    }
    cx /= npair;
    cy /= npair;
    cz /= npair;
    lsst::afw::coord::Coord cmean(lsst::afw::geom::Point3D(cx, cy, cz));
    lsst::afw::geom::PointD crval = cmean.getPosition(lsst::afw::geom::radians);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::Point2D(Sx/npair, Sy/npair);

    int order = 1;
    int ncoeff = (order+1)*(order+2)/2 - 1;
    int ndim = ncoeff * 2 + 2;

    boost::scoped_array<int> xorder(new int[ncoeff]);
    boost::scoped_array<int> yorder(new int[ncoeff]);

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    boost::scoped_array<double> a_data(new double[ndim*ndim]);
    boost::scoped_array<double> b_data(new double[ndim]);

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
        double ra = matPair[i].first->getRa().asRadians();
        double dec = matPair[i].first->getDec().asRadians();
        double xi    = calXi  (ra, dec, crval[0], crval[1]);
        double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
        double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
        double eta   = calEta  (ra, dec, crval[0], crval[1]);
        double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
        double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
        double u = matPair[i].second->getX() - crpix[0];
        double v = matPair[i].second->getY() - crpix[1];
        if (verbose) {
            std::cout << u << "," << v << " --> " << xi << "," << eta << std::endl;
        }
        int i0 = ncoeff * 2 * iexp;
        for (int k = 0; k < ncoeff; k++) {
            for (int l = 0; l < ncoeff; l++) {
                a_data[(i0+k)*ndim+i0+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
                    pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
                a_data[(i0+ncoeff+k)*ndim+i0+ncoeff+l] += pow(u, xorder[k]) * pow(v, yorder[k]) *
                    pow(u, xorder[l]) * pow(v, yorder[l]) * w1;
            }
            b_data[i0+k] += xi * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            b_data[i0+ncoeff+k] += eta * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
        }
        int j0 = ncoeff * 2 * nexp;
        for (int k = 0; k < ncoeff; k++) {
            a_data[(i0+k)*ndim+j0+iexp*2]          += -xi_A  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+k)*ndim+j0+iexp*2+1]        += -xi_D  * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+ncoeff+k)*ndim+j0+iexp*2]   += -eta_A * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1] += -eta_D * pow(u, xorder[k]) * pow(v, yorder[k]) * w1;
            a_data[(j0+iexp*2)*ndim+i0+k]          = a_data[(i0+k)*ndim+j0+iexp*2];
            a_data[(j0+iexp*2+1)*ndim+i0+k]        = a_data[(i0+k)*ndim+j0+iexp*2+1];
            a_data[(j0+iexp*2)*ndim+i0+ncoeff+k]   = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2];
            a_data[(j0+iexp*2+1)*ndim+i0+ncoeff+k] = a_data[(i0+ncoeff+k)*ndim+j0+iexp*2+1];
        }
        a_data[(j0+iexp*2)*ndim+j0+iexp*2]     += (xi_A * xi_A + eta_A * eta_A) * w1;
        a_data[(j0+iexp*2)*ndim+j0+iexp*2+1]   += (xi_D * xi_A + eta_D * eta_A) * w1;
        a_data[(j0+iexp*2+1)*ndim+j0+iexp*2]   += (xi_A * xi_D + eta_A * eta_D) * w1;
        a_data[(j0+iexp*2+1)*ndim+j0+iexp*2+1] += (xi_D * xi_D + eta_D * eta_D) * w1;
        
        b_data[j0+iexp*2]   += -(xi * xi_A + eta * eta_A) * w1;
        b_data[j0+iexp*2+1] += -(xi * xi_D + eta * eta_D) * w1;
    }
    
    gsl_matrix_view a = gsl_matrix_view_array(a_data.get(), ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data.get(), ndim);
    
    boost::shared_ptr<gsl_vector> c(gsl_vector_alloc(ndim), gsl_vector_free);

    if (verbose) {
        for (int i = 0; i < ndim; ++i) {
            for (int j = 0; j < ndim; ++j) {
                std::cout << a_data[i*ndim+j] << " ";
            }
            std::cout << std::endl;
        }
        for (int i = 0; i < ndim; ++i) {
            std::cout << b_data[i] << " ";
        }
        std::cout << std::endl;
    }

    //int s;
    /*
    boost::shared_ptr<gsl_permutation> p(gsl_permutation_alloc(ndim), gsl_permutation_free);

    gsl_linalg_LU_decomp(&a.matrix, p.get(), &s);
    gsl_linalg_LU_solve(&a.matrix, p.get(), &b.vector, c.get());
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c.get());

    if (verbose) {
	for (int i = 0; i < ncoeff; i++) {
	    printf("%2d %12.5e %12.5e\n", i, c->data[i], c->data[ncoeff+i]);
	}
	printf("\n");
	printf("   %12.5e %12.5e\n", c->data[ncoeff*2], c->data[ncoeff*2+1]);
	printf("\n");
    }

    crval[0] += c->data[ncoeff*2];
    crval[1] += c->data[ncoeff*2+1];

    Eigen::Matrix2d cd; cd << c->data[0], c->data[1], c->data[ncoeff], c->data[ncoeff+1];
    
    crval[0] *= R2D;
    crval[1] *= R2D;
    cd *= R2D;
    return boost::make_shared<lsst::afw::image::TanWcs>(crval, crpix, cd);
}

namespace {

double calXi(double a, double d, double A, double D) {
    return cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_a(double a, double d, double A, double D) {
    return cos(D)*pow(cos(d),2.)*pow(sin(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_d(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*cos(d)*sin(a-A)*sin(a-A)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -cos(d)*cos(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calXi_D(double a, double d, double A, double D) {
    return -cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.);
}

double calEta(double a, double d, double A, double D) {
    return (cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_a(double a, double d, double A, double D) {
    return cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	  +sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_d(double a, double d, double A, double D) {
    return -(sin(D)*cos(d)-cos(D)*sin(d)*cos(a-A))*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   +(cos(D)*cos(d)+sin(D)*sin(d)*cos(a-A))/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_A(double a, double d, double A, double D) {
    return -cos(D)*cos(d)*sin(a-A)*(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A))/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)
	   -sin(D)*cos(d)*sin(a-A)/(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A));
}

double calEta_D(double a, double d, double A, double D) {
    return -pow(cos(D)*sin(d)-sin(D)*cos(d)*cos(a-A),2.)/pow(sin(D)*sin(d)+cos(D)*cos(d)*cos(a-A),2.)-1.;
}

}}}} // namespace hsc::meas::astrom::<anonymous>
