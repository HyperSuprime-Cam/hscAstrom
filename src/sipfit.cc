#include "fitsio.h"

#include "hsc/meas/astrom/sipfit.h"
#include "lsst/afw/detection/Source.h"

#define D2R (M_PI/180.)
#define R2D (180./M_PI)

using namespace hsc::meas::astrom;
using namespace lsst::afw::detection;

#include <gsl/gsl_linalg.h>

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

double *sipfit(int order,
		SourceSet const &img,
		SourceSet const &cat) {
    int ncoeff = (order + 1) * (order + 2) / 2 - 1;
    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ncoeff*ncoeff];
    double *b_data = new double[ncoeff];
    double *c_data = new double[ncoeff];

    for (int i = 0; i < ncoeff; i++) {
	for (int j = 0; j < ncoeff; j++) {
	    a_data[i*ncoeff+j] = 0.0;
	    for (unsigned int k = 0; k < img.size(); k++) {
		double w = img[k]->getXAstromErr();
		if (w <= 0.0) w = 1.0;
		a_data[i*ncoeff+j] += pow(img[k]->getXAstrom(), xorder[i]) * 
		                      pow(img[k]->getYAstrom(), yorder[i]) * 
		                      pow(img[k]->getXAstrom(), xorder[j]) * 
		                      pow(img[k]->getYAstrom(), yorder[j]) * w;
	    }
	}
	b_data[i] = c_data[i] = 0.0;
	// Subtract img[k]->getXAstrom()
        //          img[k]->getYAstrom() to
	// account for Ap, Bp definition of TAN-SIP.
	//     u = U + F(U)
        //     v = V + G(V)
	for (unsigned int k = 0; k < img.size(); k++) {
	    double w = img[k]->getXAstromErr();
	    if (w <= 0.0) w = 1.0;
	    b_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		         pow(img[k]->getYAstrom(), yorder[i]) * 
		         (cat[k]->getXAstrom()-img[k]->getXAstrom()) * w;
	    c_data[i] += pow(img[k]->getXAstrom(), xorder[i]) * 
		         pow(img[k]->getYAstrom(), yorder[i]) * 
		         (cat[k]->getYAstrom()-img[k]->getYAstrom()) * w;
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
    /*
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, x);
    gsl_linalg_cholesky_solve(&a.matrix, &c.vector, y);
    */
    double *coeff = new double[ncoeff*2];
    for (int i = 0; i < ncoeff; i++) {
	coeff[i] = x->data[i];
	coeff[i+ncoeff] = y->data[i];
    }

    gsl_permutation_free(p);
    gsl_vector_free(x);
    gsl_vector_free(y);
  
    delete [] xorder;
    delete [] yorder;
    delete [] a_data;
    delete [] b_data;
    delete [] c_data;

    return coeff;
}

lsst::afw::image::Wcs::Ptr
hsc::meas::astrom::fitTANSIP(int order,
			  std::vector<SourceMatch> const &matPair,
			  lsst::afw::geom::PointD &crvalo,
			  lsst::afw::geom::PointD &crpixo,
			  bool verbose) {
    int npair = matPair.size();
    SourceSet img;
    SourceSet cat;
    std::vector<int> flag;
    double *x = new double[npair];
    double *y = new double[npair];
    double *u = new double[npair];
    double *v = new double[npair];

    double ra, dec;

    lsst::afw::geom::PointD crpix = crpixo;
    lsst::afw::geom::PointD crval = crvalo;

    crval[0] *= D2R;
    crval[1] *= D2R;

    int ncoeff = (order+1)*(order+2)/2 - 1;
    double *coeff = NULL;
    int ndim = ncoeff * 2 + 2;

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
	ra = matPair[i].first->getRa() * D2R;
	dec = matPair[i].first->getDec() * D2R;
	double xi    = calXi  (ra, dec, crval[0], crval[1]);
	double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	double eta   = calEta  (ra, dec, crval[0], crval[1]);
	double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	double u = matPair[i].second->getXAstrom() - crpix[0];
	double v = matPair[i].second->getYAstrom() - crpix[1];
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

    gsl_matrix_view a = gsl_matrix_view_array(a_data, ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data, ndim);

    gsl_vector *c = gsl_vector_alloc(ndim);

    //int s;
    /*
    gsl_permutation *p = gsl_permutation_alloc(ndim);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, c);
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c);

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

    cat.clear();
    for (int i = 0; i < npair; i++) {
	ra  = matPair[i].first->getRa() * D2R;
	dec = matPair[i].first->getDec() * D2R;
	x[i] = calXi(ra, dec, crval[0], crval[1]);
	y[i] = calEta(ra, dec, crval[0], crval[1]);
	double D = cd(0,0) * cd(1,1) - cd(0,1) * cd(1,0);
	Source::Ptr s = Source::Ptr(new Source());
	s->setXAstrom(( cd(1,1)*x[i]-cd(0,1)*y[i])/D);
	s->setYAstrom((-cd(1,0)*x[i]+cd(0,0)*y[i])/D);
	cat.push_back(s);

	u[i] = matPair[i].second->getXAstrom() - crpix[0];
	v[i] = matPair[i].second->getYAstrom() - crpix[1];
	Source::Ptr s2 = Source::Ptr(new Source());
	s2->setXAstrom(u[i]);
	s2->setYAstrom(v[i]);
	img.push_back(s2);
    }
    coeff = sipfit(order, cat, img);
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
    lsst::afw::image::TanWcs wcs(crval, crpix, cd, sipA, sipB, sipAp, sipBp);

    delete [] xorder;
    delete [] yorder;
    delete [] x;
    delete [] y;
    delete [] u;
    delete [] v;

    return wcs.clone();
}

lsst::afw::image::Wcs::Ptr
hsc::meas::astrom::fitTAN(std::vector<SourceMatch> const &matPair,
		       bool verbose) {
    int npair = matPair.size();
    SourceSet img;
    SourceSet cat;
    double *x = new double[npair];
    double *y = new double[npair];
    double *u = new double[npair];
    double *v = new double[npair];

    double Sra = 0.0;
    double Sdec = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;
    for (int i = 0; i < npair; i++) {
	Sra  += matPair[i].first->getRa();
	Sdec += matPair[i].first->getDec();
	Sx += matPair[i].second->getXAstrom();
	Sy += matPair[i].second->getYAstrom();
    }
    lsst::afw::geom::PointD crval = lsst::afw::geom::Point2D(Sra/npair, Sdec/npair);
    lsst::afw::geom::PointD crpix = lsst::afw::geom::Point2D(Sx/npair, Sy/npair);

    crval[0] *= D2R;
    crval[1] *= D2R;

    int order = 1;
    int ncoeff = (order+1)*(order+2)/2 - 1;
    int ndim = ncoeff * 2 + 2;

    int *xorder = new int[ncoeff];
    int *yorder = new int[ncoeff];

    int n = 0;
    for (int i = 1; i <= order; i++) {
	for (int k = 0; k <= i; k++) {
	    int j = i - k;
	    xorder[n] = j;
	    yorder[n] = k;
	    n++;
	}
    }

    double *a_data = new double[ndim*ndim];
    double *b_data = new double[ndim];

    for (int i = 0; i < ndim; i++) {
	for (int j = 0; j < ndim; j++) {
	    a_data[i*ndim+j] = 0.0;
	}
	b_data[i] = 0.0;
    }

    int iexp = 0; int nexp = 1;
    double w1 = 1.0;
    for (int i = 0; i < npair; i++) {
	double ra = matPair[i].first->getRa() * D2R;
	double dec = matPair[i].first->getDec() * D2R;
	double xi    = calXi  (ra, dec, crval[0], crval[1]);
	double xi_A  = calXi_A(ra, dec, crval[0], crval[1]);
	double xi_D  = calXi_D(ra, dec, crval[0], crval[1]);
	double eta   = calEta  (ra, dec, crval[0], crval[1]);
	double eta_A = calEta_A(ra, dec, crval[0], crval[1]);
	double eta_D = calEta_D(ra, dec, crval[0], crval[1]);
	double u = matPair[i].second->getXAstrom() - crpix[0];
	double v = matPair[i].second->getYAstrom() - crpix[1];
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

    gsl_matrix_view a = gsl_matrix_view_array(a_data, ndim, ndim);
    gsl_vector_view b = gsl_vector_view_array(b_data, ndim);

    gsl_vector *c = gsl_vector_alloc(ndim);

    /*
    gsl_permutation *p = gsl_permutation_alloc(ndim);

    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &b.vector, c);
    */
    gsl_linalg_cholesky_decomp(&a.matrix);
    gsl_linalg_cholesky_solve(&a.matrix, &b.vector, c);

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
    lsst::afw::image::TanWcs wcs(crval, crpix, cd);

    delete [] xorder;
    delete [] yorder;
    delete [] x;
    delete [] y;
    delete [] u;
    delete [] v;

    return wcs.clone();
}

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
