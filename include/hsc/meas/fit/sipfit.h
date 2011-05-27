// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_FIT_H)
#define HSC_MEAS_FIT_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/detection/SourceMatch.h"

namespace hsc {
    namespace meas {
        namespace fit {

	    lsst::afw::image::Wcs::Ptr
	    fitTANSIP(int order,
		      std::vector<lsst::afw::detection::SourceMatch> const &matPair,
		      lsst::afw::geom::PointD &crval,
		      lsst::afw::geom::PointD &crpix,
		      bool verbose = false);

	}
    }
}

#endif
