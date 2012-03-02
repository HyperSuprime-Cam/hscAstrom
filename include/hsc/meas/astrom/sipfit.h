// -*- lsst-c++ -*-
#if !defined(HSC_MEAS_FIT_H)
#define HSC_MEAS_FIT_H

#include <vector>
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/table/Match.h"

namespace hsc {
    namespace meas {
        namespace astrom {

	    lsst::afw::image::Wcs::Ptr
	      fitTANSIP(int order,
			lsst::afw::table::SourceMatchVector const &matPair,
			lsst::afw::coord::Coord::Ptr &crval,
			lsst::afw::geom::PointD &crpix,
			bool verbose = false);

	    lsst::afw::image::Wcs::Ptr
	      fitTAN(lsst::afw::table::SourceMatchVector const &matPair,
		     bool verbose = false);
	}
    }
}

#endif
