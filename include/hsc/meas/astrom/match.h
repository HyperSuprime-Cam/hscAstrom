// -*- lsst-c++ -*-
#if !defined(HSC_MATCH_OPM_B_MATCH_H)
#define HSC_MATCH_OPM_B_MATCH_H

#include <vector>
#include <cmath>
#include "lsst/afw/table/Source.h"
#include "lsst/afw/table/Match.h"

namespace hsc {
  namespace meas {
    namespace astrom {

	struct SourcePair {
	    CONST_PTR(lsst::afw::table::SourceRecord) first;
	    CONST_PTR(lsst::afw::table::SourceRecord) second;
	    double distance;
	    double pa;
	    double deltaD;

	    SourcePair(CONST_PTR(lsst::afw::table::SourceRecord) s1,
		       CONST_PTR(lsst::afw::table::SourceRecord) s2)
		: first(s1), second(s2) {
		double x1 = first->getX();
		double y1 = first->getY();
		double x2 = second->getX();
		double y2 = second->getY();
		distance = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		pa = atan2(y2-y1, x2-x1);
	    }
	    ~SourcePair() {}
	};
	    

	lsst::afw::table::SourceCatalog
	  selectPoint(lsst::afw::table::SourceCatalog const &a,
		      unsigned int n,
		      unsigned int start = 0);

	lsst::afw::table::SourceMatchVector
	  match(lsst::afw::table::SourceCatalog const &src,
		lsst::afw::table::SourceCatalog const &cat,
		int numBrightStars = 100,
		unsigned int minNumMatchedPair = 50,
		double matchingAllowanceInPixel = 10.,
		double offsetAllowedInPixel = 300.,
		bool verbose = false);
    }
  }
}

#endif
