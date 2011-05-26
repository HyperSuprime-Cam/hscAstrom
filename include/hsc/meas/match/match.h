// -*- lsst-c++ -*-
#if !defined(HSC_MATCH_OPM_B_MATCH_H)
#define HSC_MATCH_OPM_B_MATCH_H

#include <vector>
#include <cmath>
#include "lsst/afw/detection/Source.h"
#include "lsst/afw/detection/SourceMatch.h"

namespace hsc {
  namespace meas {
    namespace match {

	struct SourcePair {
	    lsst::afw::detection::Source::Ptr first;
	    lsst::afw::detection::Source::Ptr second;
	    double distance;
	    double pa;
	    double deltaD;

	    SourcePair(lsst::afw::detection::Source::Ptr const & s1,
		       lsst::afw::detection::Source::Ptr const & s2)
		: first(s1), second(s2) {
		double x1 = first->getXAstrom();
		double y1 = first->getYAstrom();
		double x2 = second->getXAstrom();
		double y2 = second->getYAstrom();
		distance = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		pa = atan2(y2-y1, x2-x1);
	    }
	    ~SourcePair() {}
	};
	    

	lsst::afw::detection::SourceSet selectPoint(lsst::afw::detection::SourceSet const &a,
						    int n,
						    int start = 0);

	std::vector<lsst::afw::detection::SourceMatch> match(lsst::afw::detection::SourceSet const &src,
							     lsst::afw::detection::SourceSet const &cat,
							     int numBrightStars = 100,
							     bool verbose = false);
    }
  }
}

#endif
