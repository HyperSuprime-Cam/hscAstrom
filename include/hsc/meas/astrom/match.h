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

    // Reference objects are now SimpleRecords, not Sources, and they don't have
    // x,y fields anymore (which makes sense, because they don't belong to a
    // particular frame).  But that complicates the implementation of this
    // algorithm a lot, because it operates completely in x,y space.
    
    // To deal with this, we have a RecordProxy object that we use for both
    // reference objects and sources internally; it holds a pointer to the
    // actual record and a point that it either copies from the source or
    // computes for the reference object from the initial WCS.

    struct RecordProxy {
        explicit RecordProxy() {}                // needed so we can call RecordProxy::resize()
        PTR(lsst::afw::table::SimpleRecord) record;
        lsst::afw::geom::Point2D position;

        double getX() const { return position.getX(); }
        double getY() const { return position.getY(); }

        operator PTR(lsst::afw::table::SimpleRecord) () const { return record; }

        bool operator==(RecordProxy const & other) const { return record == other.record; }
        bool operator!=(RecordProxy const & other) const { return record != other.record; }

        RecordProxy(
            PTR(lsst::afw::table::SimpleRecord) record_,
            lsst::afw::geom::Point2D const & position_
        ) : record(record_), position(position_) {}
    };

    typedef std::vector<RecordProxy> ProxyVector;

    ProxyVector makeProxies(lsst::afw::table::SourceCatalog const & sources);

    ProxyVector makeProxies(lsst::afw::table::SimpleCatalog const & refs, lsst::afw::image::Wcs const & wcs);

    struct ProxyPair {
        RecordProxy first;
        RecordProxy second;
        double distance;
        double pa;
        double deltaD;

	ProxyPair(RecordProxy const & s1, RecordProxy const & s2) : first(s1), second(s2) {
            double x1 = first.position.getX();
            double y1 = first.position.getY();
            double x2 = second.position.getX();
            double y2 = second.position.getY();
            distance = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
            pa = atan2(y2-y1, x2-x1);
        }
    };
    
    // swapped order of cat/src, because it's always cat,src everywhere else (including Match objects)
    lsst::afw::table::ReferenceMatchVector
    match(
        lsst::afw::table::SimpleCatalog const &cat,
        lsst::afw::table::SourceCatalog const &src,
        lsst::afw::image::Wcs const & wcs,
        int numBrightStars = 100,
        unsigned int minNumMatchedPair = 50,
        double matchingAllowanceInPixel = 10.,
        double offsetAllowedInPixel = 300.,
        bool verbose = false
    );

    }
  }
}

#endif
