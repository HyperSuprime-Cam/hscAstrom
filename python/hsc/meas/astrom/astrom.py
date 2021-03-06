#!/usr/bin/env python

import math, numpy, os
import lsst.daf.base as dafBase
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.meas.astrom as measAst
import lsst.pex.config as pexConfig
from . import astromLib as hscAstrom
import lsst.afw.cameraGeom as cameraGeom

class TaburAstrometryConfig(measAst.MeasAstromConfig):
    numBrightStars = pexConfig.RangeField(
        doc="Number of bright stars to use",
        dtype=int,
        default=50, min=2)
    minMatchedPairNumber = pexConfig.RangeField(
        doc="Minimum number of matched pairs",
        dtype=int,
        default=30, min=2)
    minMatchedPairFrac = pexConfig.RangeField(
        doc="Minimum number of matched pairs, expressed as a fraction of the reference catalogue size",
        dtype=float,
        default=0.3, min=0, max=1)
    pixelMargin = pexConfig.RangeField(
        doc="Padding to add to image size (pixels)",
        dtype=int,
        default=50, min=0)
    calculateSip = pexConfig.Field(
        doc='''Compute polynomial SIP distortion terms?''',
        dtype=bool,
        default=True)
    sipOrder = pexConfig.RangeField(
        doc='''Polynomial order of SIP distortion terms''',
        dtype=int,
        default=4, min=1)
    offsetAllowedInPixel = pexConfig.RangeField(
        doc="Offset between sources and catalog allowed (pixel)",
        dtype=int,
        default=300, max=4000)
    rotationAllowedInRad = pexConfig.RangeField(
        doc="Roation angle allowed between sources and catalog (radian)",
        dtype=float,
        default=0.02, max=0.1)
    angleDiffFrom90 = pexConfig.RangeField(
        doc="Difference of angle between x and y from 90 degree allowed (degree)",
        dtype=float,
        default=0.2, max=45.0)
    numPointsForShape = pexConfig.Field(
        doc="number of points to define a shape for matching",
        dtype = int,
        default = 6)
    limitOnDeterminant = pexConfig.Field(
        doc="limit on determinant of linear transforming matrix",
        dtype = float,
        default = 0.02)
    # Set these for proper operation of overridden astrometry class
    useWcsPixelScale = True
    useWcsRaDecCenter = True
    useWcsParity = True

def cleanStar(s, prefix):
    return (numpy.isfinite(s.getX()) and numpy.isfinite(s.getY()) and not s.get(prefix+"flags.pixel.edge"))

def goodStar(s, prefix):
    # FIXME: should use Key to get flag (but then we'd need schema in advance)
    return (cleanStar(s, prefix) and not s.getCentroidFlag() and not s.get(prefix+"flags.pixel.saturated.any"))

def show(debug, exposure, wcs, sources, catalog, matches=[], correctDistortion=True, frame=1, title=""):
    import lsst.afw.display.ds9 as ds9
    import numpy

    if not debug.display:
        return

    ds9.mtv(exposure, frame=frame, title=title)

    distorter = None
    if correctDistortion:
        try:
            detector = exposure.getDetector()
            distorter = detector.getDistortion()
            def toObserved(x, y):
                dist = distorter.distort(afwGeom.Point2D(x, y), detector)
                return dist.getX(), dist.getY()
        except Exception, e:
            print "WARNING: Unable to use distortion: %s" % e
            distorter = None
    if distorter is None:
        toObserved = lambda x,y: (x,y)

    with ds9.Buffering():
        if matches:
            for s in sources:
                x, y = toObserved(s.getX(), s.getY())
                ds9.dot("+", x,  y, size=10, frame=frame, ctype=ds9.GREEN)

            for s in catalog:
                x, y = wcs.skyToPixel(s.getCoord())
                x, y = toObserved(x, y)
                ds9.dot("x", x, y, size=10, frame=frame, ctype=ds9.RED)

            dr = numpy.ndarray(len(matches))

            for i, m in enumerate(matches):
                x, y = m.second.getX(), m.second.getY()
                pix = wcs.skyToPixel(m.first.getCoord())

                dr[i] = numpy.hypot(pix[0] - x, pix[1] - y)

                x, y = toObserved(x, y)
                ds9.dot("o", x,  y, size=10, frame=frame, ctype=ds9.YELLOW)
                
            print "<dr> = %.4g +- %.4g [%d matches]" % (dr.mean(), dr.std(), len(matches))
        else:
            for s in sources:
                x0, y0 = s.getX(), s.getY()
                ds9.dot("+", x0, y0, size=20, frame=frame, ctype=ds9.GREEN)
                if correctDistortion:
                    x, y = toObserved(x0, y0)
                    ds9.dot("o", x,  y,  frame=frame, ctype=ds9.GREEN)
                    ds9.line([(x0, y0), (x, y)], frame=frame, ctype=ds9.GREEN)

            for s in catalog:
                pix = wcs.skyToPixel(s.getCoord())
                ds9.dot("x", pix[0], pix[1], size=20, frame=frame, ctype=ds9.RED)

def getDistorter(exposure):
    try:
        detector = exposure.getDetector()
        distorter = detector.getDistortion()
        def toUndistort(x, y):
            dist = distorter.undistort(afwGeom.Point2D(x, y), detector)
            return dist.getX(), dist.getY()
        def toDistort(x, y):
            dist = distorter.distort(afwGeom.Point2D(x, y), detector)
            return dist.getX(), dist.getY()
    except Exception, e:
        print "WARNING: Unable to use distortion: %s" % e
        distorter = None
    if distorter is None:
        passThrough = lambda x,y: (x,y)
        toUndistort = passThrough
        toDistort = passThrough
    return toUndistort, toDistort

def getUndistortedXY0(exposure):
    toUndistort, _ = getDistorter(exposure)

    x0, y0 = toUndistort(0,0)
    x1 = x0
    y1 = y0
    for x, y in [(exposure.getWidth(), 0),
                 (exposure.getWidth(),exposure.getHeight()),
                 (0,exposure.getHeight())]:
        xx, yy = toUndistort(x, y)
        if xx < x0:
            x0 = xx
        if xx > x1:
            x1 = xx
        if yy < y0:
            y0 = yy
        if yy > y1:
            y1 = yy
    
    return x0, y0, x1-x0, y1-y0

def isAmpDeadWithSources(amp, sources, toDistort, prefix=''):
    """ Check whether specified amp is dead or not
        by checking the number of sources in that amp.
    """
    centroidKey = sources.getCentroidKey()
    dataBox = afwGeom.Box2D(amp.getDataSec(True))
    for s in sources:
        if cleanStar(s, prefix):
            center = s.get(centroidKey)
            x, y = toDistort(center.getX(), center.getY())
            if dataBox.contains(afwGeom.Point2D(x,y)):
                return False
    return True

class TaburAstrometry(measAst.Astrometry):
    """Star matching using algorithm based on V.Tabur 2007, PASA, 24, 189-198
    ("Fast Algorithms for Matching CCD Images to a Stellar Catalogue")
    """
    ConfigClass = TaburAstrometryConfig

    def determineWcs(self, sources, exposure):
        import lsstDebug
        debug = lsstDebug.Info(__name__)
        debugging = lsstDebug.Info(__name__).display

        # Get measurement prefix
        psfFluxDef = sources.getTable().getPsfFluxDefinition()
        prefix = psfFluxDef[:psfFluxDef.find('flux.psf')]

        # Distorter to de-apply distortion correction to get sources' centroid
        # in original pixel coordinate.
        _, toDistort = getDistorter(exposure)

        wcs = exposure.getWcs() # Guess WCS
        if wcs is None:
            raise RuntimeError("This matching algorithm requires an input guess WCS")

        x0, y0, w, h = getUndistortedXY0(exposure)

        if wcs.hasDistortion():
            # Need to back out the guestimated distortion already applied in the source positions
            # because we make use of the WCS which has distortion in it.
            _, toDistort = getDistorter(exposure)
            centroidKey = sources.getCentroidKey()
            for s in sources:
                center = s.get(centroidKey)
                x,y = toDistort(center.getX(), center.getY())
                s.set(centroidKey, afwGeom.Point2D(x, y))

        filterName = exposure.getFilter().getName()
        imageSize = (exposure.getWidth(), exposure.getHeight())

        # For some chips, a few amplifiers are dead.
        # Reference catalog sources which should be observed on these dead amps
        # can be contamination when finding astrometric matches.
        # The following code will find out live amps by checking number of sources
        # in each amps and determine effective area in pixel coordinate.
        bounds = None
        ccd = cameraGeom.cast_Ccd(exposure.getDetector())
        for amp in ccd:
            if not isAmpDeadWithSources(amp, sources, toDistort, prefix):
                if bounds is None:
                    bounds = amp.getDataSec(True)
                else:
                    bounds.include(amp.getDataSec(True))

        # Find reference catalog only for the effective area.
        wcs.shiftReferencePixel(-bounds.getMinX(), -bounds.getMinY())
        cat = self.getReferenceSourcesForWcs(wcs, bounds.getDimensions(), filterName, self.config.pixelMargin, x0, y0)
        wcs.shiftReferencePixel(bounds.getMinX(), bounds.getMinY())

        # Select unique objects only
        keep = type(cat)(cat.table)
        for c in cat:
            found = False
            for k in keep:
                if c.getCoord() == k.getCoord():
                    found = True
                    break
            if not found:
                keep.append(c)
        cat = keep

        # Avoid M31/M32 center
        keep = type(cat)(cat.table)
        import lsst.afw.coord as afwCoord
        m31 = afwCoord.Coord(afwGeom.Point2D(10.681, 41.269))
        m32 = afwCoord.Coord(afwGeom.Point2D(10.675, 40.864))
        for c in cat:
            if (c.getCoord().angularSeparation(m31).asArcminutes() > 2.2 and
                c.getCoord().angularSeparation(m32).asArcminutes() > 1.0):
                keep.append(c)
        cat = keep

        if self.log: self.log.info("Found %d catalog sources" % len(cat))

        #allSources = sources
        allSources = afwTable.SourceCatalog(sources.table)
        allSources.extend(s for s in sources if cleanStar(s, prefix))
        sources = afwTable.SourceCatalog(sources.table)
        sources.extend(s for s in allSources if goodStar(s, prefix))
        
        if self.log: self.log.info("Matching to %d/%d good input sources" % 
                                   (len(sources), len(allSources)))

        numBrightStars = max(self.config.numBrightStars, len(sources))
        minNumMatchedPair = min(self.config.minMatchedPairNumber,
                                int(self.config.minMatchedPairFrac * min([len(cat), len(sources)])))

        show(debug, exposure, wcs, sources, cat, frame=debug.frame1 if isinstance(debug.frame1, int) else 1,
             title="Input catalog")

        matchingRadius = self.config.catalogMatchDist / wcs.pixelScale().asArcseconds() # in pixels

        def doMatching(sources):
            """Attempt matching sources"""
            for i in range(4):
                e_dpa = self.config.rotationAllowedInRad * math.pow(2.0, 0.5*i)
                for j in range(3):
                    matchingRadius = (self.config.catalogMatchDist/wcs.pixelScale().asArcseconds()*
                                      math.pow(1.25, j))
                    for k in range(3):
                        matchList = hscAstrom.match(cat, sources, wcs, self.config.numBrightStars,
                                                    minNumMatchedPair, matchingRadius,
                                                    len(allSources)-len(sources),
                                                    self.config.offsetAllowedInPixel,
                                                    e_dpa, self.config.angleDiffFrom90*(k+1),
                                                    self.config.numPointsForShape,
                                                    self.config.limitOnDeterminant,
                                                    debug.verbose)
                        if matchList is not None and len(matchList) != 0:
                            return matchList
                    if matchList is not None and len(matchList) != 0:
                        return matchList
                if matchList is not None and len(matchList) != 0:
                    return matchList
        
        # match sources with staturated, then select unsaturated only 
        matchList = doMatching(allSources)
        if matchList is not None:
            matchList_0 = [m for m in matchList 
                           if (not m.second.getCentroidFlag() and 
                               not m.second.get(prefix+"flags.pixel.saturated.any"))]
        else:
            matchList_0 = []

        # match sources without saturated
        if len(cat) > len(allSources) - len(sources):
            matchList_1 = doMatching(sources)
        else:
            matchList_1 = []
        if matchList_1 == None:
            matchList_1 = []

        if len(matchList_0) == 0 and len(matchList_1) == 0:
            raise RuntimeError("Unable to match sources")

        # Adopt matchList with more matches
        if len(matchList_0) > len(matchList_1):
            matchList = matchList_0
        else:
            matchList = matchList_1

        if self.log: self.log.info("Matched %d sources" % len(matchList))
        if len(matchList) < minNumMatchedPair:
            self.log.warn("Number of matches is smaller than request")


        wcs = hscAstrom.fitTAN(matchList, True if debugging else False)

        if debug.showLinear:
            show(debug, exposure, wcs, sources, cat, matches=matchList,
                 frame=debug.frame2 if isinstance(debug.frame2, int) else 2, title="Linear matches")

        astrom = measAst.astrom.InitialAstrometry()
        astrom.tanMatches = matchList
        astrom.tanWcs = wcs
        astrom.solveQa = dafBase.PropertyList()

        if self.config.calculateSip:
            wcs, matchList = self._calculateSipTerms(wcs, cat, sources, matchList, imageSize)
            astrom.sipWcs = wcs
            astrom.sipMatches = matchList

        exposure.setWcs(wcs)

        astrom.matchMetadata = measAst.astrom._createMetadata(imageSize[0], imageSize[1], x0, y0,
                                                              wcs, filterName)
        astrom.sipWcs = wcs
        astrom.sipMatches = afwTable.ReferenceMatchVector()
        for m in matchList:
            astrom.sipMatches.push_back(m)

        if self.config.calculateSip:
            show(debug, exposure, wcs, sources, cat, matches=matchList,
                 frame=debug.frame3 if isinstance(debug.frame3, int) else 3, title="SIP matches")

        return astrom
