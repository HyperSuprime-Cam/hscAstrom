#!/usr/bin/env python

import numpy
import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.meas.algorithms as measAlg
import lsst.meas.astrom as measAst
from . import astromLib as hscAstrom
from lsst.pex.config import Config, Field, RangeField
from lsst.pipe.tasks.detectorUtil import getCcd

class TaburAstrometryConfig(measAst.MeasAstromConfig):
    numBrightStars = RangeField(
        doc="Number of bright stars to use",
        dtype=int,
        default=50, min=2)
    minMatchedPairNumber = RangeField(
        doc="Minimum number of matched pairs",
        dtype=int,
        default=30, min=2)
    minMatchedPairFrac = RangeField(
        doc="Minimum number of matched pairs, expressed as a fraction of the reference catalogue size",
        dtype=float,
        default=0.2, min=0, max=1)
    pixelMargin = RangeField(
        doc="Padding to add to image size (pixels)",
        dtype=int,
        default=50, min=0)
    matchingRadius = RangeField(
        doc="Radius within which to accept matches (pixels)",
        dtype=float,
        default=10.0, min=0.0)
    calculateSip = Field(
        doc='''Compute polynomial SIP distortion terms?''',
        dtype=bool,
        default=True)
    sipOrder = RangeField(
        doc='''Polynomial order of SIP distortion terms''',
        dtype=int,
        default=4, min=1)
    # Set these for proper operation of overridden astrometry class
    useWcsPixelScale = True
    useWcsRaDecCenter = True
    useWcsParity = True


def goodStar(s):
    # FIXME: should use Key to get flag (but then we'd need schema in advance)
    return numpy.isfinite(s.getX()) and numpy.isfinite(s.getY()) and not s.get("flags.pixel.saturated.any")

def show(debug, exposure, wcs, sources, catalog, matches=[], frame=1, title=""):
    import lsst.afw.display.ds9 as ds9
    import numpy

    if not debug.display:
        return

    ds9.mtv(exposure, frame=frame, title=title)

    ccd = getCcd(exposure, allowRaise=False)
    distorter = ccd.getDistortion()

    with ds9.Buffering():
        if matches:
            for s in sources:
                x, y = distorter.toObserved(s.getX(), s.getY())
                ds9.dot("+", x,  y,  frame=frame, ctype=ds9.GREEN)

            for s in catalog:
                x, y = wcs.skyToPixel(s.getCoord())
                x, y = distorter.toObserved(x, y)
                ds9.dot("x", x, y, size=3, frame=frame, ctype=ds9.RED)

            dr = numpy.ndarray(len(matches))

            for i, m in enumerate(matches):
                x, y = m.second.getX(), m.second.getY()
                pix = wcs.skyToPixel(m.first.getCoord())

                dr[i] = numpy.hypot(pix[0] - x, pix[1] - y)

                x, y = distorter.toObserved(x, y)
                ds9.dot("o", x,  y, size=4, frame=frame, ctype=ds9.YELLOW)
                
            print "<dr> = %.4g +- %.4g [%d matches]" % (dr.mean(), dr.std(), len(matches))
        else:
            for s in sources:
                x0, y0 = s.getX(), s.getY()
                x, y = distorter.toObserved(x0, y0)
                ds9.dot("+", x0, y0, size=3, frame=frame, ctype=ds9.GREEN)
                ds9.dot("o", x,  y,  frame=frame, ctype=ds9.GREEN)
                ds9.line([(x0, y0), (x, y)], frame=frame, ctype=ds9.GREEN)

            for s in catalog:
                pix = wcs.skyToPixel(s.getCoord())
                ds9.dot("x", pix[0], pix[1], size=3, frame=frame, ctype=ds9.RED)

class TaburAstrometry(measAst.Astrometry):
    """Star matching using algorithm based on V.Tabur 2007, PASA, 24, 189-198
    ("Fast Algorithms for Matching CCD Images to a Stellar Catalogue")
    """
    ConfigClass = TaburAstrometryConfig

    def determineWcs(self, sources, exposure):
        import lsstDebug
        debug = lsstDebug.Info(__name__)
        debugging = lsstDebug.Info(__name__).display
        
        wcs = exposure.getWcs() # Guess WCS
        if wcs is None:
            raise RuntimeError("This matching algorithm requires an input guess WCS")
        filterName = exposure.getFilter().getName()
        imageSize = (exposure.getWidth(), exposure.getHeight())
        cat = self.getReferenceSourcesForWcs(wcs, imageSize, filterName, self.config.pixelMargin)
        if self.log: self.log.log(self.log.INFO, "Found %d catalog sources" % len(cat))

        allSources = sources
        sources = afwTable.SourceCatalog(sources.table)
        sources.extend(s for s in allSources if goodStar(s))
        
        if self.log: self.log.log(self.log.INFO, "Matching to %d/%d good input sources" % 
                                  (len(sources), len(allSources)))

        numBrightStars = max(self.config.numBrightStars, len(sources))
        minNumMatchedPair = min(self.config.minMatchedPairNumber,
                                int(self.config.minMatchedPairFrac * len(cat)))

        show(debug, exposure, wcs, sources, cat,
             frame=debug.frame1 if isinstance(debug.frame1, int) else 1, title="Input catalog")

        try:
            matchList = hscAstrom.match(cat, sources, wcs, self.config.numBrightStars, minNumMatchedPair,
                                        self.config.matchingRadius)
            if matchList is None or len(matchList) == 0:
                raise RuntimeError("Unable to match sources")
        except:
            if self.log: self.log.info("Matching failed; retrying with saturated sources.")
            matchList = hscAstrom.match(cat, allSources, wcs, self.config.numBrightStars, minNumMatchedPair,
                                        self.config.matchingRadius)
            if matchList is None or len(matchList) == 0:
                raise RuntimeError("Unable to match sources")
        
        if self.log: self.log.log(self.log.INFO, "Matched %d sources" % len(matchList))


        wcs = hscAstrom.fitTAN(matchList, True if debugging else False)

        if debug.showLinear:
            show(debug, exposure, wcs, sources, cat, matches=matchList,
                 frame=debug.frame2 if isinstance(debug.frame2, int) else 2, title="Linear matches")

        astrom = measAst.astrom.InitialAstrometry()
        astrom.tanMatches = matchList
        astrom.tanWcs = wcs
        astrom.solveQa = dafBase.PropertyList()

        if self.config.calculateSip:
            wcs, matchList = self._calculateSipTerms(wcs, cat, sources, matchList)
            astrom.sipWcs = wcs
            astrom.sipMatches = matchList

        exposure.setWcs(wcs)

        astrom.matchMetadata = measAst.astrom._createMetadata(imageSize[0], imageSize[1], wcs, filterName)
        astrom.wcs = wcs
        astrom.matches = afwTable.ReferenceMatchVector()
        for m in matchList:
            astrom.matches.push_back(m)

        if self.config.calculateSip:
            show(debug, exposure, wcs, sources, cat, matches=matchList,
                 frame=debug.frame3 if isinstance(debug.frame3, int) else 3, title="SIP matches")

        return astrom

