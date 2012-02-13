#!/usr/bin/env python

import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
import lsst.meas.astrom.astrom as measAst
import hsc.meas.astrom.astromLib as hscAstrom
from lsst.pex.config import Config, Field, RangeField

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
    raDecSearchRadius = RangeField(
        '''When useWcsRaDecCenter=True, this is the radius, in degrees, around the RA,Dec center specified in the input exposure\'s WCS to search for a solution.''',
        float,
        default=1., min=0.)
    calculateSip = Field(
        doc='''Compute polynomial SIP distortion terms?''',
        dtype=bool,
        default=True)
    sipOrder = RangeField(
        doc='''Polynomial order of SIP distortion terms''',
        dtype=int,
        default=4, min=2)
    # Set these for proper operation of overridden astrometry class
    useWcsPixelScale = True
    useWcsRaDecCenter = True
    useWcsParity = True


def goodStar(s):
    return s.getXAstrom() == s.getXAstrom() and not (s.getFlagForDetection() & measAlg.Flags.SATUR)


class TaburAstrometry(measAst.Astrometry):
    """Star matching using algorithm based on V.Tabur 2007, PASA, 24, 189-198
    ("Fast Algorithms for Matching CCD Images to a Stellar Catalogue")
    """
    ConfigClass = TaburAstrometryConfig

    
    def determineWcs(self, sources, exposure):
        wcs = exposure.getWcs() # Guess WCS
        if wcs is None:
            raise RuntimeError("This matching algorithm requires an input guess WCS")
        filterName = exposure.getFilter().getName()
        imageSize = (exposure.getWidth(), exposure.getHeight())
        cat = self.getReferenceSourcesForWcs(wcs, imageSize, filterName, self.config.pixelMargin)
        if self.log: self.log.log(self.log.INFO, "Found %d catalog sources" % len(cat))

        sources = [s for s in sources if goodStar(s)]
        if self.log: self.log.log(self.log.INFO, "Matching to %d good input sources" % len(sources))

        numBrightStars = max(self.config.numBrightStars, len(sources))
        minNumMatchedPair = min(self.config.minMatchedPairNumber,
                                int(self.config.minMatchedPairFrac * len(cat)))

        matchList = hscAstrom.match(sources, cat, self.config.numBrightStars, minNumMatchedPair)
        if matchList is None or len(matchList) == 0:
            raise RuntimeError("Unable to match sources")
        if self.log: self.log.log(self.log.INFO, "Matched %d sources" % len(matchList))
        wcs = hscAstrom.fitTAN(matchList)

        astrom = measAst.InitialAstrometry()
        astrom.tanMatches = matchList
        astrom.tanWcs = wcs
        astrom.solveQa = dafBase.PropertyList()

        if self.config.calculateSip:
            wcs, matchList = self._calculateSipTerms(wcs, cat, sources, matchList)
            astrom.sipWcs = wcs
            astrom.sipMatches = matchList

        exposure.setWcs(wcs)
        astrom.matchMetadata = measAst._createMetadata(imageSize[0], imageSize[1], wcs, filterName)

        return astrom

