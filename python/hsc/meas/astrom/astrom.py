#!/usr/bin/env python

import lsst.daf.base as dafBase
import lsst.meas.astrom.astrom as measAst
from lsst.pex.config import Config, Field, RangeField

class TaburAstrometryConfig(Config):
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
    calculateSip = Field(
        doc='''Compute polynomial SIP distortion terms?''',
        dtype=bool,
        default=True)
    sipOrder = RangeField(
        doc='''Polynomial order of SIP distortion terms''',
        dtype=int,
        default=4, min=2)


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

        sources = [s for s in sources if self._goodStar(s)]
        if self.log: self.log.log(self.log.INFO, "Matching to %d good input sources" % len(sources))

        minNumMatchedPair = min(self.config.minMatchedPairNumber,
                                int(self.config.minMatchedPairFrac * len(cat)))

        matchList = hscAstrom.match(sources, cat, self.config.numBrightStars, minNumMatchedPair)
        if matchList is None or len(matchList) == 0:
            raise RuntimeErorr("Unable to match sources")
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

