import os
import math
from lsst.pex.logging import Log
from lsst.pex.exceptions import LsstCppException
import lsst.meas.astrom        as measAst
import lsst.meas.algorithms    as measAlg
import lsst.afw.detection      as afwDet
import lsst.daf.base           as dafBase
import lsst.afw.coord          as afwCoord
import lsst.afw.display.ds9    as ds9
import astromLib as hscAstrom

try:
    import lsstDebug

    display = lsstDebug.Info(__name__).display
except ImportError, e:
    try:
        type(display)
    except NameError:
        display = False

def goodStar(s):
    return s.getXAstrom() == s.getXAstrom() and \
           not (s.getFlagForDetection() & measAlg.Flags.SATUR)

def runMatch(solver, wcsIn, srcSet, numBrightStars, meta, policy, filterName, log=None):
    catSet = measAst.readReferenceSourcesFromMetadata(meta, log=log, policy=policy, filterName=filterName)
    if log is not None: log.log(log.INFO, "Retrieved %d reference catalog sources" % len(catSet))

    srcSet2 = [s for s in srcSet if goodStar(s)]
    if log is not None: log.log(log.INFO, "Matching to %d good input sources" % len(srcSet2))

    minNumMatchedPair = 30
    if minNumMatchedPair < len(catSet)/5:
        minNumMatchedPair = len(catSet)/5
    matchList = hscAstrom.match(srcSet2, catSet,
                                numBrightStars,
                                minNumMatchedPair)

    if len(matchList) != 0:
        wcsOut = hscAstrom.fitTAN(matchList)

    if len(matchList) != 0:
        return True, wcsOut, matchList
    else:
        return False, None, None

def determineWcs(policy, exposure, sourceSet, log=None, solver=None, doTrim=False,
                 forceImageSize=None, filterName=None):
    '''Top level function for calculating an initial (per-chip) astrometric solution.

    Get an initial World Coordinate System (WCS) from Astrometry.net,
    then calculate SIP distortion terms.

    Input:
    policy:     An lsst.pex.policy.Policy object containing the parameters for the solver
    exposure    lsst.afw.image.Exposure representation of an image and a WCS 
                this provides the initial guess at position and plate scale
    sourceSet   A list of lsst.afw.detection.Source objects, indicating the pixel positions of
                stars in the field
    log         A lsst.pex.logging.Log object (optional), used for printing progress
    doTrim      Remove sources that are not inside the image.
    solver      Optionally provide a previously created astrometry.net solver. If not provided
                one will be created.
    forceImageSize  tuple of (W,H): force this image size, rather than getting it from the Exposure.
    filterName  Use this filter name, rather than getting it from the exposure.
    '''

    astrom = measAst.InitialAstrometry()

    if log is None:
        log = Log.getDefaultLog()

    if display:
        frame = 1
        ds9.mtv(exposure, frame=frame, title="wcsDet")

    if doTrim:
        nStart = len(sourceSet)
        sourceSet = measAst.trimBadPoints(exposure, sourceSet)
        if log:
            nEnd = len(sourceSet)
            log.log(log.DEBUG, "Kept %i of %i sources after trimming" %(nEnd, nStart))

    if display:
        for s in sourceSet:
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=3, ctype=ds9.RED, frame=frame)

    if not exposure.hasWcs():
        log.log(log.WARN, "No WCS found in exposure; hsc.meas.astrom fails")
        raise RuntimeError("No WCS found in expousre")
    wcsIn = exposure.getWcs()
    
    # Setup solver
    if solver is None:
        solver = measAst.createSolver(policy, log)
    else:
        solver.reset()

    # Set solving params
    log.log(log.DEBUG, "Setting starlist")
    solver.setStarlist(sourceSet)
    log.log(log.DEBUG, "Setting numBrightObj")
    solver.setNumBrightObjects( min(policy.get('numBrightStars'), len(sourceSet)))
    if forceImageSize is not None:
        (W,H) = forceImageSize
    else:
        (W,H) = (exposure.getWidth(), exposure.getHeight())
    solver.setImageSize(W, H)
    #solver.printSolverSettings(stdout)

    filterName = measAst.chooseFilterName(exposure, policy, solver, log, filterName)
    stargalName, variableName, magerrName = measAst.getTagAlongNamesFromPolicy(policy, filterName)
    meta = measAst.createMetadata(W, H, wcsIn, filterName, stargalName, variableName, magerrName)

    isSolved, wcs, matchList = runMatch(solver, wcsIn, sourceSet,
                                        min(policy.get('numBrightStars'), len(sourceSet)),
                                        meta, policy, filterName, log=log)
    if isSolved:
        log.log(log.INFO, "Found %d matches in hscAstrom" % (0 if matchList is None else len(matchList)))

    # Did we solve?
    log.log(log.DEBUG, 'Finished astrometric solution')
    if not isSolved:
        log.log(log.WARN, "No astrometric solution found")
        return None

    astrom.tanWcs = wcs
    astrom.tanMatches = matchList

    srcids = [s.getSourceId() for s in sourceSet]
    #print 'srcids:', srcids
    for m in matchList:
        #print 'Matchlist entry ids:', m.first.getSourceId(), m.second.getSourceId()
        assert(m.second.getSourceId() in srcids)
        assert(m.second in sourceSet)

    if policy.get('calculateSip'):
        sipOrder = policy.get('sipOrder')
        wcs, matchList = measAst.calculateSipTerms(wcs, cat, sourceSet, distInArcsec, 
                                                   cleanParam, sipOrder, log)

        astrom.sipWcs = wcs
        astrom.sipMatches = matchList
    else:
        log.log(Log.DEBUG, "Updating WCS in input exposure with linear WCS")

    log.log(Log.DEBUG, "Setting exposure's WCS: to\n" + wcs.getFitsMetadata().toString())
    exposure.setWcs(wcs)

    if display:
        for s1, s2, d in matchList:
            # plot the catalogue positions
            ds9.dot("+", s1.getXAstrom(), s1.getYAstrom(), size=3, ctype=ds9.BLUE, frame=frame)

    astrom.matchMetadata = meta
    astrom.wcs = wcs
    astrom.matches = matchList

    return astrom
