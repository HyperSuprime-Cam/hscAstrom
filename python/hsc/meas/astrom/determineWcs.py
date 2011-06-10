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

def getImageSizeInDegree(srcSet, wcs):
    """ Get the approximate size of the image in arcseconds
    
    Input: 
    srcSet List of detected objects in the image (with pixel positions)
    wcs    Wcs converts pixel positions to ra/dec
    
    """
    xfunc = lambda x: x.getXAstrom()
    yfunc = lambda x: x.getYAstrom()
    
    x = map(xfunc, [s for s in srcSet if goodStar(s)])
    y = map(yfunc, [s for s in srcSet if goodStar(s)])

    minx = min(x)
    maxx = max(x)
    miny = min(y)
    maxy = max(y)
    
    llc = wcs.pixelToSky(minx, miny).getPosition()
    urc = wcs.pixelToSky(maxx, maxy).getPosition()
    
    deltaRa  = llc.getX() - urc.getX()
    deltaDec = urc.getY() - llc.getY()
    
    return deltaRa, deltaDec

def getCatalogueForField(solver, srcSet, wcsIn, imageSize, filterName, idName, margin):
    W, H = imageSize
    skyCenter = wcsIn.pixelToSky(W/2, H/2).getPosition()
    ra = skyCenter.getX()
    dec = skyCenter.getY()

    deltaRa, deltaDec = getImageSizeInDegree(srcSet, wcsIn)

    radius = max(deltaRa, deltaDec) * 3600.

    ids = solver.getIndexIdList()
    indexid_new = []
    ref_new = []
    idx_new = []
    for indexid in ids:
        if indexid % 10 != 0:
            continue;
        X = solver.getCatalogue(ra, dec, radius, '', idName, indexid)
        if len(X.refsources) != 0:
            indexid_new.append(indexid)
            cols = solver.getTagAlongColumns(indexid)
            colnames = [c.name for c in cols]

            tagdata = []
            for c in cols:
                fname = 'getTagAlong' + c.ctype
                func = getattr(solver, fname)
                data = func(indexid, c.name, X.inds)
                tagdata.append(data)

            i_filter = colnames.index(filterName)

            for i, (ref, idx) in enumerate(zip(X.refsources, X.inds)):
                ra = ref.getRa()/math.pi*180.
                dec = ref.getDec()/math.pi*180.
                pix = wcsIn.skyToPixel(ra, dec)
                if (-margin < pix[0] < W+margin and \
                    -margin < pix[1] < H+margin):
                    mag = tagdata[i_filter][i]
                    ref.setPsfFlux(math.pow(10.0, -0.4*mag))
                    ref_new.append(ref)
                    idx_new.append(idx)

    X.indexid = indexid_new[0]
    X.refsources = ref_new
    X.inds = idx_new
    return X

def queryReferenceCatalog(solver, srcSet, wcsIn, imageSize, filterName,
                          idName, margin=50):
    X = getCatalogueForField(solver, srcSet, wcsIn, imageSize, filterName, idName, margin)
    ref = X.refsources
    catSet = afwDet.SourceSet()
    if (len(ref) != 0):
        for i,r in enumerate(ref):
            ra = r.getRa()
            dec = r.getDec()
            p = wcsIn.skyToPixel(ra/math.pi*180., dec/math.pi*180.)
            mag = r.getPsfFlux()
            if (mag == mag):
                s = afwDet.Source()
                s.setId(r.getId())
                s.setXAstrom(p.getX())
                s.setYAstrom(p.getY())
                s.setRa(ra)
                s.setDec(dec)
                catSet.append(s)

    return catSet

def runMatch(solver, wcsIn, srcSet, numBrightStars, imageSize, filterName, idName):
    
    catSet = queryReferenceCatalog(solver, srcSet, wcsIn, imageSize, filterName, idName)

    srcSet2 = [s for s in srcSet if goodStar(s)]
    matchList = hscAstrom.match(srcSet2, catSet, numBrightStars)

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

    #Extract an initial guess WCS if available    
    wcsIn = exposure.getWcs() #May be None
    # Exposure uses the special object "NoWcs" instead of NULL.  Because they're special.
    haswcs = exposure.hasWcs()
    if not haswcs:
        log.log(log.WARN, "No WCS found in exposure. Doing blind solve")

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

    key = 'pixelScaleUncertainty'
    if policy.exists(key):
        dscale = float(policy.get(key))
    else:
        dscale = None

    # Do a blind solve if we're told to, or if we don't have an input WCS
    doBlindSolve = policy.get('blindSolve') or (not haswcs)
    if doBlindSolve:
        log.log(log.DEBUG, "Solving with no initial guess at position")
        isSolved = solver.solve()
    elif dscale is not None:
        isSolved = solver.solve(wcsIn, dscale)
    else:
#        isSolved = solver.solve(wcsIn)
        isSolved, wcs, matchList = runMatch(solver, wcsIn, sourceSet,
                            min(policy.get('numBrightStars'), len(sourceSet)),
                            (W,H), filterName, measAst.getIdColumn(policy))
        log.log(log.INFO, "Found %d matches in hscAstrom" % 0 if matchList is None else len(matchList))

    # Did we solve?
    log.log(log.DEBUG, 'Finished astrometric solution')
    if not isSolved:
        log.log(log.WARN, "No astrometric solution found, using input WCS")
        return astrom
#    wcs = solver.getWcs()

    # Generate a list of catalogue objects in the field.
    imgSizeInArcsec = wcs.pixelScale() * math.hypot(W,H)
    filterName = measAst.chooseFilterName(exposure, policy, solver, log, filterName)
    idName = measAst.getIdColumn(policy)
    try:
        margin = 50 # pixels
        #X = solver.getCatalogueForSolvedField(filterName, idName, margin)
        X = getCatalogueForField(solver, sourceSet, wcs, (W,H), filterName, idName, margin)
        cat = X.refsources
        indexid = X.indexid
        inds = X.inds
    except LsstCppException, e:
        log.log(Log.WARN, str(e))
        log.log(Log.WARN, "Attempting to access catalogue positions and fluxes")
        version = os.environ['ASTROMETRY_NET_DATA_DIR']
        log.log(Log.WARN, "Catalogue version: %s" %(version))
        log.log(Log.WARN, "ID column: %s" %(idName))
        log.log(Log.WARN, "Requested filter: %s" %(filterName))
        log.log(Log.WARN, "Available filters: " + str(solver.getCatalogueMetadataFields()))
        raise

    measAst.addTagAlongValuesToReferenceSources(solver, policy, log, cat, indexid, inds, filterName)
    
    if True:
        # Now generate a list of matching objects
        distInArcsec = policy.get('distanceForCatalogueMatchinArcsec')
        cleanParam = policy.get('cleaningParameter')

        matchList = measAst.matchSrcAndCatalogue(cat=cat, img=sourceSet, wcs=wcs,
            distInArcsec=distInArcsec, cleanParam=cleanParam)

        uniq = set([sm.second.getId() for sm in matchList])
        if len(matchList) != len(uniq):
            log.log(Log.WARN, "The list of matches stars contains duplicated reference sources (%i sources, %i unique ids)"
                    % (len(matchList), len(uniq)))

        if len(matchList) == 0:
            log.log(Log.WARN, "No matches found between input source and catalogue.")
            log.log(Log.WARN, "Something is wrong. Defaulting to input WCS")
            return astrom

        log.log(Log.DEBUG, "%i catalogue objects match input source list using linear WCS" %(len(matchList)))
    else:
        # Use list of matches returned by Astrometry.net
        log.log(Log.DEBUG, "Getting matched sources: Fluxes in column %s; Ids in column" % (filterName, idName))
        matchList = solver.getMatchedSources(filterName, idName)

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
        wcs, matchList = measAst.calculateSipTerms(wcs, cat, sourceSet, distInArcsec, cleanParam, sipOrder, log)

        astrom.sipWcs = wcs
        astrom.sipMatches = matchList
    else:
        log.log(Log.DEBUG, "Updating WCS in input exposure with linear WCS")

    log.log(Log.DEBUG, "Setting exposure's WCS: to\n" + wcs.getFitsMetadata().toString())
    exposure.setWcs(wcs)

    # add current EUPS astrometry_net_data setup.
    moreMeta = dafBase.PropertyList()
    andata = os.environ.get('ASTROMETRY_NET_DATA_DIR')
    if andata is None:
        moreMeta.add('ANEUPS', 'none', 'ASTROMETRY_NET_DATA_DIR')
    else:
        andata = os.path.basename(andata)
        moreMeta.add('ANEUPS', andata, 'ASTROMETRY_NET_DATA_DIR')

    # cache: field center and size.  These may be off by 1/2 or 1 or 3/2 pixels.
    # dstn does not care.
    cx,cy = W/2.,H/2.
    radec = wcs.pixelToSky(cx, cy)
    ra,dec = radec.getLongitude(afwCoord.DEGREES), radec.getLatitude(afwCoord.DEGREES)
    moreMeta.add('RA', ra, 'field center in degrees')
    moreMeta.add('DEC', dec, 'field center in degrees')
    moreMeta.add('RADIUS', imgSizeInArcsec/2./3600.,
            'field radius in degrees, approximate')
    moreMeta.add('SMATCHV', 1, 'SourceMatchVector version number')

    if display:
        for s1, s2, d in matchList:
            # plot the catalogue positions
            ds9.dot("+", s1.getXAstrom(), s1.getYAstrom(), size=3, ctype=ds9.BLUE, frame=frame)

    #matchListMeta = solver.getMatchedIndexMetadata()
    matchListMeta = dafBase.PropertyList()
    matchListMeta.add('ANINDID', X.indexid, 'Astrometry.net index id')
    moreMeta.combine(matchListMeta)

    astrom.matchMetadata = moreMeta
    astrom.wcs = wcs
    astrom.matches = matchList

    return astrom
