import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.suprimecam              as obsSc
import lsst.afw.image                   as afwImage
import lsst.afw.detection               as afwDet
import lsst.afw.display.ds9             as ds9
import hsc.meas.astrom                  as hscAst
import lsst.meas.algorithms             as measAlg
import lsst.pipette.util as pipUtil
import lsst.pipette.distortion as pipDist
import hsc.meas.match.matchLib as hscMatch
import math
import datetime
import lsst.pex.policy as policy
import lsst.meas.astrom as measAst
from lsst.pex.logging import Log

class SDSSstar:
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        self.values = {}

def queryReferenceCatalog(ra, dec, deltaRa, deltaDec):
# ra, dec, deltaRa and deltaDec in degree
    radius = max(deltaRa, deltaDec) * 3600.

    log = Log.getDefaultLog()
#    log.setThreshold(Log.DEBUG)
    pol = policy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAst.createSolver(pol, log)

#    solver.setLogLevel(3)

    ids = solver.getIndexIdList()
    catSet = afwDet.SourceSet()
    for indexid in ids:
        X = solver.getCatalogue(ra, dec, radius, '', 'id', indexid)
        ref = X.refsources
        inds = X.inds
        if (len(ref) != 0):
            cols = solver.getTagAlongColumns(indexid)
            colnames = [c.name for c in cols]

            tagdata = []
            for c in cols:
                fname = 'getTagAlong' + c.ctype
                func = getattr(solver, fname)
                data = func(indexid, c.name, inds)
                tagdata.append(data)

            stars = []
            for i,r in enumerate(ref):
                sra = r.getRa()/math.pi*180.
                sdec = r.getDec()/math.pi*180.
                if (sra > ra - 0.5 * deltaRa and sra < ra + 0.5 * deltaRa and
                    sdec > dec - 0.5 * deltaDec and sdec < dec + 0.5 * deltaDec):
                    star = SDSSstar(sra, sdec)
                    for c,d in zip(colnames, tagdata):
                        star.values[c] = d[i]
                    stars.append(star)

            print len(stars)
            return stars

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

def run(ioMgr, frameId, ccdId, display=False):

    data = {'visit': frameId, 'ccd':ccdId}
    exposure = ioMgr.read('calexp', data, ignore=True)[0]
    md = ioMgr.read('calexp_md', data, ignore=True)[0]
    wcsIn = afwImage.makeWcs(md)
    exposure.setWcs(wcsIn)

    if display:
        ds9.mtv(exposure, frame=1)

    butler = ioMgr.inButler
    srcSet = ioMgr.read('bsc', data, ignore=True)[0].getSources()

    if display:
        for s in srcSet:
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=10, frame=1)

    log = Log.getDefaultLog()
    pol = policy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAst.createSolver(pol, log)

    (W,H) = (exposure.getWidth(), exposure.getHeight())
    filterName = 'z'
    idName = 'id'

    catSet = hscAst.queryReferenceCatalog(solver, srcSet, wcsIn, (W,H),
                                          filterName, idName)

    if display:
        for s in catSet:
            ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=10, ctype=ds9.RED, frame=1)

    distConfig = {
        'radial': { 
            'coeffs': [-6.61572e-18, 5.69338e-14, 3.03146e-10, 7.16417e-08, 1.0, 0.0],
            'actualToIdeal': False,
            'step': 10.0
            }
        }
    ccd = pipUtil.getCcd(exposure)
    #print ccd.getId(), ccd.getCenter()
    dist = pipDist.createDistortion(ccd, distConfig)

    distSrc = dist.actualToIdeal(srcSet)
    #if display:
    #    for s in distSrc:
    #        ds9.dot("o", s.getXAstrom(), s.getYAstrom(), size=10, ctype=ds9.RED, frame=1)

    start = datetime.datetime.today()

    matchList = hscMatch.match(distSrc, catSet)

    end = datetime.datetime.today()

    matchList2 = []
    for m in matchList:
        x1, y1 = wcsIn.skyToPixel(m.first.getRa(), m.first.getDec())
        for s in srcSet:
            if (s.getId() == m.second.getId()):
                x0 = s.getXAstrom()
                y0 = s.getYAstrom()
                #print x0, y0, x1, y1, m.second.getId()
                if display:
                    ds9.line([[x0, y0], [x1, y1]], frame=1)
                m2 = afwDet.SourceMatch(m.first, s, 0.0)
                matchList2.append(m2)

    elapsedtime = end - start

    print "%6d %d %3d %6.2f" % (frameId, ccdId, len(matchList), elapsedtime.seconds + elapsedtime.microseconds/1000000.)
#    f = open('resultlog', 'a')
#    f.write("%6d %d %3d %6.2f\n" % (frameId, ccdId, len(matchList), elapsedtime.seconds + elapsedtime.microseconds/1000000.))
#    f.close()

if __name__ == '__main__':
    field = "ACTJ0022M0036"
    filter = "W-S-ZR"
    rerun = "yasuda-sup"

    mapper = obsSc.SuprimecamMapper(rerun=rerun)
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})
    frames = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=field,filter=filter))
    ccds = range(10)

    frames = [126968]
    ccds = [8]

    display = True

    for frameId in frames:
        for ccdId in ccds:
            if (frameId != 126963 or ccdId != 2):
                run(ioMgr, frameId, ccdId, display)
