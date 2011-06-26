import os, sys
import lsst.pipette.readwrite           as pipReadWrite
import lsst.obs.suprimecam              as obsSc
import lsst.afw.image                   as afwImage
import lsst.afw.detection               as afwDet
import lsst.afw.display.ds9             as ds9
import hsc.meas.astrom                  as hscAst
import lsst.meas.algorithms             as measAlg
import lsst.pipette.util as pipUtil
import lsst.pipette.distortion as pipDist
import math
import datetime
import lsst.pex.policy as policy
import lsst.meas.astrom as measAst
from lsst.pex.logging import Log
import lsst.pipette.config as pipConfig
import lsst.afw.geom as afwGeom

def goodStar(s):
    return s.getXAstrom() == s.getXAstrom() and \
           not (s.getFlagForDetection() & measAlg.Flags.SATUR)

def run(ioMgr, frameId, ccdId, config, display=False):

    data = {'visit': frameId, 'ccd':ccdId}

    raws = ioMgr.readRaw(data)[0]
    wcsIn = raws.getWcs()

    exposure = ioMgr.read('calexp', data, ignore=True)[0]
    wcsTrue = exposure.getWcs()
    exposure.setWcs(wcsIn)

    if display:
        ds9.mtv(exposure, frame=1)

    bscSet = ioMgr.read('bsc', data, ignore=True)[0].getSources()
    srcSet = [s for s in bscSet if goodStar(s)]
    print "# of srcSet: ", len(srcSet)

    distConfig = config['distortion']
    ccd = pipUtil.getCcd(exposure)
    dist = pipDist.createDistortion(ccd, distConfig)

    distSrc = dist.actualToIdeal(srcSet)

    if display:
        for (s, d) in zip(srcSet, distSrc):
            xs = s.getXAstrom()
            ys = s.getYAstrom()
            xd = d.getXAstrom()
            yd = d.getYAstrom()
            #ds9.dot("o", xs, ys, size=10, frame=1)
            #ds9.dot("o", xd, yd, size=10, frame=1, ctype=ds9.RED)
            #ds9.line([[xs,ys], [xd,yd]], frame=1)

    xMin, xMax, yMin, yMax = float("INF"), float("-INF"), float("INF"), float("-INF")
    for x, y in ((0.0, 0.0), (0.0, exposure.getHeight()), (exposure.getWidth(), 0.0),
                 (exposure.getWidth(), exposure.getHeight())):
        point = afwGeom.Point2D(x, y)
        point = dist.actualToIdeal(point)
        x, y = point.getX(), point.getY()
        if x < xMin: xMin = x
        if x > xMax: xMax = x
        if y < yMin: yMin = y
        if y > yMax: yMax = y
    xMin = int(xMin)
    yMin = int(yMin)
    size = (int(xMax - xMin + 0.5),
            int(yMax - yMin + 0.5))
    for s in distSrc:
        s.setXAstrom(s.getXAstrom() - xMin)
        s.setYAstrom(s.getYAstrom() - yMin)

    sortedDistSrc = sorted(distSrc, key=lambda x:x.getPsfFlux(), reverse=True)

    if display:
        for i, d in enumerate(sortedDistSrc):
            #ds9.dot('o', d.getXAstrom(), d.getYAstrom(), size=15, frame=1)
            if i >= 99: break
        ds9.line([[0, 0], [0+size[0], 0]], frame=1)
        ds9.line([[0+size[0], 0], [0+size[0], 0+size[1]]], frame=1)
        ds9.line([[0+size[0], 0+size[1]], [0, 0+size[1]]], frame=1)
        ds9.line([[0, 0+size[1]], [0, 0]], frame=1)

    log = Log.getDefaultLog()
    pol = policy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAst.createSolver(pol, log)

    (W,H) = size
    filterName = 'r'
    idName = 'id'

    wcsIn.shiftReferencePixel(-xMin, -yMin)

    catSet = hscAst.queryReferenceCatalog(solver, distSrc, wcsIn, (W,H),
                                          filterName, idName)

    print "# of catSet: ", len(catSet)

    if display:
        for c in catSet:
            ds9.dot("o", c.getXAstrom()+xMin, c.getYAstrom()+yMin, size=15, ctype=ds9.RED, frame=1)

    start = datetime.datetime.today()

    matchList = hscAst.match(distSrc, catSet)

    end = datetime.datetime.today()

    wcsIn.shiftReferencePixel(xMin, yMin)

    matchList2 = []
    for m in matchList:
        for s in srcSet:
            if (s.getId() == m.second.getId()):
                x0 = m.first.getXAstrom() + xMin
                y0 = m.first.getYAstrom() + yMin
                x1 = m.second.getXAstrom()
                y1 = m.second.getYAstrom()
                x1 = s.getXAstrom()
                y1 = s.getYAstrom()
                #print x0, y0, x1, y1, m.second.getId()
                if display:
                    ds9.dot("o", x1, y1, size=10, frame=1)
                    ds9.line([[x0, y0], [x1, y1]], ctype=ds9.RED, frame=1)
                m2 = afwDet.SourceMatch(m.first, s, 0.0)
                matchList2.append(m2)

    elapsedtime = end - start

    print "%6d %d %3d %6.2f" % (frameId, ccdId, len(matchList), elapsedtime.seconds + elapsedtime.microseconds/1000000.)
#    f = open('resultlog', 'a')
#    f.write("%6d %d %3d %6.2f\n" % (frameId, ccdId, len(matchList), elapsedtime.seconds + elapsedtime.microseconds/1000000.))
#    f.close()

if __name__ == '__main__':
    field = "DITHER3"
    filter = "W-S-R+"
    rerun = "yasuda"

    policyfile = os.path.join(os.getenv("PIPETTE_DIR"), "policy", "suprimecam.paf")
    config = pipConfig.Config(policyfile)

    mapper = obsSc.SuprimecamMapper(rerun=rerun)
    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})
    frames = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=field,filter=filter))
    ccds = range(10)

    frames = [131658]
    ccds = [6]

    display = True

    for frameId in frames:
        for ccdId in ccds:
            run(ioMgr, frameId, ccdId, config, display)
