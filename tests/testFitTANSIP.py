#!/usr/bin/env python

# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
 
"""
Tests for TanSip

Run with:
   testFitTANSIP.py
"""
import unittest

if False:
    import matplotlib
    matplotlib.use('Agg')
    import pylab
    pnum = 1
import numpy

import lsst.utils.tests as tests
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import hsc.meas.astrom as hscAstrom

from lsst.pex.logging import Log
log = Log.getDefaultLog()
#log.setThreshold(Log.DEBUG)


class FitTANSIPTestCase(unittest.TestCase):
    """A test case for fitTANSIP
    """
    def setUp(self):
        self.crCoord = afwCoord.IcrsCoord(afwGeom.PointD(44., 45.))
        self.crPix = afwGeom.PointD(0, 0)
        
        arcsecPerPixel = 1/3600.0
        CD11 = arcsecPerPixel
        CD12 = 0
        CD21 = 0
        CD22 = arcsecPerPixel
        
        self.tanWcs = afwImage.makeWcs(self.crCoord, self.crPix, CD11, CD12, CD21, CD22)

        S = 3000
        N = 4

        posRefSchema = afwTable.SimpleTable.makeMinimalSchema()
        posRefTable = afwTable.SimpleTable.make(posRefSchema)
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        srcCentroidKey = srcSchema.addField("centroid", type="PointD")
        srcSchema.addField("centroid.flags", type="Flag")
        srcCentroidErrKey = srcSchema.addField("centroid.err", type="CovPointF")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        self.matchList = afwTable.ReferenceMatchVector()

        zeroCov = numpy.zeros([2, 2], dtype=numpy.float32)
        for i in range(2):
            zeroCov[i,i] = 1
        for i in numpy.linspace(0., S, N):
            for j in numpy.linspace(0., S, N):
                src = srcTable.makeRecord()
                refObj = posRefTable.makeRecord()

                src.set(srcCentroidKey.getX(), i)
                src.set(srcCentroidKey.getY(), j)
                src.set(srcCentroidErrKey, zeroCov)

                c = self.tanWcs.pixelToSky(i, j);
                refObj.setCoord(c);
                
                if False:
                    print "x,y = (%.1f, %.1f) pixels -- RA,Dec = (%.3f, %.3f) deg" %  (
                        i, j, c.toFk5().getRa().asDegrees(),
                        c.toFk5().getDec().asDegrees())

                self.matchList.push_back(afwTable.ReferenceMatch(refObj, src, 0.0))

    def tearDown(self):
        del self.matchList
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testOffset(self):
        """Add an offset"""
        self.doTest("testOffset", lambda x, y: (x + 5, y + 7))

    def testLinearX(self):
        """Scale x, offset y"""
        self.doTest("testLinearX", lambda x, y: (2*x, y + 7))

    def testLinearXY(self):
        """Scale x and y"""
        self.doTest("testLinearXY", lambda x, y: (2*x, 3*y))

    def testLinearYX(self):
        """Add an offset to each point; scale in y and x"""
        self.doTest("testLinearYX", lambda x, y: (x + 0.2*y, y + 0.3*x))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        self.doTest("testQuadraticX", lambda x, y: (x + 1e-5*x**2, y), order=3)
        
    def assertAlmostEqualAngle(self, a1, a2, tol=3):
        self.assertAlmostEqual(a1.asArcseconds(), a2.asArcseconds(), tol)

    def checkResults(self, tanSipWcs):
        for refObj, src, d in self.matchList:
            srcX = src.getX()
            srcY = src.getY()
            refRa  = refObj.getRa();
            refDec = refObj.getDec();
            srcCoord = tanSipWcs.pixelToSky(srcX, srcY);
            srcCoord = srcCoord.toIcrs()

            if False:
                print "ref RA,Dec = (%.8f, %.8f) deg" % (refRa.asDegrees(), refDec.asDegrees())
                print "src RA,Dec = (%.8f, %.8f) deg" % (srcCoord.getRa().asDegrees(), srcCoord.getDec().asDegrees())
            self.assertAlmostEqualAngle(refRa, srcCoord.getRa())
            self.assertAlmostEqualAngle(refDec, srcCoord.getDec())

            # these are in pixels.
            refObjXY = tanSipWcs.skyToPixel(refObj.getCoord());
            if False:
                print "ref X,Y = (%.3f, %.3f)" % (refObjXY[0], refObjXY[1])
                print "src X,Y = (%.3f, %.3f)" % (srcX, srcY);
            self.assertAlmostEqual(srcX, refObjXY[0], 3) # within a milli-pixel
            self.assertAlmostEqual(srcY, refObjXY[1], 3)

    def doTest(self, name, func, order=2):
        """Apply func(x, y) to each source in self.srcCat, then fit and check the resulting WCS
        """
        print
        print name

        for refObj, src, d in self.matchList:
            x, y = func(src.getX(), src.getY())
            src.set("centroid.x", x); src.set("centroid.y", y)

        tanSipWcs = hscAstrom.fitTANSIP(order, self.matchList, self.crCoord, self.crPix)

        if False:
            md = tanSipWcs.getFitsMetadata()
            for name in md.names():
                print "%s: %s" % (name, md.get(name))

        if False:
            self.doPlot(tanSipWcs)

        self.checkResults(tanSipWcs)

    def doPlot(self, tanSipWcs):
        """Plot the reference objects, sources and WCS"""
        xs,ys, xc,yc = [],[],[],[]
        rs,ds, rc,dc = [],[],[],[]
        for refObj, src, d in self.matchList:
            xs.append(src.getX())
            ys.append(src.getY())
            refObjXY = tanSipWcs.skyToPixel(refObj.getCoord())
            xc.append(refObjXY[0])
            yc.append(refObjXY[1])
            rc.append(refObj.getRa())
            dc.append(refObj.getDec())
            srd = tanSipWcs.pixelToSky(src.getX(), src.getY()).toFk5()
            rs.append(srd.getRa())
            ds.append(srd.getDec())
        xs = numpy.array(xs)
        ys = numpy.array(ys)
        xc = numpy.array(xc)
        yc = numpy.array(yc)
            
        global pnum
        pylab.clf()
        pylab.plot(xs, ys, 'r.')
        pylab.plot(xc, yc, 'bx')
        fn = 'check-%i.png' % pnum
        pylab.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        pylab.clf()
        pylab.plot(xs, xc-xs, 'b.')
        fn = 'check-%i.png' % pnum
        pylab.xlabel('x(source)')
        pylab.ylabel('x(ref - src)')
        pylab.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        pylab.clf()
        pylab.plot(rs, ds, 'r.')
        pylab.plot(rc, dc, 'bx')
        fn = 'check-%i.png' % pnum
        pylab.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        pylab.clf()
        for y in [0., 1000., 2000., 3000., 4000.]:
            x0,y0 = [],[]
            x1,y1 = [],[]
            for x in numpy.linspace(0., 4000., 400):
                x0.append(x)
                y0.append(y)
                rd = tanSipWcs.pixelToSky(x, y)
                xy = tanSipWcs.skyToPixel(rd)
                x1.append(xy[0])
                y1.append(xy[1])
            x0 = numpy.array(x0)
            x1 = numpy.array(x1)
            pylab.plot(x0, x1-x0, 'b-')
        fn = 'check-%i.png' % pnum
        pylab.savefig(fn)
        print 'Wrote', fn
        pnum += 1

        pylab.clf()
        for y in [0., 1000., 2000., 3000., 4000.]:
            x0,y0 = [],[]
            x1,y1 = [],[]
            x2,y2 = [],[]
            for x in numpy.linspace(0., 4000., 400):
                x0.append(x)
                y0.append(y)
                xy = tanSipWcs.undistortPixel(afwGeom.Point2D(x+1,y+1))
                x1.append(xy[0]-1)
                y1.append(xy[1]-1)
                xy = tanSipWcs.distortPixel(xy)
                x2.append(xy[0]-1)
                y2.append(xy[1]-1)
            x0 = numpy.array(x0)
            x1 = numpy.array(x1)
            x2 = numpy.array(x2)
            pylab.plot(x0, x1-x0, 'b-')
            pylab.plot(x0, x1-x2, 'r-')
            pylab.plot(x0, (x2-x0)*1e5, 'g-')
        pylab.xlabel('x (orig)')
        pylab.ylabel('dx (undistorted)')
        fn = 'check-%i.png' % pnum
        pylab.savefig(fn)
        print 'Wrote', fn
        pnum += 1

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(FitTANSIPTestCase)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
