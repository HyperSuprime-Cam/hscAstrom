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
Tests for fitTANSIP

Run with:
   testFitTANSIP.py
"""
import unittest

import numpy
try:
    import matplotlib
    matplotlib.use("Agg")
    import pylab
except ImportError:
    pass

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

        refTable = afwTable.SimpleTable.make(afwTable.SimpleTable.makeMinimalSchema())
        srcSchema = afwTable.SourceTable.makeMinimalSchema()
        self.srcCentroidKey = srcSchema.addField("centroid", type="PointD")
        srcSchema.addField("centroid.flags", type="Flag")
        self.srcCentroidErrKey = srcSchema.addField("centroid.err", type="CovPointF")
        srcTable = afwTable.SourceTable.make(srcSchema)
        srcTable.defineCentroid("centroid")
        self.matchList = afwTable.ReferenceMatchVector()

        for i in numpy.linspace(0., S, N):
            for j in numpy.linspace(0., S, N):
                src = srcTable.makeRecord()
                refObj = refTable.makeRecord()

                src.set(self.srcCentroidKey.getX(), i)
                src.set(self.srcCentroidKey.getY(), j)
                src.set(self.srcCentroidErrKey[0, 0], 0.1)
                src.set(self.srcCentroidErrKey[1, 1], 0.1)

                c = self.tanWcs.pixelToSky(afwGeom.Point2D(i, j))
                refObj.setCoord(c)
                
                if False:
                    print "x,y = (%.1f, %.1f) pixels -- RA,Dec = (%.3f, %.3f) deg" % \
                        (i, j, c.toFk5().getRa().asDegrees(), c.toFk5().getDec().asDegrees())

                self.matchList.push_back(afwTable.ReferenceMatch(refObj, src, 0.0))

    def tearDown(self):
        del self.matchList
        del self.tanWcs

    def testTrivial(self):
        """Add no distortion"""
        self.doTest("testTrivial", lambda x, y: (x, y))

    def testScaleXY(self):
        """Scale x and y slightly"""
        self.doTest("testScaleXY", lambda x, y: (1.002*x, 0.998 * y))

    def testQuadraticX(self):
        """Add quadratic distortion in x"""
        self.doTest("testQuadraticX", lambda x, y: (x + 1e-5*x**2, y), order=3)

    def checkResults(self, tanSipWcs):
        for refObj, src, d in self.matchList:
            srcPixPos = afwGeom.Point2D(
                src.get(self.srcCentroidKey.getX()),
                src.get(self.srcCentroidKey.getY()),
            )
            refCoord = refObj.get("coord")
            srcCoord = tanSipWcs.pixelToSky(srcPixPos)
            srcCoord = srcCoord.toIcrs()

            if False:
                print "ref RA,Dec = (%.8f, %.8f) deg" % (refCoord.getRa().asDegrees(), refCoord.getDec().asDegrees())
                print "src RA,Dec = (%.8f, %.8f) deg" % (srcCoord.getRa().asDegrees(), srcCoord.getDec().asDegrees())
            self.assertLess(refCoord.angularSeparation(srcCoord), 0.001 * afwGeom.arcseconds)

            refPixPos = tanSipWcs.skyToPixel(refCoord)
            if False:
                print "ref X,Y = (%.3f, %.3f)" % (refPixPos[0], refPixPos[1])
                print "src X,Y = (%.3f, %.3f)" % (srcPixPos[0], srcPixPos[1])
            self.assertAlmostEqual(srcPixPos[0], refPixPos[0], 3) # within a milli-pixel
            self.assertAlmostEqual(srcPixPos[1], refPixPos[1], 3)

    def doTest(self, name, func, order=2):
        """Apply func(x, y) to each source in self.srcCat, then fit and check the resulting WCS
        """
        for refObj, src, d in self.matchList:
            x0 = src.get(self.srcCentroidKey.getX())
            y0 = src.get(self.srcCentroidKey.getY())
            x, y = func(x0, y0)
            src.set(self.srcCentroidKey.getX(), x)
            src.set(self.srcCentroidKey.getY(), y)

        tanSipWcs = hscAstrom.fitTANSIP(order, self.matchList, self.crCoord, self.crPix)

        if False:
            self.plotWcs(tanSipWcs, name=name)
        
        self.checkResults(tanSipWcs)

    def plotWcs(self, tanSipWcs, name=""):
        pnum = 1
        fileNamePrefix = "testFitTANSIP_%s" % (name,)

        xs,ys, xc,yc = [],[],[],[]
        rs,ds, rc,dc = [],[],[],[]
        for ref, src, d in self.matchList:
            xs.append(src.getX())
            ys.append(src.getY())
            refPixPos = tanSipWcs.skyToPixel(ref.getCoord())
            xc.append(refPixPos[0])
            yc.append(refPixPos[1])
            rc.append(ref.getRa())
            dc.append(ref.getDec())
            srd = tanSipWcs.pixelToSky(src.getX(), src.getY()).toFk5()
            rs.append(srd.getRa())
            ds.append(srd.getDec())
        xs = numpy.array(xs)
        ys = numpy.array(ys)
        xc = numpy.array(xc)
        yc = numpy.array(yc)
            
        pylab.clf()
        pylab.plot(xs, ys, "r.")
        pylab.plot(xc, yc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        pylab.plot(xs, xc-xs, "b.")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.xlabel("x(source)")
        pylab.ylabel("x(ref - src)")
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        pylab.plot(rs, ds, "r.")
        pylab.plot(rc, dc, "bx")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
        pnum += 1

        pylab.clf()
        for y in numpy.linspace(0, 4000, 5):
            x0,y0 = [],[]
            x1,y1 = [],[]
            for x in numpy.linspace(0., 4000., 401):
                x0.append(x)
                y0.append(y)
                rd = tanSipWcs.pixelToSky(x, y)
                xy = tanSipWcs.skyToPixel(rd)
                x1.append(xy[0])
                y1.append(xy[1])
            x0 = numpy.array(x0)
            x1 = numpy.array(x1)
            pylab.plot(x0, x1-x0, "b-")
        fileName = "%s_%i.png" % (fileNamePrefix, pnum)
        pylab.savefig(fileName)
        print "Wrote", fileName
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
