import os
import math
import struct

# Most part of this code is translated from WCStools.
# wcstools-3.7.8/libwcs/ubcread.c
#
# Written  by N.Yasuda 2009/12/03 v1.0
# Modified by N.Yasuda 2009/12/24 v1.1
#          Add errors of proper motions and flags

class UBCstar:
    def __init__(self, name, b):
        self.name = name
        vals = struct.unpack('20i', b)
        self.rasec = vals[0]
        self.decsec = vals[1]
        self.pm = vals[2]
        self.pmerr = vals[3]
        self.poser = vals[4]
        self.mag0 = [0] * 5
        self.magerr = [0] * 5
        self.index = [0] * 5
        for i in range(5):
            self.mag0[i] = vals[5+i]
            self.magerr[i] = vals[10+i]
            self.index[i] = vals[15+i]

    def __cmp__(self, other):
        return cmp(self.mag(0), other.mag(0))
        #return cmp(self.dist, other.dist)

    def ra(self):                       # R.A. in degree
        return self.rasec / 360000.0

    def dec(self):                      # Decl. in degree
        return (self.decsec - 32400000) / 360000.0
        
    def rastr(self):                    # R.A. in string HH:MM:SS.SSS
        return ra2str(self.rasec / 360000.0)

    def decstr(self):                   # Decl. in string +DD:MM:SS.SS
        return dec2str((self.decsec - 32400000) / 360000.0)

    def era(self):                      # Error of R.A. in arcsec
        err = self.poser % 1000
        return err * 0.001

    def edec(self):                     # Error of Decl. in arcsec
        err = (self.poser % 1000000) / 1000
        return err * 0.001

    def epoch(self):                    # Mean epoch
        epoch = (self.poser % 1000000000) / 1000000
        return epoch * 0.1 + 1950.0
    
    def pra(self):                      # Proper motion on R.A. in arcsec/year
        pm = self.pm % 10000
        pm = pm * 0.002 - 10.0
        return pm

    def pdec(self):                     # Proper motion on Decl. in arcsec/year
        pm = (self.pm % 100000000) / 10000
        pm = pm * 0.002 - 10.0
        return pm

    def epra(self):                     # Error of proper motion on R.A. in arcsec/year
        epm = self.pmerr % 1000
        epm = epm * 0.001
        return epm
    
    def epdec(self):                    # Error of proper motion on Decl. in arcsec/year
        epm = (self.pmerr % 1000000) / 1000
        epm = epm * 0.001
        return epm
    
    def ndet(self):                     # Number of detections
        n = (self.pmerr % 1000000000) / 100000000
        return n

    def flags(self):                     
        Mflag = self.pm    / 1000000000 # Motion catalog flag
        sflag = self.pmerr / 1000000000 # Diffraction spike flag
        Yflag = self.poser / 1000000000 # YS4.0 correlation flag
        flags = ''
        if (Mflag):
            flags += 'M'
        else:
            flags += '.'
        if (sflag):
            flags += 's'
        else:
            flags += '.'
        if (Yflag):
            flags += 'Y'
        else:
            flags += '.'
        return flags
    
    def mag(self, i):                   # Magnitudes B1, R1, B2, R2, N
        xmag = (self.mag0[i] % 10000) * 0.01 #   i = 0,  1,  2,  3,  4
        if xmag == 0.0:
            xmag = 99.99
        return (xmag)

def degrad(x):
    return x*math.pi/180.0

def wcsdist(x1, y1, x2, y2):
    pos1 = []
    pos1.append(math.cos(degrad(x1)) * math.cos(degrad(y1)))
    pos1.append(math.sin(degrad(x1)) * math.cos(degrad(y1)))
    pos1.append(math.sin(degrad(y1)))
    pos2 = []
    pos2.append(math.cos(degrad(x2)) * math.cos(degrad(y2)))
    pos2.append(math.sin(degrad(x2)) * math.cos(degrad(y2)))
    pos2.append(math.sin(degrad(y2)))

    w = 0.0
    d1 = 0.0
    d2 = 0.0
    for i in range(3):
        w = w + (pos1[i] * pos2[i])
        d1 = d1 + (pos1[i] * pos1[i])
        d2 = d2 + (pos2[i] * pos2[i])
    diff = math.acos(w / (math.sqrt(d1) * math.sqrt(d2))) / math.pi * 180.0
    return diff
    
def RefLim (cra, cdec, dra, ddec, verbose):
    # Deal with all or nearly all of the sky
    if (ddec > 80.0 and dra > 150.0):
        ramin = 0.0
        ramax = 360.0
        decmin = -90.0
        decmax = 90.0
        wrap = 0
        if verbose:
            print "RefLim: RA: 0.0 - 360.0  Dec: -90.0 - 90.0"
        return ramin, ramax, decmin, decmax, wrap


    # Set declination limits for search
    dec1 = cdec - ddec
    dec2 = cdec + ddec

    # dec1 is always the smallest declination
    if (dec1 > dec2):
        dec = dec1
        dec1 = dec2
        dec2 = dec

    # Adjust width in right acension to that at max absolute declination
    adec1 = abs(dec1)
    adec2 = abs(dec2)
    if (adec1 > adec2):
        adec = adec1
    else:
        adec = adec2
    acdec = abs(cdec)
    if (adec < 90.0 and adec > acdec):
        dra = dra * (math.cos(degrad(acdec)) / math.cos(degrad(adec)))

    # Set right ascension limits for search
    ra1 = cra - dra
    ra2 = cra + dra

    # Keep right ascension limits between 0 and 360 degrees
    if (ra1 < 0.0):
        nrot = 1 - int(ra1 / 360.0)
        ra1 = ra1 + 360.0 * nrot
    if (ra1 > 360.0):
        nrot = int(ra1 / 360.0)
        ra1 = ra1 - 360.0 * nrot
    if (ra2 < 0.0):
        nrot = 1 - int(ra2 / 360.0)
        ra2 = ra2 + 360.0 * nrot
    if (ra2 > 360.0):
        nrot = int(ra2 / 360.0)
        ra2 = ra2 - 360.0 * nrot

    if (ra1 > ra2):
        wrap = 1
    else:
        wrap = 0

    ramin = ra1
    ramax = ra2

    decmin = dec1
    decmax = dec2

    # Check for pole
    dist = wcsdist(cra, cdec, ramax, decmax)
    if (cdec + dist > 90.0):
        ramin = 0.0
        ramax = 359.99999
        decmax = 90.0
    elif (cdec - dist < -90.0):
        ramin = 0.0
        ramax = 359.99999
        decmin = -90.0
    elif (decmin < -90.0):
        decmin = -90.0
        ramin = 0.0
        ramax = 359.99999
    elif (decmax > 90.0):
        decmax = 90.0
        ramin = 0.0
        ramax = 359.99999

    if verbose:
        print "RefLim: RA: %s - %s  Dec: %s - %s" % (ra2str(ramin), ra2str(ramax), dec2str(decmin), dec2str(decmax)),
        if (wrap == 1):
            print " wrap"
        else:
            print ""

    return ramin, ramax, decmin, decmax, wrap

def ra2str(ra):
    a = ra / 15.0
    h = int(math.floor(a))
    b = (a - h) * 60.0
    m = int(math.floor(b))
    s = (b - m) * 60

    return '%02d:%02d:%06.3f' % (h, m, s)

def dec2str(dec):
    sign = 1
    if (dec < 0.0):
        sign = -1
    a = sign * dec
    d = int(math.floor(a))
    b = (a - d) * 60.0
    m = int(math.floor(b))
    s = (b - m) * 60.0

    if (sign == 1):
        return '+%02d:%02d:%05.2f' % (d, m, s)
    else:
        return '-%02d:%02d:%05.2f' % (d, m, s)

class UBCread:
    upath = "/data/USNO-B1.0"
    NZONES = 1800
    nbent = 80
    zonesize = 0.1                      # number of degrees per declination zone

    def ubczone (self, dec):

        zone = int((dec + 90.0) / self.zonesize)
        if (zone > 1799):
            zone = 1799
        elif (zone < 0):
            zone = 0

        return zone

    def ubczones (self, ra1, ra2, dec1, dec2, verbose = False):
        iz1 = self.ubczone(dec1)
        iz2 = self.ubczone(dec2)

        zones = []
        if (iz2 >= iz1):
            for iz in range(iz1, iz2+1):
                zones.append(iz)
        else:
            for iz in range(iz2, iz1+1):
                zones.append(iz)
            
        if verbose:
            print "UBCZONES: %d zones: %d - %d" % (len(zones), zones[0], zones[len(zones)-1])
            print "UBCZONES: RA: %.5f - %.5f, Dec: %.5f - %.5f" % (ra1, ra2, dec1, dec2)

        return zones

    def ubcpath(self, zn):
        # Return error code and null path if zone is out of range
        if (zn < 0 or zn > 1799):
            print "UBCPATH: zone %d out of range 0-1799" % (zn)
            return ""

        # Set path for USNO-B1.0 zone catalog
        path = "%s/%03d/b%04d.cat" % (self.upath, zn/10, zn)
        return path

    def ubcopen(self, znum):
        # Get path to zone catalog
        zonepath = self.ubcpath(znum)
        if (zonepath == ""):
            print "UBCOPEN: Cannot find zone catalog for %d" % (znum)
            return 0

        # Find number of stars in zone catalog by its length
        lfile = os.path.getsize(zonepath)
        if (lfile < 2):
            print "UB zone catalog %s has no entries"
            return 0
        else:
            nstars = lfile / self.nbent

        # Open zone catalog
        self.fcat = open(zonepath, "rb")
    
        return nstars

    def ubcsra(self, znum, ra1, ra2):
        path = "%s/%03d/b%04d.acc" % (self.upath, znum/10, znum)
        f = open(path, "rt")
        ra = [0.0] * 96
        ncum = [0] * 96
        nbin = [0] * 96
        i = 0
        for line in f:
            itemList = line[:-1].split()
            ra[i] = float(itemList[0])
            ncum[i] = int(itemList[1])
            nbin[i] = int(itemList[2])
            i = i + 1
        f.close()
    
        istar1 = -1
        istar2 = -1
        for i in range(96):
            if (istar1 == -1 and ra[i] > ra1 / 15.):
                istar1 = ncum[i-1]
            if (istar2 == -1 and ra[i] > ra2 / 15.):
                istar2 = ncum[i]

        if istar1 == -1 and istar2 == -1:
            istar1 = ncum[95]
            istar2 = ncum[95] + nbin[95]

        return istar1, istar2

    def __init__(self, cra, cdec, rad, verbose=False):
        # cra  : center R.A. (degree)
        # cdec : center Decl. (degree)
        # rad  : search radius (arcsec)
        self.cra  = cra
        self.cdec = cdec
        self.rad  = rad
        self.verbose = verbose
        self.dra  = self.rad / 3600.0 / math.cos(degrad(self.cdec))
        self.ddec = self.rad / 3600.0

        # Find RA and Dec limits in catalog coordinate system
        rra1, rra2, rdec1, rdec2, wrap = RefLim(self.cra, self.cdec, self.dra, self.ddec, self.verbose)
       
        # Find declination zones to search
        zones = self.ubczones(rra1, rra2, rdec1, rdec2, self.verbose)

        # Convert RA and Dec limits to same units as catalog for quick filter
        ubra1 = int(rra1 * 360000.0 + 0.5);
        ubra2 = int(rra2 * 360000.0 + 0.5);
        ubdec1 = int((rdec1 * 360000.0) + 32400000.5)
        ubdec2 = int((rdec2 * 360000.0) + 32400000.5)

        # Loop through region list
        self.stars = []
        for iz in range(len(zones)):
            # Get path to zone catalog
            znum = zones[iz]
            nstars = self.ubcopen(znum)
            if (nstars != 0):
                for iwrap in range(wrap+1):
                    if (wrap == 0):
                        istar1, istar2 = self.ubcsra(znum, rra1, rra2)
                    else:
                        if (iwrap == 0):
                            istar1, istar2 = self.ubcsra(znum, rra1, 360.0)
                        else:
                            istar1, istar2 = self.ubcsra(znum, 0.0, rra2)
                    #print istar1, istar2
                    nread = istar2 - istar1
                    self.fcat.seek((istar1-1)*self.nbent)
                    
                    # Loop through zone catalog for this region
                    for i in range(nread):
                        name = "%04d-%07d" % (znum, istar1+i)
                        b = self.fcat.read(self.nbent)
                        star = UBCstar(name, b)

                        # Extract selected fields

                        # Check rough position limits
                        if ((star.decsec >= ubdec1 and star.decsec <= ubdec2) and
                            ((wrap and (star.rasec >= ubra1 or star.rasec <= ubra2)) or
                             (not (wrap) and (star.rasec >= ubra1 and star.rasec <= ubra2)))):
                            ra = star.rasec / 360000.0
                            dec = (star.decsec - 32400000) / 360000.0

                            # Test spatial limits
                            if (wcsdist(self.cra, self.cdec, ra, dec) < self.ddec):
                                star.dist = wcsdist(cra, cdec, ra, dec) * 3600.0
                                self.stars.append(star)

            self.fcat.close()

        self.stars.sort()

    def __init__(self, cra, cdec, dra, ddec, verbose=False):
        # cra  : center R.A. (degree)
        # cdec : center Decl. (degree)
        # rad  : search radius (arcsec)
        self.cra  = cra
        self.cdec = cdec
        self.verbose = verbose
        self.dra  = dra / math.cos(degrad(self.cdec))
        self.ddec = ddec

        # Find RA and Dec limits in catalog coordinate system
        rra1, rra2, rdec1, rdec2, wrap = RefLim(self.cra, self.cdec, self.dra, self.ddec, self.verbose)
       
        # Find declination zones to search
        zones = self.ubczones(rra1, rra2, rdec1, rdec2, self.verbose)

        # Convert RA and Dec limits to same units as catalog for quick filter
        ubra1 = int(rra1 * 360000.0 + 0.5);
        ubra2 = int(rra2 * 360000.0 + 0.5);
        ubdec1 = int((rdec1 * 360000.0) + 32400000.5)
        ubdec2 = int((rdec2 * 360000.0) + 32400000.5)

        # Loop through region list
        self.stars = []
        for iz in range(len(zones)):
            # Get path to zone catalog
            znum = zones[iz]
            nstars = self.ubcopen(znum)
            if (nstars != 0):
                for iwrap in range(wrap+1):
                    if (wrap == 0):
                        istar1, istar2 = self.ubcsra(znum, rra1, rra2)
                    else:
                        if (iwrap == 0):
                            istar1, istar2 = self.ubcsra(znum, rra1, 360.0)
                        else:
                            istar1, istar2 = self.ubcsra(znum, 0.0, rra2)
                    #print istar1, istar2
                    nread = istar2 - istar1
                    self.fcat.seek((istar1-1)*self.nbent)
                    
                    # Loop through zone catalog for this region
                    for i in range(nread):
                        name = "%04d-%07d" % (znum, istar1+i)
                        b = self.fcat.read(self.nbent)
                        star = UBCstar(name, b)

                        # Extract selected fields

                        # Check rough position limits
                        if ((star.decsec >= ubdec1 and star.decsec <= ubdec2) and
                            ((wrap and (star.rasec >= ubra1 or star.rasec <= ubra2)) or
                             (not (wrap) and (star.rasec >= ubra1 and star.rasec <= ubra2)))):
                            ra = star.rasec / 360000.0
                            dec = (star.decsec - 32400000) / 360000.0

                            # Test spatial limits
                            if (ra  > cra  - 0.5 * dra  and ra  < cra  + 0.5 * dra and
                                dec > cdec - 0.5 * ddec and dec < cdec + 0.5 * ddec):
                                star.dist = wcsdist(cra, cdec, ra, dec) * 3600.0
                                self.stars.append(star)

            self.fcat.close()

        self.stars.sort()

if __name__ == '__main__':
    #ubread = UBCread(359.97, 30.0, 180.0)
    #ubread = UBCread(150.0, 30.0, 180.0)
    #ubread = UBCread(34.5, -5.0, 180.0)
    ubread = UBCread(352.034, -2.922, 0.114911476026, 0.234367790703)
    for star in ubread.stars:
        print "%s %10.6f %+10.6f %5.3f %5.3f %6.1f %6.3f %6.3f %5.3f %5.3f %d %s %6.2f %6.2f %6.2f %6.2f %6.2f %7.2f" % \
        (star.name, star.ra(), star.dec(), star.era(), star.edec(),
         star.epoch(), star.pra(), star.pdec(), star.epra(), star.epdec(),
         star.ndet(), star.flags(),
         star.mag(0), star.mag(1), star.mag(2), star.mag(3), star.mag(4),
         star.dist)
