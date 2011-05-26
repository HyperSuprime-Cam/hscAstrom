// -*- lsst-c++ -*-
%define fitLib_DOCSTRING
"
Python interface to hsc::meas::fit
"
%enddef

%feature("autodoc", "1");
%module(package="hsc.meas.fit", docstring=fitLib_DOCSTRING) fitLib

%{
#include "lsst/afw/image.h"
%}

%include "lsst/p_lsstSwig.i"

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: ssh://hsc-gw2.mtk.nao.ac.jp/ana/hgrepo/hscAstrom/python/hsc/meas/fit/fitLib.i $"):
    version_svn = lsst.utils.guessSvnVersion(HeadURL)

    try:
	import eups
    except ImportError:
        return version_svn
    else:
	try:
	    version_eups = eups.setup("hscAstrom")
	except AttributeError:
	    return version_svn

    if version_eups == version_svn:
        return version_svn
    else:
	return "%s (setup: %s)" % (version_svn, version_eups)
%}

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/detection/detectionLib.i"

%include "fit.i"
