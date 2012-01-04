// -*- lsst-c++ -*-
%define astromLib_DOCSTRING
"
Python interface to hsc::meas::match
"
%enddef

%feature("autodoc", "1");
%module(package="hsc.meas.astrom", docstring=astromLib_DOCSTRING) astromLib

%{
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/AperturePhotometry.h"
%}

%include "lsst/p_lsstSwig.i"

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: ssh://hsc-gw2.mtk.nao.ac.jp/ana/hgrepo/hscAstrom/python/hsc/meas/astrom/astromLib.i $"):
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

%import "lsst/afw/detection/detectionLib.i"

%include "match.i"
%include "sipfit.i"
