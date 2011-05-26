// -*- lsst-c++ -*-
%define matchLib_DOCSTRING
"
Python interface to hsc::meas::match
"
%enddef

%feature("autodoc", "1");
%module(package="hsc.meas.match", docstring=matchLib_DOCSTRING) matchLib

%{
#include "lsst/afw/image.h"
%}

%include "lsst/p_lsstSwig.i"

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: ssh://hsc-gw2.mtk.nao.ac.jp/ana/hgrepo/hscAstrom/python/hsc/meas/match/matchLib.i $"):
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
