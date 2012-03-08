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
#include "lsst/afw/table.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/TanWcs.h"
%}

%include "lsst/p_lsstSwig.i"

%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/wcs.i"

%include "match.i"
%include "sipfit.i"
