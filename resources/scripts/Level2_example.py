#!/usr/bin/env python

# This example script is for converting level 2 2008 data to IceEvent structures
# M.R.Duvoort 2009

import os, sys
from stat import *
from os.path import expandvars
from math import pi
from optparse import OptionParser
from I3Tray import *

usage = "usage: %prog [options]"

# option parsing
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="",type = "string",
                  dest="OUTFILE", help="write output to OUTFILE (.IceEvent.root format)")
parser.add_option("-i", "--infiles", default = "", type = "string",
                  dest="INFILE", help="read from these input files, seperated by commas")
parser.add_option("-g", "--gcd", default = "", type = "string",
                  dest="GCDFILE", help="read the GCD frames from this file")
parser.add_option("-n","--nevents",default="-1",type="int",
                  dest="NEVENTS", help="the number of events to be processed [default: all events]")

(options,args) = parser.parse_args()

gcdfile = options.GCDFILE
outfile = options.OUTFILE
nevt = options.NEVENTS
infiles = [""] #the input files is a vector of strings
infiles = options.INFILE.split(",") #split the string into a vector of strings
infiles.insert(0,gcdfile) #the GCDfile comes first!!!

##################################################################
# LOAD LIBS, CREATE A TRAY
##################################################################

load("libicetray")
load("libdataclasses")
load("libdataio")
load("libDOMcalibrator")
load("libIceDwalk")
load("libIcePacker")

tray = I3Tray()

##################################################################
# SERVICES
##################################################################

# read input files
tray.AddService("I3ReaderServiceFactory", "reader")(
    ("FileNameList",infiles),
    )

##################################################################
# MODULES
##################################################################

tray.AddModule("I3Muxer","muxer")

tray.AddModule("I3DOMcalibrator", "calibrate-inice")(
    ("CalibrationMode", 0),
    ("OutputToFile", False),
    ("InputRawDataName", "InIceRawData"),
    ("CorrectPedestalDroopDualTau", False),
    ("OutputATWDDataName", "InIceCalibratedATWD"),
    ("OutputFADCDataName", "InIceCalibratedFADC")
    )

tray.AddModule("I3DOMcalibrator", "calibrate-icetop")(
    ("CalibrationMode", 0),
    ("OutputToFile", False),
    ("InputRawDataName", "IceTopRawData"),
    ("CorrectPedestalDroopDualTau", False),
    ("OutputATWDDataName", "IceTopCalibratedATWD"),
    ("OutputFADCDataName", "IceTopCalibratedFADC")
    )

#Standard IceDwalk fit
tray.AddModule("I3IceDwalkFit","IceDwalk")(
    ("detector","I"),
    ("asqualitytypei",3),
    ("dminhiti",75),
    ("dtmargi",0),
    ("mergejetangi",7.5),
    ("mergejetdisti",30),
    ("jetinvoli",True),
    ("jetiteratei",True),
    ("maxdhiti",3.07126),
    ("maxhitsperdom",1),
    ("maximumdom",9999),
    ("minimumdom",0),
    ("mergecandangi",15),
    ("mergecanddisti",20),
    ("trackinvoli",True),
    ("usevgroupi",1),
    #  ("inputhitreadout","InitialHitSeriesReco"),
    ("inputpulsereadout","InitialPulseSeriesReco"),
    ("name","IceDwalkFit"),
    ("multisolutions",True)
    )

# The IcePack converter
tray.AddModule("IcePacker","icepacker")(
    ("RootFileName",outfile),
    ("TrackExtract",True),
    ("DOMmapInIceName","InIceRawData"),
    ("DOMmapIceTopName","IceTopRawData"),
    ("WaveformmapInIceATWDName","InIceCalibratedATWD"),
    ("WaveformmapInIceFADCName","InIceCalibratedFADC"),
    ("WaveformmapIceTopATWDName","IceTopCalibratedATWD"),
    ("WaveformmapIceTopFADCName","IceTopCalibratedFADC"),
    ("InIcePulseSeriesName","InitialPulseSeriesReco"),
    ("IceTopPulseSeriesName","IceTopVemPulseSeries"),
#    ("InIcePulseSeriesNames",["InitialPulseSeriesReco","InitialSLCPulseSeriesReco"]), # If this parameter is set, the above single pulse-extraction is overruled
#    ("IceTopPulseSeriesNames",["IceTopVemPulseSeries","IceTopSLCPulseSeries"]),# If this parameter is set, the above single pulse-extraction is overruled
    ("MCTreeName","I3MCTree"),
    ("MCWeightInfo","I3MCWeightDict"),
    ("EventHeaderName","I3EventHeader"),
    ("TriggerName","I3TriggerHierarchy"),
    ("setmaxtreesize",200000000)
    )

# count events (and stop when we've enough)
tray.AddModule("I3EventCounter","nprocessed")(
    ("NEvents",nevt),
    ("CounterStep",10),
    )

# from bits thou art made, to bits thou shalt return
tray.AddModule("TrashCan","destination")
    
# DO IT
tray.Execute()
tray.Finish()
