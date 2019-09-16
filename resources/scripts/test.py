#!/usr/bin/env python

import os, sys
from stat import *
from os.path import expandvars
from math import pi
from optparse import OptionParser
from I3Tray import *

infilename = expandvars("/net/local/icecube/I3_PORTS/RHEL4-x86_64/test-data/2007data/2007_I3Only_Run109732_Nch20.i3.gz")

##################################################################
# LOAD LIBS, CREATE A TRAY
##################################################################

load("libicetray")
load("libdataclasses")
load("libdataio")
load("libDOMcalibrator")
load("libFeatureExtractor")
load("libIceDwalk")
load("liblilliput")
load("libgulliver")
load("libgulliver-modules")
load("libparaboloid")
load("libIcePacker")

tray = I3Tray()

##################################################################
# SERVICES
##################################################################

# read input files
tray.AddService("I3ReaderServiceFactory", "reader")(
    ("FileNameList", [infilename] ),
)

# services to do Gulliver reconstruction
tray.AddService("I3SimpleParametrizationFactory","simpletrack")(
    ("StepX",20*I3Units.m),
    ("StepY",20*I3Units.m),
    ("StepZ",20*I3Units.m),
    ("StepZenith",0.1*I3Units.radian),
    ("StepAzimuth",0.2*I3Units.radian),
    ("BoundsX",[-2000*I3Units.m,2000*I3Units.m]),
    ("BoundsY",[-2000*I3Units.m,2000*I3Units.m]),
    ("BoundsZ",[-2000*I3Units.m,2000*I3Units.m]),
)

# minimizer (simplex, as implemented in ROOT's Minuit library)
tray.AddService("I3GulliverMinuitFactory","minuit")(
    ("Algorithm","SIMPLEX"),
)

# Gauss convoluted Pandel: standard workhorse for AMANDA/IceCube fitting
tray.AddService("I3GulliverIPDFPandelFactory","gpandel")(
    ("InputReadout", "FEPulses"),
    ("Likelihood", "SPE1st"),
    ("PEProb", "GaussConvoluted"),
    ("IceModel", 2),
    ("AbsorptionLength", 98.0*I3Units.m ),
    ("JitterTime", 4.0*I3Units.ns ),
)

# TFirst seems to work better than TMean, in IceCube
tray.AddService("I3BasicSeedServiceFactory","seedprep")(
    ("InputReadout", "FEPulses"),
    ("TimeShiftType", "TFirst"),
    ("FirstGuesses", ["IceDwalkFit"])
)

# when feeding a seed to paraboloid: don't change the time
tray.AddService("I3BasicSeedServiceFactory","paraseedprep")(
    ("InputReadout", "FEPulses"),
    ("TimeShiftType", "TNone"),
    ("FirstGuesses", ["pandelfit"])
)

##################################################################
# MODULES
##################################################################

tray.AddModule("I3Muxer","muxer")

tray.AddModule("I3FeatureExtractor","fe")(
    ("InitialHitSeriesReco","FEHits"),
    ("InitialPulseSeriesReco","FEPulses"),
    ("CalibratedFADCWaveforms","CalibratedFADC"),
    ("CalibratedATWDWaveforms","CalibratedATWD"),
    ("RawReadoutName","InIceRawData"),
    ("FastFirstPeak",3),
    ("FastPeakUnfolding",0),
    ("MaxNumHits", 0)            #set to 1-20 to enable slow root fit feature extraction
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
    ("inputpulsereadout","FEPulses"),
    ("name","IceDwalkFit"),
    ("multisolutions",True)
    )

# basic llh fit, nothing fancy
tray.AddModule("I3SimpleFitter","pandelfit")(
    ("SeedService","seedprep"),
    ("Parametrization","simpletrack"),
    ("LogLikelihood","gpandel"),
    ("Minimizer","minuit")
)

tray.AddModule("I3ParaboloidFitter","paraboloid")(
    ("SeedService","paraseedprep"),
    ("LogLikelihood","gpandel"),
    ("VertexStepSize",5.0*I3Units.m),
    ("MaxMissingGridPoints",1),
    ("GridpointVertexCorrection","seedprep"),
    ("Minimizer","minuit"),
)

# The IcePack converter
tray.AddModule("IcePacker","icepacker")(
    ("RootFileName","TestingIcePacker.root"),
    ("TrackExtract",True),
    ("DOMmapInIceName","InIceRawData"),
    ("DOMmapIceTopName","IceTopRawData"),
    ("WaveformmapInIceATWDName","CalibratedATWD"),
    ("WaveformmapInIceFADCName","CalibratedFADC"),
    ("WaveformmapIceTopATWDName","IceTopCalibratedATWD"),
    ("WaveformmapIceTopFADCName","noFADC"),
    ("InIcePulseSeriesName","FEPulses"),
    ("IceTopPulseSeriesName","IceTopFEPulses"),
#    ("InIcePulseSeriesNames",["FEPulses","SLCPulseSeriesReco"]), # If this parameter is set, the above single pulse-extraction is overruled
#    ("IceTopPulseSeriesNames",["IceTopFEPulses","IceTopSLCPulseSeries"]),# If this parameter is set, the above single pulse-extraction is overruled
    ("MCTreeName","I3MCTree"),
    ("MCWeightInfo","I3MCWeightDict"),
    ("EventHeaderName","I3EventHeader"),
    ("TriggerName","I3TriggerHierarchy"),
    ("setmaxtreesize",200000000)
    )

# write everything to an i3 file
tray.AddModule("I3Writer","i3out")(
    ("Filename","test.i3"),
)

# count events (and stop when we've enough)
tray.AddModule("I3EventCounter","nprocessed")(
    ("NEvents",10),
    ("CounterStep",10),
)

# from bits thou art made, to bits thou shalt return
tray.AddModule("TrashCan","destination")
    
# DO IT
tray.Execute()
tray.Finish()
