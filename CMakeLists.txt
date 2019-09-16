#
# $Id$
#
i3_project(IcePacker)

i3_add_library(IcePacker
	IcePackBase/*.cxx
	private/IcePacker/*.cxx 
	USE_TOOLS python boost root #log4cplus :Uncomment for backwards compatibility
	USE_PROJECTS icetray serialization dataclasses recclasses phys-services dataio gulliver paraboloid linefit cramer-rao finiteReco CommonVariables
	ROOTCINT IcePackBase/IPBHeaders.h
)
