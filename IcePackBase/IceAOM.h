#ifndef IceAOM_h
#define IceAOM_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id: IceAOM.h 3 2007-12-27 13:56:06Z nick $

#include "IceGOM.h"

class IceAOM : public IceGOM
{
 public:
  IceAOM();                                          // Default constructor
  virtual ~IceAOM();                                 // Default destructor
  IceAOM(const IceAOM& m);                           // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(IceAOM,1) // Signal (Hit) handling of a generic Amanda Optical Module (AOM).
};
#endif
