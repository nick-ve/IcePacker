#ifndef IceDCDOM_h
#define IceDCDOM_h

// Copyright(c) 2003, IceCube Experiment at the South Pole, All rights reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "IceIDOM.h"

class IceDCDOM : public IceIDOM
{
 public:
  IceDCDOM();                                         // Default constructor
  virtual ~IceDCDOM();                                // Default destructor
  IceDCDOM(const IceDCDOM& m);                        // Copy constructor
  virtual TObject* Clone(const char* name="") const;  // Make a deep copy and provide its pointer

 ClassDef(IceDCDOM,1) // Signal (Hit) handling of an In-Ice DeepCore Digital Optical Module (DCDOM).
};
#endif
