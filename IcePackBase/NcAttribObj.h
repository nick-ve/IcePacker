#ifndef NcAttribObj_h
#define NcAttribObj_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"

#include "NcAttrib.h"

class NcAttribObj : public TObject,public NcAttrib
{
 public:
  NcAttribObj();                                   // Default constructor
  NcAttribObj(NcAttrib& a);                        // Constructor
  virtual ~NcAttribObj();                          // Destructor
  NcAttribObj(const NcAttribObj& a);               // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer

 ClassDef(NcAttribObj,1) // Generic handling of detector signal (calibration) attributes.
};
#endif
