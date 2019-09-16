#ifndef NcPositionObj_h
#define NcPositionObj_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"

#include "NcPosition.h"
 
class NcPositionObj : public TObject,public NcPosition
{
 public:
  NcPositionObj();                        // Default constructor
  NcPositionObj(NcPosition& p);           // Constructor
  virtual ~NcPositionObj();               // Destructor
  NcPositionObj(const NcPositionObj& p);  // Copy constructor

 ClassDef(NcPositionObj,1) // Handling of positions in various reference frames.
};
#endif
