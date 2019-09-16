#ifndef Nc4VectorObj_h
#define Nc4VectorObj_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"

#include "Nc4Vector.h"
 
class Nc4VectorObj : public TObject,public Nc4Vector
{
 public:
  Nc4VectorObj();                       // Default constructor
  Nc4VectorObj(Nc4Vector& q);           // Constructor
  virtual ~Nc4VectorObj();              // Destructor
  Nc4VectorObj(Nc4VectorObj& q);        // Copy constructor

 ClassDef(Nc4VectorObj,1) // Handling of Lorentz 4-vectors in various reference frames.
};
#endif
