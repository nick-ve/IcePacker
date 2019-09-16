#ifndef Nc3VectorObj_h
#define Nc3VectorObj_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"

#include "Nc3Vector.h"
 
class Nc3VectorObj : public TObject,public Nc3Vector
{
 public:
  Nc3VectorObj();                       // Default constructor
  Nc3VectorObj(Nc3Vector& q);           // Constructor
  virtual ~Nc3VectorObj();              // Destructor
  Nc3VectorObj(Nc3VectorObj& q);        // Copy constructor

 ClassDef(Nc3VectorObj,1) // Handling of 3-vectors in various reference frames.
};
#endif
