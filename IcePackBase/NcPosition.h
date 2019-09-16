#ifndef NcPosition_h
#define NcPosition_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include <stdlib.h>
#include <math.h>
 
#include "TObject.h"
#include "TString.h"

#include "Nc3Vector.h"
#include "NcTimestamp.h"
 
class NcPosition : public Nc3Vector
{
 public:
  NcPosition();                                         // Default constructor
  virtual ~NcPosition();                                // Destructor
  NcPosition(const NcPosition& p);                      // Copy constructor
  void SetPosition(Double_t* r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void SetPosition(Double_t r1,Double_t r2,Double_t r3,TString f,TString u="rad"); // Store position r
  void GetPosition(Double_t* r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  void SetPosition(Float_t*  r,TString f,TString u="rad");       // Store position r in frame f with ang units u
  void GetPosition(Float_t*  r,TString f,TString u="rad",Float_t s=-1) const; // Provide position r in frame f in ang units u
  NcPosition& GetPosition();                            // Provide position
  void SetPosition(Nc3Vector& r);                       // Store position r
  Double_t GetDistance(NcPosition& p,Float_t scale=-1); // Provide distance to position p
  Double_t GetDistance(NcPosition* p,Float_t scale=-1) { return GetDistance(*p,scale); }
  void SetPositionErrors(Double_t* e,TString f,TString u="rad"); // Store position errors in frame f with ang units u
  void SetPositionErrors(Double_t e1,Double_t e2,Double_t e3,TString f,TString u="rad"); // Store errors of position r
  void GetPositionErrors(Double_t* e,TString f,TString u="rad",Float_t s=-1) const; // Provide position errors in frame f in ang units u
  void SetPositionErrors(Float_t*  e,TString f,TString u="rad");       // Store position errors in frame f with ang units u
  void GetPositionErrors(Float_t*  e,TString f,TString u="rad",Float_t s=-1) const; // Provide position errors in frame f in ang units u
  void ResetPosition();                                  // Reset position and errors
  void SetUnitScale(Float_t s);                          // Set metric unit scale for the position coordinates
  Float_t GetUnitScale() const;                          // Provide metric unit scale for the position coordinates
  void SetTimestamp(NcTimestamp& t);                     // Set the timestamp for this position
  NcTimestamp* GetTimestamp();                           // Provide the timestamp for this position
  void RemoveTimestamp();                                // Remove the timestamp from this position
  virtual void Data(TString f="car",TString u="rad") const; // Print position/time info for frame f and ang units u

 protected:
  Float_t fScale;       // The unit scale used for the position coordinates
  NcTimestamp* fTstamp; // The timestamp for this position

 ClassDef(NcPosition,1) // Handling of positions (with timestamps) in various reference frames.
};
#endif
