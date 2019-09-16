#ifndef NcBoost_h
#define NcBoost_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include <stdlib.h>
#include <math.h>
 
#include "TObject.h"

#include "Nc4Vector.h" 

class NcBoost : public TObject
{
 public:
  NcBoost();                             // Default constructor
  virtual ~NcBoost();                    // Default destructor
  NcBoost(const NcBoost& b);             // Copy constructor
  void SetBeta(Nc3Vector& b);            // Set boost parameters by beta 3-vector
  void Set4Momentum(Nc4Vector& p);       // Set boost parameters by 4-momentum
  Nc3Vector GetBetaVector() const;       // Provide the beta 3-vector
  Double_t GetBeta();                    // Provide norm of beta 3-vector
  Double_t GetGamma();                   // Provide gamma value
  void Data(TString f="car");            // Print boost parameter info in coord. frame f
  Nc4Vector Boost(Nc4Vector& v);         // Perform Lorentz boost on 4-vector v
  Nc4Vector Inverse(Nc4Vector& v);       // Perform inverse Lorentz boost on 4-vector v
  Double_t GetResultError() const;       // Provide error on scalar result
 
 protected:
  Nc3Vector fBeta;     // The beta 3-vector
  Double32_t fGamma;   // The gamma factor
  Double32_t fDgamma;  // Error on the gamma value
  Double32_t fDresult; //! Error on scalar result
 
 ClassDef(NcBoost,1) // Perform various Lorentz transformations.
};
#endif
