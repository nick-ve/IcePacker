#ifndef Nc3Vector_h
#define Nc3Vector_h
// Copyright(c) 1997-2019, NCFS/IIHE, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id: Nc3Vector.h 5 2010-03-19 10:10:02Z nickve $

#include <math.h>
 
#include "TObject.h"
#include "TString.h"
#include "TRotMatrix.h"
 
class Nc3Vector
{
 public:
  Nc3Vector();                                  // Default constructor
  virtual ~Nc3Vector();                         // Destructor
  Nc3Vector(const Nc3Vector& v);                // Copy constructor
  Nc3Vector& operator=(const Nc3Vector& q);     // Set equal to vector q
  virtual void Load(Nc3Vector& q);              // Load all attributes of input Nc3Vector
  virtual void SetZero();                       // (Re)set all attributes to zero.
  void SetVector(Double_t* v,TString f,TString u="rad");       // Store vector v in frame f with ang units u
  void SetVector(Double_t v1,Double_t v2,Double_t v3,TString f,TString u="rad"); // Store vector v
  void GetVector(Double_t* v,TString f,TString u="rad") const; // Provide vector v in frame f in ang units u
  void SetVector(Float_t*  v,TString f,TString u="rad");       // Store vector v in frame f with ang units u
  void GetVector(Float_t*  v,TString f,TString u="rad") const; // Provide vector v in frame f in ang units u
  void SetErrors(Double_t* e,TString f,TString u="rad");       // Store errors of vector in frame f with ang units u
  void SetErrors(Double_t e1,Double_t e2,Double_t e3,TString f,TString u="rad"); // Store errors of vector v
  void GetErrors(Double_t* e,TString f,TString u="rad") const; // Provide errors of vector in frame f in ang units u
  void SetErrors(Float_t*  e,TString f,TString u="rad");       // Store errors of vector in frame f with ang units u
  void GetErrors(Float_t*  e,TString f,TString u="rad") const; // Provide errors of vector in frame f in ang units u
  virtual void Data(TString f="car",TString u="rad") const;    // Print vector components in frame f in ang units u
  Int_t HasVector() const;                      // Check whether vector components have been set
  Int_t HasErrors() const;                      // Check whether errors on the vector components have been set
  Double_t GetNorm();                           // Provide norm of the vector
  Double_t Dot(Nc3Vector& q);                   // Provide dot product with q
  Double_t GetPseudoRapidity();                 // Provide the pseudorapidity w.r.t z-axis
  Double_t GetResultError() const;              // Provide error on scalar result (e.g. norm)
  Nc3Vector Cross(Nc3Vector& q) const;          // Provide cross product with q
  Nc3Vector operator+(Nc3Vector& q) const;      // Add vector q
  Nc3Vector operator-(Nc3Vector& q) const;      // Subtract vector q
  Nc3Vector operator*(Double_t s) const;        // Multiply vector with scalar s
  Nc3Vector operator/(Double_t s) const;        // Divide vector by scalar s
  Nc3Vector& operator+=(Nc3Vector& q);          // Add vector q
  Nc3Vector& operator-=(Nc3Vector& q);          // Subtract vector q
  Nc3Vector& operator*=(Double_t s);            // Multiply with scalar s
  Nc3Vector& operator/=(Double_t s);            // Divide by scalar s
  Nc3Vector GetVecTrans() const;                // Provide transverse vector w.r.t. z-axis
  Nc3Vector GetVecLong() const;                 // Provide longitudinal vector w.r.t. z-axis
  Nc3Vector GetPrimed(TRotMatrix* m) const;     // Provide vector components in a rotated frame
  Nc3Vector GetUnprimed(TRotMatrix* m) const;   // Provide original vector components from a rotated one
  Double_t GetX(Int_t i,TString f,TString u="rad"); // Provide i-th vector component in frame f in units u
  virtual Double_t GetOpeningAngle(Nc3Vector& q,TString u="rad"); // Provide opening angle with q in units u

 protected:
  Int_t fNv;           // The number of components of the array fV
  Double32_t* fV;      //[fNv] Vector in spherical (and errors in Cartesian) coordinates
  Double32_t fDresult; //! Error on scalar result (e.g. norm or dotproduct)

 ClassDef(Nc3Vector,1) // Handling of 3-vectors in various reference frames.
};
#endif
