#ifndef NcTrack_h
#define NcTrack_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayD.h"
#include "TArrayI.h"
 
#include "NcDevice.h"
#include "NcBoost.h"
#include "NcPositionObj.h"
#include "NcTimestamp.h"
 
class NcTrack : public TNamed,public Nc4Vector
{
 public:
  NcTrack();                              // Default constructor
  virtual ~NcTrack();                     // Destructor
  NcTrack(const NcTrack& t);              // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  virtual void Reset();                   // Reset all values to 0
  void Set4Momentum(Nc4Vector& p);        // Set track 4-momentum
  void Set3Momentum(Nc3Vector& p);        // Set track 3-momentum
  void SetMass(Double_t m,Double_t dm=0); // Set particle mass and error
  void SetMass();                         // Set mass and error to the values of the hyp. with highest prob.
  void SetCharge(Float_t q);              // Set particle charge
  virtual void Data(TString f="car",TString u="rad"); // Print track information for frame f and ang units u
  virtual void List(TString f="car",TString u="rad"); // Track and decay level 1 info for frame f and ang units u
  virtual void ListAll(TString f="car",TString u="rad");// Track and all decay level info for frame f and ang units u
  Nc3Vector Get3Momentum(Float_t scale=-1) const; // Provide track 3-momentum
  Double_t GetMomentum(Float_t scale=-1); // Provide value of track 3-momentum
  Double_t GetMass(Float_t scale=-1);     // Provide particle mass
  Float_t GetCharge() const;              // Provide particle charge
  Double_t GetEnergy(Float_t scale=-1);   // Provide particle total energy
  void Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms); // Perform 2-body decay
  Int_t GetNdecay() const;                // Provide number of decay products
  NcTrack* GetDecayTrack(Int_t j) const;  // Access to decay produced track number j
  void RemoveDecays();                    // Remove all the decay products of this track
  void AddSignal(NcSignal& s,Int_t mode=0);     // Relate an NcSignal to this track
  void RemoveSignal(NcSignal& s,Int_t mode=1);  // Remove related NcSignal from this track
  void RemoveSignals(Int_t mode=1);             // Remove all related NcSignals from this track
  Int_t GetNsignals() const;              // Provide number of related NcSignals
  Int_t GetNsignals(const char* classname,Int_t par=0) const; // Provide the number of selected signals
  NcSignal* GetSignal(Int_t j) const;    // Access to the related NcSignal number j
  void SetBeginPoint(NcPosition& p);     // Set the track begin-point
  NcPosition* GetBeginPoint();           // Provide the track begin-point
  void SetEndPoint(NcPosition& p);       // Set the track end-point
  NcPosition* GetEndPoint();             // Provide the track end-point
  void SetReferencePoint(NcPosition& p); // Set the track reference-point for the 3-momentum vector
  NcPosition* GetReferencePoint();       // Provide the track reference-point for the 3-momentum vector
  void AddTrackHypothesis(NcTrack& t);   // Add track hypothesis
  void AddTrackHypothesis(Double_t prob,Double_t m,Double_t dm=0); // Add track hypothesis with mass data
  Int_t GetNhypotheses() const;           // Provide number of track hypotheses
  NcTrack* GetTrackHypothesis(Int_t j=0) const; // Provide the j-th track hypothesis 
  void RemoveTrackHypothesis(NcTrack& t);// Remove the specified track hypothesis 
  void RemoveTrackHypotheses();           // Remove all track hypotheses 
  Double_t GetPt(Float_t scale=-1);       // Provide trans. momentum w.r.t. z-axis
  Double_t GetPl(Float_t scale=-1);       // Provide long. momentum w.r.t. z-axis
  Double_t GetEt(Float_t scale=-1);       // Provide trans. energy w.r.t. z-axis
  Double_t GetEl(Float_t scale=-1);       // Provide long. energy w.r.t. z-axis
  Double_t GetMt(Float_t scale=-1);       // Provide trans. mass w.r.t. z-axis
  Double_t GetRapidity();                 // Provide rapidity value w.r.t. z-axis
  void SetImpactPoint(NcPosition& p,TString q); // Set the impact-point in plane "q=0"
  NcPosition* GetImpactPoint(TString q);        // Provide the impact-point in plane "q=0"
  void SetId(Int_t id);                   // Set the user defined unique track identifier
  Int_t GetId() const;                    // Provide the user defined unique track identifier
  void SetClosestPoint(NcPosition& p);    // Set position p as point of closest approach w.r.t. some reference
  NcPosition* GetClosestPoint();          // Provide point of closest approach w.r.t. some reference
  void SetParticleCode(Int_t code);       // Set the user defined particle id code (e.g. the PDF convention)
  Int_t GetParticleCode() const;          // Provide the user defined particle id code
  void SetParentTrack(NcTrack* t);        // Set pointer to the parent track
  NcTrack* GetParentTrack();              // Provide pointer to the parent track
  void SetProb(Double_t prob);            // Set the hypothesis probability for this track
  Float_t GetProb() const;                // Provide the hypothesis probability for this track
  void SetFitDetails(TObject* obj);       // Enter the object containing the fit details
  void SetFitDetails(TObject& obj) { SetFitDetails(&obj); }
  TObject* GetFitDetails();               // Provide pointer to the object containing the fit details
  void SetTimestamp(NcTimestamp& t);      // Set the track timestamp
  NcTimestamp* GetTimestamp();            // Provide the track timestamp
  void RemoveTimestamp();                 // Remove timestamp from this track
  Double_t GetDistance(NcPosition* p,Float_t scale=-1);   // Provide distance to position p
  Double_t GetDistance(NcPosition& p,Float_t scale=-1) { return GetDistance(&p,scale); }
  Double_t GetDistance(NcTrack* t,Float_t scale=-1);      // Provide distance to track t
  Double_t GetDistance(NcTrack& t,Float_t scale=-1) { return GetDistance(&t,scale); }
  void SetEscale(Float_t scale);          // Set the scale of the energy-momentum units of the track
  Float_t GetEscale() const;              // Provide the scale of the energy-momentum units of the track
 
 protected:
  void Init();               // Initialisation of pointers etc...
  Float_t fQ;                // The charge of the particle
  TObjArray* fDecays;        // The array of decay produced tracks
  TObjArray* fSignals;       // The array of related NcSignals
  TObjArray* fHypotheses;    // The array of track hypotheses
  NcPositionObj* fBegin;     // The begin-point of the track 
  NcPositionObj* fEnd;       // The end-point of the track 
  NcPositionObj* fRef;       // The reference-point of the track for the 3-momentum vector
  NcPositionObj* fImpactXY;  // The (extrapolated) impact-point in the plane z=0
  NcPositionObj* fImpactXZ;  // The (extrapolated) impact-point in the plane y=0
  NcPositionObj* fImpactYZ;  // The (extrapolated) impact-point in the plane x=0
  Int_t fUserId;             // The user defined identifier
  NcPositionObj* fClosest;   // The (extrapolated) point of closest approach w.r.t some reference
  Int_t fCode;               // The user defined particle id code
  NcTrack* fParent;          // Pointer to the parent track
  Float_t fProb;             // Probability for this track as a hypothesis
  TObject* fFit;             // Object containing details of the fit
  NcTimestamp* fTstamp;      // The track timestamp
  Float_t fEscale;           // The scale of the energy-momentum units of the track

 private:
  void Dumps(NcTrack* t,Int_t n,TString f,TString u); // Recursively print all decay levels
 
 ClassDef(NcTrack,1) // Handling of the attributes of a reconstructed particle track.
};
#endif
