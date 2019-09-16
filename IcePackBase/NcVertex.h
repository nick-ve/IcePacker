#ifndef NcVertex_h
#define NcVertex_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

// $Id$

#include <stdlib.h>
#include <math.h>
 
#include "TObject.h"
#include "TObjArray.h"
#include "TPolyLine3D.h"
 
#include "NcJet.h"
#include "NcPosition.h"
 
class NcVertex : public NcJet,public NcPosition
{
 public:
  NcVertex();                              // Default constructor
  NcVertex(Int_t n);                       // Create a vertex to hold initially n tracks
  virtual ~NcVertex();                     // Default destructor
  NcVertex(const NcVertex& v);             // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  virtual void SetOwner(Bool_t own=kTRUE); // Set ownership of all added objects
  virtual void Reset();                    // Reset all values and stored vertex and jet lists
  void ResetVertices();                    // Reset stored vertex list
  void AddJet(NcJet& j,Int_t tracks=1);    // Add a jet (and its tracks) to the vertex
  void AddVertex(NcVertex& v,Int_t connect=1); // Add (and connect) a (sec.) vertex to the current vertex
  void AddJet(NcJet* j,Int_t tracks=1)    { AddJet(*j,tracks); }
  void AddVertex(NcVertex* v,Int_t connect=1) { AddVertex(*v,connect); }
  using NcPosition::Data;
  virtual void Data(TString f="car",TString u="rad"); // Print the vertex info within frame f and ang units u
  virtual void List(TString f="car",TString u="rad"); // Vertex prim. track info for frame f and ang units u
  virtual void ListAll(TString f="car",TString u="rad");// Prim.+sec. vtx full track info for frame f and ang units u
  Int_t GetNvertices() const;              // Return the number of (secondary) vertices
  NcVertex* GetVertex(Int_t i) const;      // Provide i-th (secondary) vertex
  NcVertex* GetIdVertex(Int_t id) const;   // Provide the vertex with user identifier "id"
  void SetNvmax(Int_t n=2);                // Set the initial max. number of (secondary) vertices
  void SetVertexCopy(Int_t j);             // (De)activate creation of private copies in fVertices
  Int_t GetVertexCopy() const;             // Provide VertexCopy flag value      
  Int_t GetNjets() const;                  // Return the number of jets
  NcJet* GetJet(Int_t i) const;            // Provide i-th jet
  NcJet* GetIdJet(Int_t id) const;         // Provide the jet with user identifier "id"
  void SetNjmax(Int_t n=2);                // Set the initial max. number of jets
  void SetJetCopy(Int_t j);                // (De)activate creation of private copies in fJets
  Int_t GetJetCopy() const;                // Provide JetCopy flag value      
  Int_t IsConnectTrack(NcTrack* t) const;  // Indicate if track is created by vertex connection
  Int_t IsJetTrack(NcTrack* t) const;      // Indicate if track is resulting from jet addition
  virtual void Draw(Option_t*) { Draw(1,1,0); } // Override TObject::Draw for default event display
  virtual void Draw(Int_t secs,Int_t cons=1,Int_t jets=0); // Draw the vertex in an event display
  TObjArray* SortJets(Int_t mode=-1,TObjArray* jets=0); // Sort jets by a certain observable

 protected:
  void Init();          // Initialisation of pointers etc... 
  Int_t fNvmax;         // The maximum number of (secondary) vertices
  Int_t fNvtx;          // The number of (secondary) vertices
  TObjArray* fVertices; // Array to hold the pointers to the (secondary) vertices
  TObjArray* fConnects; // Array to hold the pointers to the auto-generated connecting tracks
  Int_t fVertexCopy;    // Flag to denote creation of private copies in fVertices
  Int_t fNjmax;         // The maximum number of jets
  Int_t fNjets;         // The number of jets
  TObjArray* fJets;     // Array to hold the pointers to the jets
  TObjArray* fJetTracks;// Array to hold the pointers to tracks introduced by jet addition
  Int_t fJetCopy;       // Flag to denote creation of private copies in fJets
  TObjArray* fLines;    //! Array to (temporarily) store the 3D lines for the event display 

 private:
  void Dumps(NcVertex* v,Int_t n,TString f,TString u); // Recursively print all sec. vertices
 
 ClassDef(NcVertex,1) // Creation and investigation of an NcVertex.
};
#endif
