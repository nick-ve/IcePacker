/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright(c) 1997-2009, NCFS, All Rights Reserved.                          *
 *                                                                             *
 * Author: The Netherlands Center for Fundamental Studies (NCFS).              *
 *         http://sites.google.com/site/ncfsweb ncfs.nl@gmail.com              *
 *                                                                             *
 * Contributors are mentioned in the code where appropriate.                   *
 *                                                                             * 
 * No part of this software may be used, copied, modified or distributed       *
 * by any means nor transmitted or translated into machine language without    *
 * written permission by the NCFS.                                             *
 * Permission to use the documentation strictly for non-commercial purposes    *
 * is hereby granted without fee, provided that the above copyright notice     *
 * appears in all copies and that both the copyright notice and this           *
 * permission notice appear in the supporting documentation.                   *
 * This software is provided "as is" without express or implied warranty.      *
 * The authors make no claims that this software is free of error, or is       *
 * consistent with any particular standard of merchantability, or that it      *
 * will meet your requirements for any particular application, other than      *
 * indicated in the corresponding documentation.                               *
 * This software should not be relied on for solving a problem whose           *
 * incorrect solution could result in injury to a person or loss of property.  *
 * If you do use this software in such a manner, it is at your own risk.       *
 * The authors disclaim all liability for direct or consequential damage       *
 * resulting from your use of this software.                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class NcVertex
// Creation and investigation of an NcVertex.
// An NcVertex can be constructed by adding NcTracks and/or NcJets.
//
// Note : Also (secondary) vertices can be added to a vertex.
//
// To provide maximal flexibility to the user, two modes of vertex storage
// are provided by means of the memberfunction SetVertexCopy().
// The same holds for the storage of jets via SetJetCopy().
//
// a) SetVertexCopy(0) (which is the default).
//    Only the pointers of the 'added' vertices are stored.
//    This mode is typically used by making vertex studies based on a fixed list
//    of vertices which stays under user control or is contained for instance
//    in an NcEvent.  
//    In this way the NcVertex just represents a 'logical structure' for the
//    physics analysis which can be embedded in e.g. an NcEvent or NcVertex.
//
//    Note :
//    Modifications made to the original vertices also affect the NcVertex objects
//    which are stored.
// 
// b) SetVertexCopy(1).
//    Of every 'added' vertex a private copy will be made of which the pointer
//    will be stored.
//    In this way the NcVertex represents an entity on its own and modifications
//    made to the original vertices do not affect the NcVertex objects which are
//    stored. 
//    This mode will allow 'adding' many different NcVertex objects by
//    creating only one NcVertex instance in the main programme and using the
//    NcVertex::Reset, NcVertex::AddTrack and parameter setting memberfunctions.
//
// See also the documentation provided for the memberfunction SetOwner(). 
//
// Coding example to make 3 vertices v1, v2 and v3.
// ------------------------------------------------
// v1 contains the tracks 1,2,3 and 4
// v2 contains many different tracks
// v3 contains the jets 1 and 2
//
//        NcTrack t1,t2,t3,t4;
//         ...
//         ... // code to fill the track data
//         ...
//
//        NcJet j1,j2;
//         ...
//         ... // code to fill the jet data
//         ...
//
//        NcVertex v1;
//        v1.SetVertexCopy(1);
//
//        v1.AddTrack(t1);
//        v1.AddTrack(t2);
//        v1.AddTrack(t3);
//        v1.AddTrack(t4);
//
//        Float_t r1[3]={2.4,0.1,-8.5};
//        v1.SetPosition(r1,"car");
//
//        NcVertex v2;
//        v2.SetTrackCopy(1);
//
//        NcTrack* tx=new NcTrack();
//        for (Int_t i=0; i<10; i++)
//        {
//         ...
//         ... // code to fill the track data
//         ...
//         v2.AddTrack(tx);
//         tx->Reset(); 
//        }
//
//        Float_t r2[3]={1.6,-3.2,5.7};
//        v2.SetPosition(r2,"car");
//
//        NcVertex v3;
//
//        v3.AddJet(j1);
//        v3.AddJet(j2);
//
//        Float_t r3[3]={6.2,4.8,1.3};
//        v3.SetPosition(r3,"car");
//
//        v1.Data("sph");
//        v2.ListAll();
//        v3.List("cyl");
//
//        Float_t e1=v1.GetEnergy();
//        Nc3Vector p1=v1.Get3Momentum();
//        Float_t loc[3];
//        v1.GetPosition(loc,"sph");
//        NcPosition r=v2.GetPosition();
//        r.Data(); 
//        Int_t nt=v2.GetNtracks();
//        NcTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
// Specify the vertices v2 and v3 as secondary vertices of v1
//
//        v1.AddVertex(v2);
//        v1.AddVertex(v3);
//
//        v1.List();
//
//        Int_t nv=v1.GetNvtx();
//        NcVertex* vx=v1.GetVertex(1); // Access 1st secondary vertex of v1
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=v1.GetInvmass(); 
//
// Reconstruct Vertex v1 from scratch
//
//        v1.Reset();
//        v1.SetNvmax(25); // Increase initial no. of sec. vertices
//        v1.AddTrack(t3);
//        v1.AddTrack(t4);
//        v1.AddJet(j2);
//        Float_t pos[3]={7,9,4};
//        v1.SetPosition(pos,"car");
//
// Note : By default all quantities are in meter, GeV, GeV/c or GeV/c**2
//        but the user can indicate the usage of a different scale for
//        the metric and/or energy-momentum units via the SetUnitScale()
//        and SetEscale() memberfunctions, respectively.
//        The actual metric and energy-momentum unit scales in use can be
//        obtained via the GetUnitScale() and GetEscale() memberfunctions.
//
//--- Author: Nick van Eijndhoven 04-apr-1998 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcVertex.h"
#include "Riostream.h"
 
ClassImp(NcVertex) // Class implementation to enable ROOT I/O
 
NcVertex::NcVertex() : NcJet(),NcPosition()
{
// Default constructor.
// All variables initialised to 0.
// Initial maximum number of tracks is set to the default value.
// Initial maximum number of sec. vertices is set to the default value.
 Init();
 Reset();
 SetNvmax();
 SetNjmax();
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::Init()
{
// Initialisation of pointers etc...
 fNvmax=0;
 fVertices=0;
 fConnects=0;
 fVertexCopy=0;
 fNjmax=0;
 fJets=0;
 fJetTracks=0;
 fJetCopy=0;
 fLines=0;
}
///////////////////////////////////////////////////////////////////////////
NcVertex::NcVertex(Int_t n) : NcJet(n),NcPosition()
{
// Create a vertex to hold initially a maximum of n tracks
// All variables initialised to 0
 if (n<=0)
 {
  cout << " *** This NcJet initialisation was invoked via the NcVertex ctor." << endl;
 }
 Init();
 Reset();
 SetNvmax();
 SetNjmax();
}
///////////////////////////////////////////////////////////////////////////
NcVertex::~NcVertex()
{
// Default destructor
 if (fVertices)
 {
  delete fVertices;
  fVertices=0;
 }
 if (fConnects)
 {
  delete fConnects;
  fConnects=0;
 }
 if (fJets)
 {
  delete fJets;
  fJets=0;
 }
 if (fJetTracks)
 {
  delete fJetTracks;
  fJetTracks=0;
 }
 if (fLines)
 {
  delete fLines;
  fLines=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::SetOwner(Bool_t own)
{
// Set ownership of all added objects. 
// The default parameter is own=kTRUE.
//
// Invokation of this memberfunction also sets all the copy modes
// (e.g. TrackCopy & co.) according to the value of own.
//
// This function (with own=kTRUE) is particularly useful when reading data
// from a tree/file, since Reset() will then actually remove all the
// added objects from memory irrespective of the copy mode settings
// during the tree/file creation process. In this way it provides a nice way
// of preventing possible memory leaks in the reading/analysis process.
//
// In addition this memberfunction can also be used as a shortcut to set all
// copy modes in one go during a tree/file creation process.
// However, in this case the user has to take care to only set/change the
// ownership (and copy mode) for empty objects (e.g. newly created objects
// or after invokation of the Reset() memberfunction) otherwise it will
// very likely result in inconsistent destructor behaviour.

 Int_t mode=1;
 if (!own) mode=0;
 if (fVertices) fVertices->SetOwner(own);
 fVertexCopy=mode;
 if (fJets) fJets->SetOwner(own);
 fJetCopy=mode;

 NcJet::SetOwner(own);
}
///////////////////////////////////////////////////////////////////////////
NcVertex::NcVertex(const NcVertex& v) : NcJet(v.fNtinit),NcPosition(v)
{
// Copy constructor
 Init();
 fNvtx=0;
 fNjets=0;
 SetNvmax(v.fNvmax);
 SetNjmax(v.fNjmax);
 SetTrackCopy(v.GetTrackCopy());
 SetVertexCopy(v.GetVertexCopy());
 SetJetCopy(v.GetJetCopy());
 SetId(v.GetId());

 // Copy all tracks except the ones coming from jets
 NcTrack* tx=0;
 Int_t jetflag=0,connect=0;
 NcTrack* tx2=0;
 for (Int_t it=1; it<=v.fNtrk; it++)
 {
  tx=v.GetTrack(it);
  if (tx)
  {
   jetflag=v.IsJetTrack(tx);
   connect=v.IsConnectTrack(tx);

   if (!jetflag && !connect) AddTrack(tx);

   if (connect)
   {
    if (!fConnects)
    {
     fConnects=new TObjArray(fNvmax);
     if (!fTrackCopy) fConnects->SetOwner();
    }
    tx2=new NcTrack(*tx);
    fConnects->Add(tx2);
    AddTrack(tx2,0);
   } 
  }
 }

 // Copy all the (secondary) vertices without re-creating connecting tracks
 // The connecting tracks have already been copied above
 NcVertex* vx=0;
 for (Int_t iv=1; iv<=v.GetNvertices(); iv++)
 {
  vx=v.GetVertex(iv);
  if (vx) AddVertex(vx,0); 
 }

 // Copy all the jets including the jet tracks for these jets for which
 // this was also the case in the original vertex
 NcJet* jx=0;
 for (Int_t ij=1; ij<=v.GetNjets(); ij++)
 {
  jx=v.GetJet(ij);
  if (jx)
  {
   jetflag=0;
   if (jx->GetNtracks())
   {
    tx=jx->GetTrack(1);
    if (tx)
    {
     jetflag=v.IsJetTrack(tx);
    }
   }
   AddJet(jx,jetflag);
  } 
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::SetNvmax(Int_t n)
{
// Set the initial maximum number of (secondary) vertices
 if (n > 0)
 {
  fNvmax=n;
 }
 else
 {
  fNvmax=1;
 }
 if (fVertices)
 {
  delete fVertices;
  fVertices=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::SetNjmax(Int_t n)
{
// Set the initial maximum number of jets
 if (n > 0)
 {
  fNjmax=n;
 }
 else
 {
  fNjmax=1;
 }
 if (fJets)
 {
  delete fJets;
  fJets=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::Reset()
{
// Reset all variables to 0 and reset all stored vertex and jet lists.
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
// The max. number of jets is set to the default value again

 NcJet::Reset();

 ResetPosition();

 fNvtx=0;
 if (fNvmax>0) SetNvmax(fNvmax);
 if (fConnects)
 {
  delete fConnects;
  fConnects=0;
 }

 fNjets=0;
 if (fNjmax>0) SetNjmax(fNjmax);
 if (fJetTracks)
 {
  delete fJetTracks;
  fJetTracks=0;
 }
 
 if (fLines)
 {
  delete fLines;
  fLines=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::ResetVertices()
{
// Reset the stored vertex list and delete all connecting tracks which
// were generated automatically via connect=1 in AddVertex().
// The max. number of vertices is set to the default value again.
// All physics quantities are updated according to the removal of the
// connecting tracks.
 NcTrack* t;
 if (fConnects)
 {
  for (Int_t i=0; i<=fConnects->GetLast(); i++)
  {
   t=(NcTrack*)fConnects->At(i);
   NcTrack* test=(NcTrack*)fTracks->Remove(t);
   if (test)
   {
    fNtrk--;
    (Nc4Vector&)(*this)-=(Nc4Vector&)(*t);
    fQ-=t->GetCharge();
    if (fTrackCopy) delete t;
   }
  }
  fTracks->Compress();
 }

 fNvtx=0;
 if (fNvmax>0) SetNvmax(fNvmax);
 if (fConnects)
 {
  delete fConnects;
  fConnects=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::AddJet(NcJet& j,Int_t tracks)
{
// Add a jet (and its tracks) to the vertex
// In case the maximum number of jets has been reached,
// the array space will be extended automatically
//
// Note : By default the tracks of the jet are added to the current (primary)
//        vertex.
//        The automatic addition of the tracks of the jet can be suppressed
//        by specifying tracks=0. In this case only the NcJet object will
//        be stored according to the mode specified by SetJetCopy().
//        The latter will enable jet studies based on a fixed list of tracks
//        as contained e.g. in an NcVertex or NcEvent. 
//
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All NcJet and derived classes have the default TObject::Clone() memberfunction.
// However, derived classes generally contain an internal data structure which may
// include pointers to other objects. Therefore it is recommended to provide
// for all derived classes a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from NcJet.   

 if (!fJets)
 {
  fJets=new TObjArray(fNjmax);
  if (fJetCopy) fJets->SetOwner();
 }
 if (fNjets == fNjmax) // Check if maximum jet number is reached
 {
  fNjmax++;
  fJets->Expand(fNjmax);
 }

 // Add the jet to the list 
 NcJet* jx=&j;
 if (fJetCopy) jx=(NcJet*)j.Clone();

 if (jx)
 {
  fNjets++;
  fJets->Add(jx); 
 }

 // Add the tracks of the jet to this vertex
 if (tracks)
 {
  if (!fJetTracks)
  {
   fJetTracks=new TObjArray();
  }
  Int_t copy=1-(jx->GetTrackCopy());
  NcTrack* tj;
  for (Int_t i=1; i<=jx->GetNtracks(); i++)
  {
   tj=jx->GetTrack(i);
   if (tj)
   {
    AddTrack(tj,copy);
    fJetTracks->Add(tj);
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::AddVertex(NcVertex& v,Int_t connect)
{
// Add a (secondary) vertex to the current vertex.
// In case the maximum number of (secondary) vertices has been reached,
// the array space will be extended automatically
//
// Note : By default the 4-momentum and charge of the current (primary) vertex
//        are updated by automatically creating the track connecting
//        both vertices. The track parameters are taken from the
//        4-momentum and charge of the secondary vertex.
//        The automatic creation of the connecting track and updating
//        of the (primary) vertex 4-momentum and charge can be suppressed
//        by specifying connect=0. In this case, however, the user
//        has to introduce the connecting track lateron by hand
//        explicitly in order to match the kinematics and charge.
//
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All NcVertex and derived classes have the default TObject::Clone() memberfunction.
// However, derived classes generally contain an internal data structure which may
// include pointers to other objects. Therefore it is recommended to provide
// for all derived classes a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from NcVertex.   

 if (!fVertices)
 {
  fVertices=new TObjArray(fNvmax);
  if (fVertexCopy) fVertices->SetOwner();
 }
 if (fNvtx == fNvmax) // Check if maximum vertex number is reached
 {
  fNvmax++;
  fVertices->Expand(fNvmax);
 }

 // Add the linked (secondary) vertex to the list 
 NcVertex* vx=&v;
 if (fVertexCopy) vx=(NcVertex*)v.Clone();

 if (vx)
 {
  fNvtx++;
  fVertices->Add(vx);
 } 

 // Create connecting track and update 4-momentum and charge for current vertex
 if (connect)
 {
  NcTrack* t=new NcTrack();
  t->SetBeginPoint(GetPosition());
  t->SetEndPoint(v.GetPosition());
  t->SetCharge(v.GetCharge());
  t->Set4Momentum((Nc4Vector&)v);

  AddTrack(t,0);

  if (!fConnects)
  {
   fConnects=new TObjArray(fNvmax);
   if (!fTrackCopy) fConnects->SetOwner();
  }
  fConnects->Add(t);
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::Data(TString f,TString u)
{
// Provide vertex information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 const char* name=GetName();
 const char* title=GetTitle();
 cout << " *NcVertex::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;
 cout << " Id : " << fUserId << " Invmass : " << GetInvmass()
      << " Charge : " << GetCharge() << " Momentum : " << GetMomentum()
      << " Ntracks : " << GetNtracks() << endl;
 cout << " Nvertices : " << fNvtx << " Njets : " << fNjets
      << " Energy scale : " << fEscale << " GeV" << endl;
 cout << " ";
 Nc4Vector::Data(f,u);
 cout << "  Position";
 NcPosition::Data(f,u); 
} 
///////////////////////////////////////////////////////////////////////////
void NcVertex::List(TString f,TString u)
{
// Provide primary track and sec. vertex information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current vertex

 // The tracks of this vertex
 NcTrack* t; 
 for (Int_t it=1; it<=GetNtracks(); it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->Data(f,u); 
  }
  else
  {
   cout << " *NcVertex::List* Error : No track present." << endl; 
  }
 }

 // The secondary vertices of this vertex
 NcVertex* v; 
 for (Int_t iv=1; iv<=GetNvertices(); iv++)
 {
  v=GetVertex(iv);
  if (v)
  {
   cout << "  ---Level 1 sec. vertex no. " << iv << endl;
   cout << " ";
   v->Data(f,u); 
  }
  else
  {
   cout << " *NcVertex::List* Error : No sec. vertex present." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void NcVertex::ListAll(TString f,TString u)
{
// Provide complete (sec) vertex and (decay) track info within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current vertex

 // The tracks of this vertex
 NcTrack* t; 
 for (Int_t it=1; it<=GetNtracks(); it++)
 {
  t=GetTrack(it);
  if (t)
  {
   cout << "  ---Track no. " << it << endl;
   cout << " ";
   t->ListAll(f,u); 
  }
  else
  {
   cout << " *NcVertex::ListAll* Error : No track present." << endl; 
  }
 }

 NcVertex* v=this;
 Dumps(v,1,f,u); // Information of all sec. vertices
}
//////////////////////////////////////////////////////////////////////////
void NcVertex::Dumps(NcVertex* v,Int_t n,TString f,TString u)
{
// Recursively provide the info of all secondary vertices of this vertex
 NcVertex* vs; 
 for (Int_t iv=1; iv<=v->GetNvertices(); iv++)
 {
  vs=v->GetVertex(iv);
  if (vs)
  {
   cout << "  ---Level " << n << " sec. vertex no. " << iv << endl;
   cout << " ";
   vs->Data(f,u); 

   // The tracks of this vertex
   NcTrack* t; 
   for (Int_t it=1; it<=vs->GetNtracks(); it++)
   {
    t=vs->GetTrack(it);
    if (t)
    {
     cout << "  ---Track no. " << it << endl;
     cout << " ";
     t->ListAll(f,u); 
    }
    else
    {
     cout << " *NcVertex::Dumps* Error : No track present." << endl; 
    }
   }

   // Go for next sec. vertex level of this sec. vertex recursively
   Dumps(vs,n+1,f,u);
  }
  else
  {
   cout << " *NcVertex::Dumps* Error : No sec. vertex present." << endl; 
  }
 }
} 
//////////////////////////////////////////////////////////////////////////
Int_t NcVertex::GetNvertices() const
{
// Return the current number of (secondary) vertices
 return fNvtx;
}
///////////////////////////////////////////////////////////////////////////
NcVertex* NcVertex::GetVertex(Int_t i) const
{
// Return the i-th (secondary) vertex of the current vertex
 if (!fVertices)
 {
  cout << " *NcVertex*::GetVertex* No (secondary) vertices present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNvtx)
  {
   cout << " *NcVertex*::GetVertex* Invalid argument i : " << i
        << " Nvtx = " << fNvtx << endl;
   return 0;
  }
  else
  {
   return (NcVertex*)fVertices->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
NcVertex* NcVertex::GetIdVertex(Int_t id) const
{
// Return the (sec.) vertex with user identifier "id"
 NcVertex* vx=0;
 NcVertex* v=0;
 if (!fVertices)
 {
  cout << " *NcVertex*::GetIdVertex* No (secondary) vertices present." << endl;
  return 0;
 }
 else
 {
  for (Int_t i=0; i<fNvtx; i++)
  {
   vx=(NcVertex*)fVertices->At(i);
   if (id == vx->GetId()) v=vx;
  }
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::SetVertexCopy(Int_t j)
{
// (De)activate the creation of private copies of the added vertices.
// j=0 ==> No private copies are made; pointers of original vertices are stored.
// j=1 ==> Private copies of the vertices are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to NcVertex objects one cannot
//        change the VertexCopy mode anymore.
//        To change the VertexCopy mode for an existing NcVertex containing
//        vertices one first has to invoke Reset().
 if (!fVertices)
 {
  if (j==0 || j==1)
  {
   fVertexCopy=j;
  }
  else
  {
   cout << "*NcVertex::SetVertexCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*NcVertex::SetVertexCopy* Storage already contained vertices."
       << "  ==> VertexCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcVertex::GetVertexCopy() const
{
// Provide value of the VertexCopy mode.
// 0 ==> No private copies are made; pointers of original vertices are stored.
// 1 ==> Private copies of the vertices are made and these pointers are stored.
 return fVertexCopy;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcVertex::GetNjets() const
{
// Return the current number of jets
 return fNjets;
}
///////////////////////////////////////////////////////////////////////////
NcJet* NcVertex::GetJet(Int_t i) const
{
// Return the i-th jet of the current vertex
 if (!fJets)
 {
  cout << " *NcVertex*::GetJet* No jets present." << endl;
  return 0;
 }
 else
 {
  if (i<=0 || i>fNjets)
  {
   cout << " *NcVertex*::GetJet* Invalid argument i : " << i
        << " Njets = " << fNjets << endl;
   return 0;
  }
  else
  {
   return (NcJet*)fJets->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
NcJet* NcVertex::GetIdJet(Int_t id) const
{
// Return the jet with user identifier "id"
 NcJet* jx=0;
 NcJet* j=0;
 if (!fJets)
 {
  cout << " *NcVertex*::GetIdJet* No jets present." << endl;
  return 0;
 }
 else
 {
  for (Int_t i=0; i<fNjets; i++)
  {
   jx=(NcJet*)fJets->At(i);
   if (id == jx->GetId()) j=jx;
  }
  return j;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::SetJetCopy(Int_t j)
{
// (De)activate the creation of private copies of the added jets.
// j=0 ==> No private copies are made; pointers of original jets are stored.
// j=1 ==> Private copies of the jets are made and these pointers are stored.
//
// Note : Once the storage contains pointer(s) to NcJet objects one cannot
//        change the JetCopy mode anymore.
//        To change the JetCopy mode for an existing NcVertex containing
//        jets one first has to invoke Reset().
 if (!fJets)
 {
  if (j==0 || j==1)
  {
   fJetCopy=j;
  }
  else
  {
   cout << "*NcVertex::SetJetCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << "*NcVertex::SetJetCopy* Storage already contained jets."
       << "  ==> JetCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcVertex::GetJetCopy() const
{
// Provide value of the JetCopy mode.
// 0 ==> No private copies are made; pointers of original jets are stored.
// 1 ==> Private copies of the jets are made and these pointers are stored.
 return fJetCopy;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcVertex::IsConnectTrack(NcTrack* t) const
{
// Indicate whether a track from the tracklist was created via the
// connection of a (secondary) vertex or not.
// In case the track was the result of (secondary) vertex addition the
// return value is 1, otherwise the value 0 will be returned.
 Int_t connect=0;
 if (fConnects)
 {
  if (fConnects->FindObject(t)) connect=1;
 }
 return connect;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcVertex::IsJetTrack(NcTrack* t) const
{
// Indicate whether a track from the tracklist was created via the
// addition of a jet or not.
// In case the track was the result of jet addition the return value is 1,
// otherwise the value 0 will be returned.
 Int_t jetflag=0;
 if (fJetTracks)
 {
  if (fJetTracks->FindObject(t)) jetflag=1;
 }
 return jetflag;
}
///////////////////////////////////////////////////////////////////////////
void NcVertex::Draw(Int_t secs,Int_t cons,Int_t jets)
{
// 3-Dimensional visualisation of an NcVertex with its attributes.
// The displayed tracklength is proportional to the momentum of the track.
//
// Color conventions :
// -------------------
// positive track : red
// neutral  track : green
// negative track : blue
// jet-track      : magenta (if explicit marking selected)
//
// secs = 1 --> Draw secondary vertices. (Default)
//        0 --> Don't draw secondary vertices.
//
// cons = 1 --> Draw (auto generated) connecting tracks. (Default) 
//        0 --> Don't draw (auto generated) connecting tracks.
//                  
// jets = 1 --> Mark tracks belonging to jets.
//        0 --> Don't mark jet-tracks. (Default)
//
// Notes :
// -------
// Auto generated connecting tracks will be drawn as thin lines.
// Tracks belonging to jets will be marked as somewhat thinner magenta lines.
// This memberfunction is used recursively.
//
 Double_t vec[3]={0,0,0};
 NcTrack* tx=0;
 NcVertex* vx=0;
 NcPosition* r=0;
 Nc3Vector p;
 Float_t charge;

 NcPosition dummy;

 if (fLines) delete fLines;
 fLines=new TObjArray();
 fLines->SetOwner();

 Int_t ntk=GetNtracks();
 for (Int_t jtk=1; jtk<=ntk; jtk++)
 {
  tx=GetTrack(jtk);

  if (!tx) continue;

  charge=tx->GetCharge();

  TPolyLine3D* line=new TPolyLine3D();
  fLines->Add(line);

  if (IsConnectTrack(tx))
  {
   if (cons==1)
   {
    r=tx->GetBeginPoint();
    if (!r) r=&dummy;
    r->GetPosition(vec,"car");
    line->SetNextPoint(vec[0],vec[1],vec[2]);
    r=tx->GetEndPoint();
    if (!r) r=&dummy;
    r->GetPosition(vec,"car");
    line->SetNextPoint(vec[0],vec[1],vec[2]);
    line->SetLineWidth(1);
   }
  }
  else
  {
   r=tx->GetClosestPoint();
   if (!r) r=&dummy;
   r->GetPosition(vec,"car");
   line->SetNextPoint(vec[0],vec[1],vec[2]);
   p=tx->Get3Momentum();
   p=p+(*r);
   p.GetVector(vec,"car");
   line->SetNextPoint(vec[0],vec[1],vec[2]);
   line->SetLineWidth(3);
  }

  line->SetLineColor(kGreen);              // Neutral track
  if (charge>0) line->SetLineColor(kRed);  // Positive track
  if (charge<0) line->SetLineColor(kBlue); // Negative track
 
  // Mark tracks belonging to jets
  if (IsJetTrack(tx))
  {
   if (jets==1)
   {
    line->SetLineWidth(2);
    line->SetLineColor(kMagenta);
   }
  }
    
  line->Draw();
 }

 // Go for secondary vertices if selected
 if (secs==1)
 {
  Int_t nvtx=GetNvertices();
  for (Int_t jvtx=1; jvtx<=nvtx; jvtx++)
  {
   vx=GetVertex(jvtx);
   if (vx) vx->Draw(secs,cons,jets);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcVertex::SortJets(Int_t mode,TObjArray* jets)
{
// Order the references to an array of jets by looping over the input array "jets"
// and checking the value of a certain observable.
// The ordered array is returned as a TObjArray.
// In case jets=0 (default), the registered jets of the current vertex are used. 
// Note that the original jet array is not modified.
// Via the "mode" argument the user can specify the observable to be checked upon
// and specify whether sorting should be performed in decreasing order (mode<0)
// or in increasing order (mode>0).
//
// The convention for the observable selection is the following :
// mode : 1 ==> Number of tracks in the jet
//        2 ==> Jet energy
//        3 ==> Jet momentum
//        4 ==> Invariant mass of the jet
//        5 ==> Transverse momentum of the jet
//        6 ==> Longitudinal momentum of the jet
//        7 ==> Transverse energy of the jet
//        8 ==> Longitudinal energy of the jet
//        9 ==> Transverse mass of the jet
//       10 ==> Jet rapidity
//       11 ==> Pseudo-rapidity of the jet
//       12 ==> Number of associated signals
//       13 ==> Total charge of the jet
//
// The default is mode=-1.
//
// Note : This sorting routine uses a common area in memory, which is used
//        by various other sorting facilities as well.
//        This means that the resulting sorted TObjArray may be overwritten
//        when another sorting is invoked.
//        To retain the sorted list of pointers, the user is advised to copy
//        the pointers contained in the returned TObjArray into a private
//        TObjArray instance.

 if (fSelected)
 {
  delete fSelected;
  fSelected=0;
 }

 if (!jets) jets=fJets;
 
 if (!mode || abs(mode)>13 || !jets) return fSelected;

 Int_t njets=jets->GetEntries();
 if (!njets)
 {
  return fSelected;
 }
 else
 {
  fSelected=new TObjArray(njets);
 }

 Double_t val1,val2; // Values of the observable to be tested upon
 
 Int_t nord=0;
 for (Int_t i=0; i<njets; i++) // Loop over all jets of the array
 {
  NcJet* jx=(NcJet*)jets->At(i);

  if (!jx) continue;
 
  if (nord == 0) // store the first jet at the first ordered position
  {
   nord++;
   fSelected->AddAt(jx,nord-1);
   continue;
  }
 
  for (Int_t j=0; j<=nord; j++) // put jet in the right ordered position
  {
   if (j == nord) // jet has smallest (mode<0) or largest (mode>0) observable value seen so far
   {
    nord++;
    fSelected->AddAt(jx,j); // add jet at the end
    break; // go for next jet
   }

   val1=0;
   val2=0;
   
   switch (abs(mode))
   {
    case 1:
     val1=jx->GetNtracks();
     val2=((NcJet*)fSelected->At(j))->GetNtracks();
     break;
    case 2:
     val1=jx->GetEnergy(1);
     val2=((NcJet*)fSelected->At(j))->GetEnergy(1);
     break;
    case 3:
     val1=jx->GetMomentum(1);
     val2=((NcJet*)fSelected->At(j))->GetMomentum(1);
     break;
    case 4:
     val1=jx->GetInvmass(1);
     val2=((NcJet*)fSelected->At(j))->GetInvmass(1);
     break;
    case 5:
     val1=jx->GetPt(1);
     val2=((NcJet*)fSelected->At(j))->GetPt(1);
     break;
    case 6:
     val1=jx->GetPl(1);
     val2=((NcJet*)fSelected->At(j))->GetPl(1);
     break;
    case 7:
     val1=jx->GetEt(1);
     val2=((NcJet*)fSelected->At(j))->GetEt(1);
     break;
    case 8:
     val1=jx->GetEl(1);
     val2=((NcJet*)fSelected->At(j))->GetEl(1);
     break;
    case 9:
     val1=jx->GetMt(1);
     val2=((NcJet*)fSelected->At(j))->GetMt(1);
     break;
    case 10:
     val1=jx->GetRapidity();
     val2=((NcJet*)fSelected->At(j))->GetRapidity();
     break;
    case 11:
     val1=jx->GetPseudoRapidity();
     val2=((NcJet*)fSelected->At(j))->GetPseudoRapidity();
     break;
    case 12:
     val1=jx->GetNsignals();
     val2=((NcJet*)fSelected->At(j))->GetNsignals();
     break;
    case 13:
     val1=jx->GetCharge();
     val2=((NcJet*)fSelected->At(j))->GetCharge();
     break;
   }

   if (mode<0 && val1 <= val2) continue;
   if (mode>0 && val1 >= val2) continue;
 
   nord++;
   for (Int_t k=nord-1; k>j; k--) // create empty position
   {
    fSelected->AddAt(fSelected->At(k-1),k);
   }
   fSelected->AddAt(jx,j); // put jet at empty position
   break; // go for next jet
  }
 }
 return fSelected;
}
///////////////////////////////////////////////////////////////////////////
TObject* NcVertex::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. NcEvent to store either NcVertex objects or
// objects derived from NcVertex via the AddVertex memberfunction, provided
// these derived classes also have a proper Clone memberfunction. 

 NcVertex* vtx=new NcVertex(*this);
 if (name)
 {
  if (strlen(name)) vtx->SetName(name);
 }
 return vtx;
}
///////////////////////////////////////////////////////////////////////////
