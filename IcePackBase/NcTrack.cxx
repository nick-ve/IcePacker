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
// Class NcTrack
// Handling of the attributes of a reconstructed particle track.
//
// Coding example :
// ----------------
//
// Float_t a[4]={195.,1.2,-0.04,8.5};
// Nc4Vector pmu;
// pmu.SetVector(a,"car");
// NcTrack t1;
// t1.Set4Momentum(pmu);
//
// Float_t b[3]={1.2,-0.04,8.5};
// Nc3Vector p;
// p.SetVector(b,"car");
// NcTrack t2;
// t2.Set3Momentum(p);
// t2.SetCharge(0);
// t2.SetMass(1.115);
//
// t1.Data();
// t2.Data();
//
// Float_t pi=acos(-1.);
// Float_t thcms=0.2*pi; // decay theta angle in cms
// Float_t phicms=pi/4.; // decay theta angle in cms
// Float_t m1=0.938;
// Float_t m2=0.140;
// t2.Decay(m1,m2,thcms,phicms); // Track t2 decay : Lambda -> proton + pion
//
// t2.List();
//
// Int_t ndec=t2.GetNdecay();
// NcTrack* d1=t2.GetDecayTrack(1); // Access to decay track number 1
// NcTrack* d2=t2.GetDecayTrack(2); // Access to decay track number 2
//
// NcSignal s1,s2,s3,s4;
//
// .... // Code (e.g. detector readout) to fill NcSignal data
//
// NcTrack trec; // Track which will be reconstructed from signals
// trec.AddSignal(s1);
// trec.AddSignal(s3);
// trec.AddSignal(s4);
//
// Nc3Vector P;
// Float_t Q,M;
//
// ... // Code which accesses signals from trec and reconstructs
//        3-momentum P, charge Q, mass M etc...
//
// trec.Set3Momentum(P);
// trec.SetCharge(Q);
// trec.SetMass(M);
//
// Float_t r1[3]={1.6,-3.8,25.7};
// Float_t er1[3]={0.2,0.5,1.8};
// Float_t r2[3]={8.6,23.8,-6.7};
// Float_t er2[3]={0.93,1.78,0.8};
// NcPosition begin,end;
// begin.SetPosition(r1,"car");
// begin.SetPositionErrors(er1,"car");
// end.SetPosition(r2,"car");
// end.SetPositionErrors(er2,"car");
// trec.SetBeginPoint(begin);
// trec.SetEndPoint(end);
// 
// Note : By default all quantities are in GeV, GeV/c or GeV/c**2
//        but the user can indicate the usage of a different scale
//        for the energy-momentum units via the SetEscale() memberfunction.
//        The actual energy-momentum unit scale can be obtained via the
//        GetEscale() memberfunction.
//
//--- Author: Nick van Eijndhoven 10-jul-1997 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcTrack.h"
#include "Riostream.h"
 
ClassImp(NcTrack) // Class implementation to enable ROOT I/O
 
NcTrack::NcTrack() : TNamed(),Nc4Vector()
{
// Default constructor
// All variables initialised to 0
 Init();
 Reset();
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Init()
{
// Initialisation of pointers etc...
 fDecays=0;
 fSignals=0;
 fHypotheses=0;
 fBegin=0;
 fEnd=0;
 fRef=0;
 fImpactXY=0;
 fImpactXZ=0;
 fImpactYZ=0;
 fClosest=0;
 fParent=0;
 fFit=0;
 fTstamp=0;
 fEscale=1;
}
///////////////////////////////////////////////////////////////////////////
NcTrack::~NcTrack()
{
// Destructor to delete memory allocated for decay tracks array.
// This destructor automatically cleares all references to this NcTrack
// from all the related NcSignal objects.
 Int_t nsig=GetNsignals();
 for (Int_t i=1; i<=nsig; i++)
 {
  NcSignal* s=GetSignal(i);
  if (s) s->RemoveTrack(*this,0);
 }
 
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  delete fSignals;
  fSignals=0;
 }
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
 if (fBegin)
 {
  delete fBegin;
  fBegin=0;
 }
 if (fEnd)
 {
  delete fEnd;
  fEnd=0;
 }
 if (fRef)
 {
  delete fRef;
  fRef=0;
 }
 if (fImpactXY)
 {
  delete fImpactXY;
  fImpactXY=0;
 }
 if (fImpactXZ)
 {
  delete fImpactXZ;
  fImpactXZ=0;
 }
 if (fImpactYZ)
 {
  delete fImpactYZ;
  fImpactYZ=0;
 }
 if (fClosest)
 {
  delete fClosest;
  fClosest=0;
 }
 if (fFit)
 {
  delete fFit;
  fFit=0;
 }
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
NcTrack::NcTrack(const NcTrack& t) : TNamed(t),Nc4Vector(t)
{
// Copy constructor
 Init();

 fQ=t.fQ;
 fProb=t.fProb;
 if (t.fBegin) fBegin=new NcPositionObj(*(t.fBegin));
 if (t.fEnd) fEnd=new NcPositionObj(*(t.fEnd));
 if (t.fRef) fRef=new NcPositionObj(*(t.fRef));
 if (t.fImpactXY) fImpactXY=new NcPositionObj(*(t.fImpactXY));
 if (t.fImpactXZ) fImpactXZ=new NcPositionObj(*(t.fImpactXZ));
 if (t.fImpactYZ) fImpactYZ=new NcPositionObj(*(t.fImpactYZ));
 if (t.fClosest) fClosest=new NcPositionObj(*(t.fClosest));
 if (t.fFit) fFit=t.fFit->Clone();
 if (t.fTstamp) fTstamp=new NcTimestamp(*(t.fTstamp));
 fUserId=t.fUserId;
 fEscale=t.fEscale;
 fCode=t.fCode;
 fParent=t.fParent;

 Int_t ndec=t.GetNdecay();
 if (ndec)
 {
  fDecays=new TObjArray(ndec);
  fDecays->SetOwner();
  for (Int_t it=1; it<=ndec; it++)
  {
   NcTrack* tx=t.GetDecayTrack(it);
   fDecays->Add(new NcTrack(*tx));
  }
 }

 Int_t nsig=t.GetNsignals();
 if (nsig)
 {
  fSignals=new TObjArray(nsig);
  for (Int_t is=1; is<=nsig; is++)
  {
   NcSignal* sx=t.GetSignal(is);
   fSignals->Add(sx);
  }
 }

 Int_t nhyp=t.GetNhypotheses();
 if (nhyp)
 {
  fHypotheses=new TObjArray(nhyp);
  fHypotheses->SetOwner();
  for (Int_t ih=1; ih<=nhyp; ih++)
  {
   NcTrack* tx=t.GetTrackHypothesis(ih);
   fHypotheses->Add(new NcTrack(*tx));
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Reset()
{
// Reset all variables to 0 and delete all auto-generated decay tracks.
// Note : The scale for the energy/momentum units will not be changed.
 fQ=0;
 fUserId=0;
 fCode=0;
 fProb=0;
 SetZero();
 fParent=0;
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 if (fSignals)
 {
  delete fSignals;
  fSignals=0;
 }
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
 if (fBegin)
 {
  delete fBegin;
  fBegin=0;
 }
 if (fEnd)
 {
  delete fEnd;
  fEnd=0;
 }
 if (fRef)
 {
  delete fRef;
  fRef=0;
 }
 if (fImpactXY)
 {
  delete fImpactXY;
  fImpactXY=0;
 }
 if (fImpactXZ)
 {
  delete fImpactXZ;
  fImpactXZ=0;
 }
 if (fImpactYZ)
 {
  delete fImpactYZ;
  fImpactYZ=0;
 }
 if (fClosest)
 {
  delete fClosest;
  fClosest=0;
 }
 if (fFit)
 {
  delete fFit;
  fFit=0;
 }
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Set3Momentum(Nc3Vector& p)
{
// Set the track parameters according to the 3-momentum p.
// In case the mass was not yet set, the energy is set to correspond to m=0. 
 Set3Vector(p);
 Double_t inv=GetInvariant();
 if (inv<0) SetMass(0.);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Set4Momentum(Nc4Vector& p)
{
// Set the track parameters according to the 4-momentum p
 Double_t E=p.GetScalar();
 Double_t dE=p.GetResultError();
 Nc3Vector pv=p.Get3Vector();
 SetVector(E,pv);
 SetScalarError(dE);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetMass(Double_t m,Double_t dm)
{
// Set the particle mass
// The default value for the error dm is 0.
 Double_t inv=pow(m,2);
 Double_t dinv=fabs(2.*m*dm);
 SetInvariant(inv,dinv);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetCharge(Float_t q)
{
// Set the particle charge
 fQ=q;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Data(TString f,TString u)
{
// Provide track information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Double_t m=GetMass();
 Double_t dm=GetResultError();
 const char* name=GetName();
 const char* title=GetTitle();

 cout << " *" << ClassName() << "::Data*";
 if (strlen(name))  cout << " Name : " << name;
 if (strlen(title)) cout << " Title : " << title;
 cout << endl;
 if (fTstamp) fTstamp->Date(1);
 cout << " Id : " << fUserId << " Code : " << fCode
      << " m : " << m << " dm : " << dm << " Charge : " << fQ
      << " p : " << GetMomentum() << endl;
 cout << " Nhypotheses : " << GetNhypotheses() << " Ndecay-tracks : " << GetNdecay()
      << " Nsignals : " << GetNsignals() << " Energy scale : " << fEscale << " GeV" << endl;
 if (fParent)
 {
  cout << " Parent track Id : " << fParent->GetId() << " Code : " << fParent->GetParticleCode()
       << " m : " << fParent->GetMass() << " Q : " << fParent->GetCharge()
       << " p : " << fParent->GetMomentum();
  const char* pname=fParent->GetName();
  const char* ptitle=fParent->GetTitle();
  if (strlen(pname))  cout << " Name : " << pname;
  if (strlen(ptitle)) cout << " Title : " << ptitle;
  cout << endl;
 }
 if (fFit)
 {
  cout << " Fit details present in object of class " << fFit->ClassName() << endl; 
  if (fFit->InheritsFrom("NcSignal")) ((NcSignal*)fFit)->List(-1);
 }
 Nc4Vector::Data(f,u); 
} 
///////////////////////////////////////////////////////////////////////////
void NcTrack::List(TString f,TString u)
{
// Provide current track and decay level 1 information within coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current track
 if (fBegin) { cout << " Begin-point :"; fBegin->Data(f,u); }
 if (fEnd)   { cout << " End-point   :"; fEnd->Data(f,u); }
 if (fRef)   { cout << " Ref-point   :"; fRef->Data(f,u); }

 // Decay products of this track
 NcTrack* td; 
 for (Int_t id=1; id<=GetNdecay(); id++)
 {
  td=GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level 1 sec. track no. " << id << endl;
   td->Data(f,u); 
  }
  else
  {
   cout << " *NcTrack::List* Error : Empty decay track slot." << endl; 
  }
 }
} 
///////////////////////////////////////////////////////////////////////////
void NcTrack::ListAll(TString f,TString u)
{
// Provide complete track and decay information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Data(f,u); // Information of the current track
 if (fBegin) { cout << " Begin-point :"; fBegin->Data(f,u); }
 if (fEnd)   { cout << " End-point   :"; fEnd->Data(f,u); }
 if (fRef)   { cout << " Ref-point   :"; fRef->Data(f,u); }

 Int_t nhyp=GetNhypotheses();
 if (nhyp)
 {
  cout << " List of the " << nhyp << " track hypotheses : " << endl;
  for (Int_t ih=1; ih<=nhyp; ih++)
  {
   NcTrack* tx=GetTrackHypothesis(ih);
   if (tx) tx->Data(f,u);
  }
 }

 Int_t nsig=GetNsignals();
 if (nsig)
 {
  cout << " List of the corresponding slots for the " << nsig
       << " related signals : " << endl;
  NcPosition r;
  Int_t nrefs,jslot;
  TArrayI slotarr;
  for (Int_t is=1; is<=nsig; is++)
  {
   NcSignal* sx=GetSignal(is);
   if (sx)
   {
    nrefs=sx->GetIndices(this,slotarr,0);
    for (Int_t jref=0; jref<nrefs; jref++)
    {
     jslot=slotarr.At(jref);
     sx->List(jslot);
    }
    r=sx->GetPosition();
    cout << "   Position";
    r.Data(f,u);
   }
  }
 }

 NcTrack* t=this;
 Dumps(t,1,f,u); // Information of all decay products
}
//////////////////////////////////////////////////////////////////////////
void NcTrack::Dumps(NcTrack* t,Int_t n,TString f,TString u)
{
// Recursively provide the info of all decay levels of this track
 NcTrack* td; 
 for (Int_t id=1; id<=t->GetNdecay(); id++)
 {
  td=t->GetDecayTrack(id);
  if (td)
  {
   cout << "  ---Level " << n << " sec. track no. " << id << endl;
   td->Data(f,u); 

   Int_t nhyp=td->GetNhypotheses();
   if (nhyp)
   {
    cout << " List of the " << nhyp << " track hypotheses : " << endl;
    for (Int_t ih=1; ih<=nhyp; ih++)
    {
     NcTrack* tx=td->GetTrackHypothesis(ih);
     if (tx) tx->Data(f,u);
    }
   }

   Int_t nsig=td->GetNsignals();
   if (nsig)
   {
    cout << " List of the " << nsig << " related signals : " << endl;
    for (Int_t is=1; is<=nsig; is++)
    {
     NcSignal* sx=td->GetSignal(is);
     if (sx) sx->Data(f,u);
    }
   }

   // Go for next decay level of this decay track recursively
   Dumps(td,n+1,f,u);
  }
  else
  {
   cout << " *NcTrack::Dumps* Error : Empty decay track slot." << endl; 
  }
 }
} 
//////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetMomentum(Float_t scale)
{
// Provide the value of the track 3-momentum.
// By default the momentum is returned in the units as it was stored in the track
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the momentum in MeV/c.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMomentum().

 Double_t norm=fV.GetNorm();
 fDresult=fV.GetResultError();
 if (scale>0)
 {
  norm*=fEscale/scale;
  fDresult*=fEscale/scale;
 }
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Nc3Vector NcTrack::Get3Momentum(Float_t scale) const
{
// Provide the track 3-momentum.
// By default the components of the 3-momentum are returned in the units
// as they were stored in the track structure.
// However, the user can select a different momentum unit scale for the
// components by specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the 3-momentum in MeV/c.

 Nc3Vector p=Get3Vector();
 if (scale>0) p*=fEscale/scale;
 return p;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetMass(Float_t scale)
{
// Provide the particle mass.
// By default the mass is returned in the units as it was stored in the track
// structure. However, the user can select a different mass unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c**2, so specification
// of scale=0.001 will provide the mass in MeV/c**2.
// The error can be obtained by invoking GetResultError() after
// invokation of GetMass().

 Double_t inv=GetInvariant();
 Double_t dinv=GetResultError();
 Double_t dm=0;
 if (inv >= 0)
 {
 Double_t m=sqrt(inv);
 if (m) dm=dinv/(2.*m);
 if (scale>0)
 {
  m*=fEscale/scale;
  dm*=fEscale/scale;
 }
 fDresult=dm;
 return m;
 }
 else
 {
  cout << "*NcTrack::GetMass* Unphysical situation m**2 = " << inv << endl;
  cout << " Value 0 will be returned." << endl;
  fDresult=dm;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t NcTrack::GetCharge() const
{
// Provide the particle charge
 return fQ;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetEnergy(Float_t scale)
{
// Provide the particle's energy.
// By default the energy is returned in the units as it was stored in the track
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the energy in MeV.
// The error can be obtained by invoking GetResultError() after
// invokation of GetEnergy().
 Double_t E=GetScalar();
 if (E>0)
 {
  if (scale>0)
  {
   E*=fEscale/scale;
   fDresult*=fEscale/scale;
  }
  return E;
 }
 else
 {
  cout << "*NcTrack::GetEnergy* Unphysical situation E = " << E << endl;
  cout << " Value 0 will be returned." << endl;
  return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::Decay(Double_t m1,Double_t m2,Double_t thcms,Double_t phicms)
{
// Perform 2-body decay of current track
// m1     : mass of decay product 1
// m2     : mass of decay product 2
// thcms  : cms theta decay angle (in rad.) of m1
// phicms : cms phi decay angle (in rad.) of m1
 
 Double_t M=GetMass();
 
// Compute the 4-momenta of the decay products in the cms
// Note : p2=p1=pnorm for a 2-body decay
 Double_t e1=0;
 if (M) e1=((M*M)+(m1*m1)-(m2*m2))/(2.*M);
 Double_t e2=0;
 if (M) e2=((M*M)+(m2*m2)-(m1*m1))/(2.*M);
 Double_t pnorm=(e1*e1)-(m1*m1);
 if (pnorm>0.)
 {
  pnorm=sqrt(pnorm);
 }
 else
 {
  pnorm=0;
 }
 
 Double_t a[3];
 a[0]=pnorm;
 a[1]=thcms;
 a[2]=phicms;
 Nc3Vector p;
 p.SetVector(a,"sph");

 Nc4Vector pprim1;
 pprim1.SetVector(e1,p);
 pprim1.SetInvariant(m1*m1);

 Nc4Vector pprim2;
 p*=-1;
 pprim2.SetVector(e2,p);
 pprim2.SetInvariant(m2*m2);

 // Determine boost parameters from the parent particle
 Double_t E=GetEnergy();
 p=Get3Vector();
 Nc4Vector pmu;
 pmu.SetVector(E,p);

 NcBoost q;
 q.Set4Momentum(pmu);
 
 Nc4Vector p1=q.Inverse(pprim1); // Boost decay product 1
 Nc4Vector p2=q.Inverse(pprim2); // Boost decay product 2
 
 // Enter the boosted data into the decay tracks array
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
 fDecays=new TObjArray(2);
 fDecays->SetOwner();

 fDecays->Add(new NcTrack);
 ((NcTrack*)fDecays->At(0))->Set4Momentum(p1);
 ((NcTrack*)fDecays->At(0))->SetMass(m1);
 fDecays->Add(new NcTrack);
 ((NcTrack*)fDecays->At(1))->Set4Momentum(p2);
 ((NcTrack*)fDecays->At(1))->SetMass(m2);
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetNdecay() const
{
// Provide the number of decay produced tracks
 Int_t ndec=0;
 if (fDecays) ndec=fDecays->GetEntries();
 return ndec;
}
///////////////////////////////////////////////////////////////////////////
NcTrack* NcTrack::GetDecayTrack(Int_t j) const
{
// Provide decay produced track number j
// Note : j=1 denotes the first decay track
 if (!fDecays)
 {
  cout << " *NcTrack::GetDecayTrack* No tracks present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= GetNdecay()))
  {
   return (NcTrack*)fDecays->At(j-1);
  }
  else
  {
   cout << " *NcTrack* decay track number : " << j << " out of range."
        << " Ndec = " << GetNdecay() << endl;
   return 0;  
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveDecays()
{
// Remove all decay tracks from this track.
 if (fDecays)
 {
  delete fDecays;
  fDecays=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::AddSignal(NcSignal& s,Int_t mode)
{
// Relate an NcSignal object to this track.
//
// mode = 0 : Only the reference to the specified signal is stored in
//            the current track, without storing the (backward) reference
//            to this track into the NcSignal structure. 
//        1 : The (backward) reference to the current track is also automatically
//            stored into the NcSignal (or derived) object specified in the
//            input argument.
//
// The default is mode=0.

 if (!fSignals) fSignals=new TObjArray(1);

 // Check if this signal is already stored for this track
 Int_t nsig=GetNsignals();
 for (Int_t i=0; i<nsig; i++)
 {
  if (&s==fSignals->At(i)) return; 
 }

 fSignals->Add(&s);
 if (mode==1) s.AddTrack(*this,0);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveSignal(NcSignal& s,Int_t mode)
{
// Remove related NcSignal object from this track.
//
// mode = 0 : Only the reference to the specified signal is removed from
//            the current track, without removing the (backward) reference(s)
//            to this track from the NcSignal structure. 
//        1 : The (backward) reference(s) to the current track are also automatically
//            removed from the NcSignal (or derived) object specified in the
//            input argument.
//
// The default is mode=1.

 if (fSignals)
 {
  NcSignal* test=(NcSignal*)fSignals->Remove(&s);
  if (test) fSignals->Compress();
 }
 if (mode==1) s.RemoveTrack(*this,0);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveSignals(Int_t mode)
{
// Remove all related NcSignal objects from this track.
//
// mode = 0 : All signal references are removed from the current track,
//            without removing the (backward) references to this track from
//            the corresponding NcSignal objects.
//        1 : The (backward) references to the current track are also automatically
//            removed from the corresponding NcSignal (or derived) objects.
//
// The default is mode=1.

 if (!fSignals) return;

 Int_t ns=GetNsignals();
 for (Int_t i=0; i<ns; i++)
 {
  NcSignal* sx=(NcSignal*)fSignals->At(i);
  if (sx && mode==1) sx->RemoveTrack(*this,0);
 }

 delete fSignals;
 fSignals=0;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetNsignals() const
{
// Provide the number of related NcSignals.
 Int_t nsig=0;
 if (fSignals) nsig=fSignals->GetEntries();
 return nsig;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetNsignals(const char* classname,Int_t par) const
{
// Provide the number of stored signals (derived) of the specified class.
//
// par = 0 ==> The signal itself has to be (derived) of the specified class
//       1 ==> The parent device of the signal has to be (derived) of the specified class
//       2 ==> The signal or the parent device has to be (derived) of the specified class
//
// The default is par=0.

 Int_t nsigs=0;
 for (Int_t isig=1; isig<=GetNsignals(); isig++)
 {
  NcSignal* sx=GetSignal(isig);
  if (!sx) continue;

  if ((par==0 || par==2) && sx->InheritsFrom(classname))
  {
   nsigs++;
   continue;
  }
  
  NcDevice* parent=sx->GetDevice();
  if (!parent) continue;

  if ((par==1 || par==2) && parent->InheritsFrom(classname)) nsigs++;
 }
 return nsigs;
} 
///////////////////////////////////////////////////////////////////////////
NcSignal* NcTrack::GetSignal(Int_t j) const
{
// Provide the related NcSignal number j.
// Note : j=1 denotes the first signal.
 if (!fSignals)
 {
  cout << " *NcTrack::GetSignal* No signals present." << endl;
  return 0;
 }
 else
 {
  if ((j >= 1) && (j <= GetNsignals()))
  {
   return (NcSignal*)fSignals->At(j-1);
  }
  else
  {
   cout << " *NcTrack* signal number : " << j << " out of range."
        << " Nsig = " << GetNsignals() << endl;
   return 0;
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::AddTrackHypothesis(NcTrack& t)
{
// Relate a track hypothesis to this track.
// Note : a private copy of the input track will be made via the Clone()
//        facility.
 if (!fHypotheses)
 {
  fHypotheses=new TObjArray(1);
  fHypotheses->SetOwner();
 }
 fHypotheses->Add(t.Clone());
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::AddTrackHypothesis(Double_t prob,Double_t m,Double_t dm)
{
// Add a track hypothesis by explicitly setting the mass and probability.
// This will affect e.g. the hypothesis track's energy, since the momentum
// and all other attributes will be copied from the current track.
//
// Input arguments :
// ----------------- 
// prob=probalility  m=mass value  dm=error on the mass value.
// The default value for the mass error dm is 0.

 NcTrack t(*this);
 t.RemoveDecays();
 t.RemoveTrackHypotheses();
 t.RemoveSignals();
 t.SetTitle("Mass hypothesis");
 t.SetMass(m,dm);
 t.SetProb(prob);
 AddTrackHypothesis(t);
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveTrackHypothesis(NcTrack& t)
{
// Remove the specified track hypothesis from this track.
 if (fHypotheses)
 {
  NcTrack* test=(NcTrack*)fHypotheses->Remove(&t);
  if (test) fHypotheses->Compress();
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveTrackHypotheses()
{
// Remove all track hypotheses from this track.
 if (fHypotheses)
 {
  delete fHypotheses;
  fHypotheses=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetNhypotheses() const
{
// Provide the number of track hypotheses.
 Int_t nhyp=0;
 if (fHypotheses) nhyp=fHypotheses->GetEntries();
 return nhyp;
}
///////////////////////////////////////////////////////////////////////////
NcTrack* NcTrack::GetTrackHypothesis(Int_t j) const
{
// Provide the j-th track hypothesis.
// Note : j=1 denotes the first hypothesis.
// Default : j=0 ==> Hypothesis with highest probability.

 if (!fHypotheses) return 0;

 Int_t nhyp=GetNhypotheses();

 // Check validity of index j
 if (j<0 || j>nhyp)
 {
   cout << " *NcTrack* hypothesis number : " << j << " out of range."
        << " Nhyp = " << nhyp << endl;
   return 0;
 } 

 NcTrack* t=0;

 if (j==0) // Provide track hypothesis with highest probability
 {
  Float_t prob=0;   
  t=(NcTrack*)fHypotheses->At(0);
  if (t) prob=t->GetProb();
  Float_t probx=0;
  for (Int_t ih=1; ih<nhyp; ih++)
  {
   NcTrack* tx=(NcTrack*)fHypotheses->At(ih);
   if (tx)
   {
    probx=tx->GetProb();
    if (probx > prob) t=tx; 
   }
  }
  return t;
 }
 else // Provide requested j-th track hypothesis
 {
  return (NcTrack*)fHypotheses->At(j-1);
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetBeginPoint(NcPosition& p)
{
// Store the position of the track begin-point.
 if (fBegin) delete fBegin;
 fBegin=new NcPositionObj(p);
}
///////////////////////////////////////////////////////////////////////////
NcPosition* NcTrack::GetBeginPoint()
{
// Provide the position of the track begin-point.
 return fBegin;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetEndPoint(NcPosition& p)
{
// Store the position of the track end-point.
 if (fEnd) delete fEnd;
 fEnd=new NcPositionObj(p);
}
///////////////////////////////////////////////////////////////////////////
NcPosition* NcTrack::GetEndPoint()
{
// Provide the position of the track end-point.
 return fEnd;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetReferencePoint(NcPosition& p)
{
// Store the position of the track reference-point.
// The reference-point is the point on the track in which the 
// 3-momentum vector components have been defined.
// This reference point is the preferable point to start track extrapolations
// etc... which are sensitive to the components of the 3-momentum vector.
 if (fRef) delete fRef;
 fRef=new NcPositionObj(p);
}
///////////////////////////////////////////////////////////////////////////
NcPosition* NcTrack::GetReferencePoint()
{
// Provide the position of the track reference-point.
// The reference-point is the point on the track in which the 
// 3-momentum vector components have been defined.
// This reference point is the preferable point to start track extrapolations
// etc... which are sensitive to the components of the 3-momentum vector.
 return fRef;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetMass()
{
// Set the mass and error to the value of the hypothesis with highest prob.

 Double_t m=0,dm=0;

 // Select mass hypothesis with highest probability
 NcTrack* t=GetTrackHypothesis(0);
 if (t) 
 {
  m=t->GetMass();
  dm=t->GetResultError();
  SetMass(m,dm);
 }
 else
 {
  cout << " *NcTrack::SetMass()* No hypothesis present => No action." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetPt(Float_t scale)
{
// Provide the transverse momentum value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the track
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the transverse momentum in MeV/c.
// The error on the value can be obtained by GetResultError()
// after invokation of GetPt().
 Nc3Vector v;
 v=GetVecTrans();
 Double_t norm=v.GetNorm();
 fDresult=v.GetResultError();
 if (scale>0)
 {
  norm*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return norm;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetPl(Float_t scale)
{
// Provide the longitudinal momentum value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the track
// structure. However, the user can select a different momentum unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV/c, so specification
// of scale=0.001 will provide the longitudinal momentum in MeV/c.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetPl().

 Nc3Vector v;
 v=GetVecLong();

 Double_t pl=v.GetNorm();
 fDresult=v.GetResultError();

 Double_t a[3];
 v.GetVector(a,"sph");
 if (cos(a[1])<0) pl=-pl;
 if (scale>0)
 {
  pl*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return pl;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetEt(Float_t scale)
{
// Provide transverse energy value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the track
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the transverse energy in MeV.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEt().

 Double_t et=GetScaTrans();
 if (scale>0)
 {
  et*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return et;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetEl(Float_t scale)
{
// Provide longitudinal energy value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the track
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the longitudinal energy in MeV.
// Note : the returned value can also be negative.
// The error on the value can be obtained by GetResultError()
// after invokation of GetEl().

 Double_t el=GetScaLong();
 if (scale>0)
 {
  el*=fEscale/scale;
  fDresult*=fEscale/scale;
 }

 return el;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetMt(Float_t scale)
{
// Provide transverse mass value w.r.t. z-axis.
// By default the value is returned in the units as it was stored in the track
// structure. However, the user can select a different energy unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to GeV, so specification
// of scale=0.001 will provide the transverse mass in MeV.
// The error on the value can be obtained by GetResultError()
// after invokation of GetMt().
 Double_t pt=GetPt();
 Double_t dpt=GetResultError();
 Double_t m=GetMass();
 Double_t dm=GetResultError();

 Double_t mt=sqrt(pt*pt+m*m);
 Double_t dmt2=0;
 if (mt) dmt2=(pow((pt*dpt),2)+pow((m*dm),2))/(mt*mt);

 fDresult=sqrt(dmt2);
 if (scale>0)
 {
  mt*=fEscale/scale;
  fDresult*=fEscale/scale;
 }
 return mt;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetRapidity()
{
// Provide rapidity value w.r.t. z-axis.
// The error on the value can be obtained by GetResultError()
// after invokation of GetRapidity().
// Note : Also GetPseudoRapidity() is available since this class is
//        derived from Nc4Vector.
 Double_t e=GetEnergy();
 Double_t de=GetResultError();
 Double_t pl=GetPl();
 Double_t dpl=GetResultError();
 Double_t sum=e+pl;
 Double_t dif=e-pl;

 Double_t y=9999,dy2=0;
 if (sum && dif) y=0.5*log(sum/dif);

 if (sum*dif) dy2=(1./(sum*dif))*(pow((pl*de),2)+pow((e*dpl),2));

 fDresult=sqrt(dy2);
 return y;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetImpactPoint(NcPosition& p,TString q)
{
// Store the position of the impact-point in the plane "q=0".
// Here q denotes one of the axes X, Y or Z.
// Note : The character to denote the axis may be entered in lower or
//        in uppercase.
 Int_t axis=0;
 if (q=="x" || q=="X") axis=1;
 if (q=="y" || q=="Y") axis=2;
 if (q=="z" || q=="Z") axis=3;

 switch (axis)
 {
  case 1: // Impact-point in the plane X=0
   if (fImpactYZ) delete fImpactYZ;
   fImpactYZ=new NcPositionObj(p);
   break;

  case 2: // Impact-point in the plane Y=0
   if (fImpactXZ) delete fImpactXZ;
   fImpactXZ=new NcPositionObj(p);
   break;

  case 3: // Impact-point in the plane Z=0
   if (fImpactXY) delete fImpactXY;
   fImpactXY=new NcPositionObj(p);
   break;

  default: // Unsupported axis
   cout << "*NcTrack::SetImpactPoint* Unsupported axis : " << q << endl
        << " Possible axes are 'X', 'Y' and 'Z'." << endl; 
   break;
 }
}
///////////////////////////////////////////////////////////////////////////
NcPosition* NcTrack::GetImpactPoint(TString q)
{
// Provide the position of the impact-point in the plane "q=0".
// Here q denotes one of the axes X, Y or Z.
// Note : The character to denote the axis may be entered in lower or
//        in uppercase.
 Int_t axis=0;
 if (q=="x" || q=="X") axis=1;
 if (q=="y" || q=="Y") axis=2;
 if (q=="z" || q=="Z") axis=3;

 switch (axis)
 {
  case 1: // Impact-point in the plane X=0
   return fImpactYZ;

  case 2: // Impact-point in the plane Y=0
   return fImpactXZ;

  case 3: // Impact-point in the plane Z=0
   return fImpactXY;

  default: // Unsupported axis
   cout << "*NcTrack::GetImpactPoint* Unsupported axis : " << q << endl
        << " Possible axes are 'X', 'Y' and 'Z'." << endl; 
   return 0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetId(Int_t id)
{
// Set a user defined unique identifier for this track.
 fUserId=id;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetId() const
{
// Provide the user defined unique identifier of this track.
 return fUserId;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetClosestPoint(NcPosition& p)
{
// Set position p as the point of closest approach w.r.t. some reference
 if (fClosest) delete fClosest;
 fClosest=new NcPositionObj(p);
}
///////////////////////////////////////////////////////////////////////////
NcPosition* NcTrack::GetClosestPoint()
{
// Provide the point of closest approach w.r.t. some reference
 return fClosest;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetEscale(Float_t scale)
{
// Indicate the energy/momentum scale as used by the user.
// The convention is that scale=1 indicates values in units
// of GeV, GeV/c or GeV/c**2.
// So, in case one decides to store values in units of MeV, MeV/c or MeV/c**2
// the scale indicator should be set to scale=0.001.
//
// By default scale=1 is set in the constructor.

 if (scale>0)
 {
  fEscale=scale;
 }
 else
 {
  cout << " *NcTrack::SetEscale* Invalid scale value : " << scale << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t NcTrack::GetEscale() const
{
// Provide the energy/momentum scale as used by the user.
// The convention is that scale=1 indicates values in units
// of GeV, GeV/c or GeV/c**2.
// So, a value of scale=0.001 indicates that energy/momentum values are
// stored in units of MeV, MeV/c or MeV/c**2.
 return fEscale;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetParticleCode(Int_t code)
{
// Set the user defined particle id code (e.g. the PDF convention).
 fCode=code;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcTrack::GetParticleCode() const
{
// Provide the user defined particle id code.
 return fCode;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetParentTrack(NcTrack* t)
{
// Set pointer to the parent track.
 fParent=t;
}
///////////////////////////////////////////////////////////////////////////
NcTrack* NcTrack::GetParentTrack()
{
// Provide pointer to the parent track.
 return fParent;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetProb(Double_t prob)
{
// Set hypothesis probability for this track.
 fProb=prob;
}
///////////////////////////////////////////////////////////////////////////
Float_t NcTrack::GetProb() const
{
// Provide the hypothesis probability for this track.
 return fProb;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetFitDetails(TObject* obj)
{
// Enter the object containing the fit details.
// In case an object to hold fit details was already present, this
// will be deleted first before the new one is stored.
// This means that SetFitDetails(0) can be used to just remove the
// existing object with the fit details.
// All objects derived from TObject can be entered in this way.
// Obvious candidates for objects containing detailed fit information
// are functions (e.g. TF1) and histograms (e.g. TH1F).
// However, using an NcDevice object provides a very versatile facility
// to store the parameters of various fit procedures.
// In such a case the NcDevice can be used to provide the various fit
// definitions and the corresponding fit parameters can be entered as
// separate NcSignal objects which are stored as hits to the NcDevice.
// In addition various functions and histograms can be linked to the
// various NcSignal instances
// The latter procedure is based on the original idea of Adam Bouchta.
//
// Note : The entered object is owned by this NcTrack instance.
//        As such, a private copy of obj will be stored using the Clone()
//        memberfunction.
//        In case the entered object contains pointers to other objects,
//        the user has to provide the appropriate Clone() memberfunction
//        for the class to which the entered object belongs.
//        An example can be seen from NcTrack::Clone().   
//
 if (fFit)
 {
  delete fFit;
  fFit=0;
 }

 if (obj) fFit=obj->Clone();
}
///////////////////////////////////////////////////////////////////////////
TObject* NcTrack::GetFitDetails()
{
// Provide the pointer to the object containing the fit details.
 return fFit;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::SetTimestamp(NcTimestamp& t)
{
// Store the timestamp for this track.
 if (fTstamp) delete fTstamp;
 fTstamp=new NcTimestamp(t);
}
///////////////////////////////////////////////////////////////////////////
NcTimestamp* NcTrack::GetTimestamp()
{
// Provide the timestamp of this track.
 return fTstamp;
}
///////////////////////////////////////////////////////////////////////////
void NcTrack::RemoveTimestamp()
{
// Remove the timestamp from this track.
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetDistance(NcPosition* p,Float_t scale)
{
// Provide distance of the current track to the position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the NcPosition p.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the track parameters have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
// Note : In case of incomplete information, a distance value of -1 is
//        returned.
 
 Double_t dist=-1.;
 fDresult=0.;

 if (!p) return dist;

 // Obtain a defined position on this track
 NcPosition* rx=fRef;
 if (!rx) rx=fBegin;
 if (!rx) rx=fEnd;

 if (!rx) return dist;

 Nc3Vector p1=Get3Momentum();

 if (p1.GetNorm() <= 0.) return dist;

 Nc3Vector r0=(Nc3Vector)(*rx);

 Float_t tscale=rx->GetUnitScale();
 Float_t pscale=p->GetUnitScale();
 if ((tscale/pscale > 1.1) || (pscale/tscale > 1.1)) r0=r0*(tscale/pscale);
 
 // Obtain the direction unit vector of this track
 Double_t vec[3];
 Double_t err[3];
 p1.GetVector(vec,"sph");
 p1.GetErrors(err,"sph");
 vec[0]=1.;
 err[0]=0.;
 p1.SetVector(vec,"sph");
 p1.SetErrors(err,"sph");

 Nc3Vector q=(Nc3Vector)(*p);
 Nc3Vector r=q-r0;
 Nc3Vector d=r.Cross(p1);
 dist=d.GetNorm();
 fDresult=d.GetResultError();
 if (scale>0)
 {
  dist*=pscale/scale;
  fDresult*=pscale/scale;
 }
 return dist;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcTrack::GetDistance(NcTrack* t,Float_t scale)
{
// Provide distance of the current track to the track t.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the current track.
// This implies that the results of t1.GetDistance(t2) and t2.GetDistance(t1)
// may be numerically different in case t1 and t2 have different metric units.
// However, the user can specify a required metric unit scale by specification
// of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the track parameters have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
// Note : In case of incomplete information, a distance value of -1 is
//        returned.
 
 Double_t dist=-1.;
 fDresult=0.;

 if (!t) return dist;

 // Obtain a defined position on this track
 NcPosition* rx=fRef;
 if (!rx) rx=fBegin;
 if (!rx) rx=fEnd;

 if (!rx) return dist;

 // Obtain a defined position on track t
 NcPosition* ry=t->GetReferencePoint();
 if (!ry) ry=t->GetBeginPoint();
 if (!ry) ry=t->GetEndPoint();

 if (!ry) return dist;
 
 Nc3Vector p1=Get3Momentum();
 Nc3Vector p2=t->Get3Momentum();

 if (p1.GetNorm() <= 0. || p2.GetNorm() <= 0.) return dist;

 // The vector normal to both track directions
 Nc3Vector n=p1.Cross(p2);

 Float_t scalex=rx->GetUnitScale();
 Float_t scaley=ry->GetUnitScale();

 if (n.GetNorm() > 1.e-10)
 {
  // Normalise n to a unit vector
  Double_t vec[3];
  Double_t err[3];
  n.GetVector(vec,"sph");
  n.GetErrors(err,"sph");
  vec[0]=1.;
  err[0]=0.;
  n.SetVector(vec,"sph");
  n.SetErrors(err,"sph");
  Nc3Vector r1=(Nc3Vector)(*rx);
  Nc3Vector r2=(Nc3Vector)(*ry);
  // Correct components of r2 in case of different unit scales
  if ((scaley/scalex > 1.1) || (scalex/scaley > 1.1)) r2=r2*(scaley/scalex);
  Nc3Vector r=r1-r2;
  dist=fabs(r.Dot(n));
  fDresult=r.GetResultError();
 }
 else // Parallel tracks
 {
  dist=t->GetDistance(rx);
  fDresult=t->GetResultError();
 }

 if (scale>0)
 {
  dist*=scalex/scale;
  fDresult*=scalex/scale;
 }
 return dist;
}
///////////////////////////////////////////////////////////////////////////
TObject* NcTrack::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows e.g. NcJet to store either NcTrack objects or
// objects derived from NcTrack via the AddTrack memberfunction, provided
// these derived classes also have a proper Clone memberfunction. 

 NcTrack* trk=new NcTrack(*this);
 if (name)
 {
  if (strlen(name)) trk->SetName(name);
 }
 return trk;
}
///////////////////////////////////////////////////////////////////////////
