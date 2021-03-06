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
// Class NcEvent
// Creation and investigation of an NCFS generic event structure.
// An NcEvent can be constructed by adding NcTracks, NcVertices, NcJets
// and/or devices like NcCalorimeters or NcDevice (derived) objects.
//
// All objects which are derived from TObject can be regarded as a device.
// However, NcDevice (or derived) objects profit from additional hit
// handling facilities.
// A "hit" is a generic name indicating an NcSignal (or derived) object.
// Note that NcEvent does NOT own hits; it only provides references to hits
// obtained from the various devices.
// This implies that hits should be owned by the devices themselves.
//
// The basic functionality of NcEvent is identical to the one of NcVertex.
// So, an NcEvent may be used as the primary vertex with some additional
// functionality compared to NcVertex.
//
// To provide maximal flexibility to the user, the two modes of track/jet/vertex
// storage as described in NcJet and NcVertex can be used.
// In addition an identical structure is provided for the storage of devices like
// NcCalorimeter objects, which can be selected by means of the memberfunction
// SetDevCopy().
//
// a) SetDevCopy(0) (which is the default).
//    Only the pointers of the 'added' devices are stored.
//    This mode is typically used by making studies based on a fixed set
//    of devices which stays under user control or is kept on an external
//    file/tree. 
//    In this way the NcEvent just represents a 'logical structure' for the
//    physics analysis.
//
//    Note :
//    Modifications made to the original devices also affect the device
//    objects which are stored in the NcEvent. 
//
// b) SetDevCopy(1).
//    Of every 'added' device a private copy will be made of which the pointer
//    will be stored.
//    In this way the NcEvent represents an entity on its own and modifications
//    made to the original calorimeters do not affect the NcCalorimeter objects
//    which are stored in the NcEvent. 
//    This mode will allow 'adding' many different devices into an NcEvent by
//    creating only one device instance in the main programme and using the
//    Reset() and parameter setting memberfunctions of the object representing the device.
//
//    Note :
//    The copy is made using the Clone() memberfunction.
//    All devices (i.e. classes derived from TObject) have the default TObject::Clone() 
//    memberfunction.
//    However, devices generally contain an internal (signal) data structure
//    which may include pointers to other objects. Therefore it is recommended to provide
//    for all devices a specific copy constructor and override the default Clone()
//    memberfunction using this copy constructor.
//    Examples for this may be seen from NcCalorimeter, NcSignal and NcDevice.   
//
// See also the documentation provided for the memberfunction SetOwner(). 
//
// Coding example to make an event consisting of a primary vertex,
// 2 secondary vertices and a calorimeter.
// --------------------------------------------------------------
// vp contains the tracks 1,2,3 and 4 (primary vertex)
// v1 contains the tracks 5,6 and 7   (sec. vertex)
// v2 contains the jets 1 and 2       (sec. vertex)
//
//        NcEvent evt;
//
// Specify the event object as the repository of all objects
// for the event building and physics analysis.
// 
//        evt.SetDevCopy(1);
//        evt.SetTrackCopy(1);
//
// Fill the event structure with the basic objects
// 
//        NcCalorimeter emcal1;
//        NcCalorimeter emcal2;
//         ...
//         ... // code to fill the emcal1 and emcal2 calorimeter data
//         ...
//
//        evt.AddDevice(emcal1);
//        evt.AddDevice(emcal2);
//
//        // Assume NcTOF has been derived from NcDevice
//        NcTOF tof1;
//        NcTOF tof2;
//         ...
//         ... // code to fill the tof1 and tof2 data
//         ...
//
//        evt.AddDevice(tof1);
//        evt.AddDevice(tof2);
//
//        NcTrack* tx=new NcTrack();
//        for (Int_t i=0; i<10; i++)
//        {
//         ...
//         ... // code to fill the track data
//         ...
//         evt.AddTrack(tx);
//         tx->Reset(); 
//        }
//
//        if (tx)
//        {
//         delete tx;
//         tx=0;
//        }
//
// Order and investigate all the hits of all the TOF devices
//
//        TObjArray* hits=evt.GetHits("NcTOF");
//        TObjArray* orderedtofs=evt.SortHits(hits);
//        Int_t nhits=0;
//        if (orderedtofs) nhits=orderedtofs->GetEntries();
//        for (Int_t i=0; i<nhits; i++)
//        {
//         NcSignal* sx=(NcSignal*)orderedtofs->At(i);
//         if (sx) sx->Data();
//        }
//
// Order and investigate all the hits of all the calorimeter devices
//
//        TObjArray* hits=evt.GetHits("NcCalorimeter");
//        TObjArray* orderedcals=evt.SortHits(hits);
//        Int_t nhits=0;
//        if (orderedcals) nhits=orderedcals->GetEntries();
//        for (Int_t i=0; i<nhits; i++)
//        {
//         NcSignal* sx=(NcSignal*)orderedcals->At(i);
//         if (sx) sx->Data();
//        }
//
// Build the event structure (vertices, jets, ...) for physics analysis
// based on the basic objects from the event repository.
//
//        NcJet j1,j2;
//        for (Int_t i=0; i<evt.GetNtracks(); i++)
//        {
//         tx=evt.GetTrack(i);
//         ...
//         ... // code to fill the jet data
//         ...
//        }
//
//        NcVertex vp;
//        tx=evt.GetTrack(1);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(2);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(3);
//        vp.AddTrack(tx);
//        tx=evt.GetTrack(4);
//        vp.AddTrack(tx);
//
//        Float_t rp[3]={2.4,0.1,-8.5};
//        vp.SetPosition(rp,"car");
//
//        NcVertex v1;
//        tx=evt.GetTrack(5);
//        v1.AddTrack(tx);
//        tx=evt.GetTrack(6);
//        v1.AddTrack(tx);
//        tx=evt.GetTrack(7);
//        v1.AddTrack(tx);
//
//        Float_t r1[3]={1.6,-3.2,5.7};
//        v1.SetPosition(r1,"car");
//
//
//        NcVertex v2;
//        v2.SetJetCopy(1);
//        v2.AddJet(j1);
//        v2.AddJet(j2);
//
//        Float_t r2[3]={6.2,4.8,1.3};
//        v2.SetPosition(r2,"car");
//
// Specify the vertices v1 and v2 as secondary vertices of the primary
//
//        vp.SetVertexCopy(1);
//        vp.AddVertex(v1);
//        vp.AddVertex(v2);
//
// Enter the physics structures into the event
//        evt.SetVertexCopy(1);
//        evt.AddVertex(vp,0);
//
// The jets j1 and j2 are already available via sec. vertex v2,
// but can be made available also from the event itself if desired.
//        NcJet* jx;
//        jx=v2.GetJet(1);
//        evt.AddJet(jx,0); 
//        jx=v2.GetJet(2);
//        evt.AddJet(jx,0); 
// 
//        evt.Data("sph");
//        v1.ListAll();
//        v2.List("cyl");
//
//        Float_t etot=evt.GetEnergy();
//        Nc3Vector ptot=evt.Get3Momentum();
//        Float_t loc[3];
//        evt.GetPosition(loc,"sph");
//        NcPosition r=v1.GetPosition();
//        r.Data(); 
//        Int_t nt=v2.GetNtracks();
//        NcTrack* tv=v2.GetTrack(1); // Access track number 1 of Vertex v2
//
//        evt.List();
//
//        Int_t nv=evt.GetNvtx();
//        NcVertex* vx=evt.GetVertex(1); // Access primary vertex
//        Float_t e=vx->GetEnergy();
//
//        Float_t M=evt.GetInvmass(); 
//
// Reconstruct the event from scratch
//
//        evt.Reset();
//        evt.SetNvmax(25); // Increase initial no. of sec. vertices
//        ...
//        ... // code to create tracks etc... 
//        ...
//
// Note : By default all quantities are in meter, GeV, GeV/c or GeV/c**2
//        but the user can indicate the usage of a different scale for
//        the metric and/or energy-momentum units via the SetUnitScale()
//        and SetEscale() memberfunctions, respectively.
//        The actual metric and energy-momentum unit scales in use can be
//        obtained via the GetUnitScale() and GetEscale() memberfunctions.
//
//--- Author: Nick van Eijndhoven 27-may-2001 Utrecht University
//- Modified: NvE $Date$
///////////////////////////////////////////////////////////////////////////

#include "NcEvent.h"
#include "Riostream.h"
 
ClassImp(NcEvent) // Class implementation to enable ROOT I/O
 
NcEvent::NcEvent() : NcVertex(),NcTimestamp()
{
// Default constructor.
// All variables initialised to default values.
 fRun=0;
 fEvent=0;
 fWeight=1;
 fDevices=0;
 fDevCopy=0;
 fHits=0;
 fOrdered=0;
 fDisplay=0;
 fDevs=0;
}
///////////////////////////////////////////////////////////////////////////
NcEvent::NcEvent(Int_t n) : NcVertex(n),NcTimestamp()
{
// Create an event to hold initially a maximum of n tracks
// All variables initialised to default values
 if (n<=0)
 {
  cout << " *** This NcVertex initialisation was invoked via the NcEvent ctor." << endl;
 }
 fRun=0;
 fEvent=0;
 fWeight=1;
 fDevices=0;
 fDevCopy=0;
 fHits=0;
 fOrdered=0;
 fDisplay=0;
 fDevs=0;
}
///////////////////////////////////////////////////////////////////////////
NcEvent::~NcEvent()
{
// Default destructor
 if (fDevices)
 {
  delete fDevices;
  fDevices=0;
 }
 if (fHits)
 {
  delete fHits;
  fHits=0;
 }
 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }
 if (fDisplay)
 {
  delete fDisplay;
  fDisplay=0;
 }
 if (fDevs)
 {
  delete fDevs;
  fDevs=0;
 }
}
///////////////////////////////////////////////////////////////////////////
NcEvent::NcEvent(const NcEvent& evt) : NcVertex(evt),NcTimestamp(evt)
{
// Copy constructor.
 fRun=evt.fRun;
 fEvent=evt.fEvent;
 fWeight=evt.fWeight;
 fDevCopy=evt.fDevCopy;

 fHits=0;
 fOrdered=0;
 fDisplay=0;
 fDevs=0;

 fDevices=0;
 Int_t ndevs=evt.GetNdevices();
 if (ndevs)
 {
  fDevices=new TObjArray(ndevs);
  if (fDevCopy) fDevices->SetOwner();
  for (Int_t i=1; i<=ndevs; i++)
  {
   TObject* dev=evt.GetDevice(i);
   if (dev)
   {
    if (fDevCopy)
    {
     fDevices->Add(dev->Clone());
    }
    else
    {
     fDevices->Add(dev);
    }
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::Reset()
{
// Reset all variables to default values
// The max. number of tracks is set to the initial value again
// The max. number of vertices is set to the default value again
// The event weight is set to 1 again
// Note : The DevCopy mode is maintained as it was set by the user before.

 NcVertex::Reset();

 Set();
 fRun=0;
 fEvent=0;
 fWeight=1;

 if (fDevices)
 {
  delete fDevices;
  fDevices=0;
 }
 if (fHits)
 {
  delete fHits;
  fHits=0;
 }
 if (fOrdered)
 {
  delete fOrdered;
  fOrdered=0;
 }
 if (fDisplay)
 {
  delete fDisplay;
  fDisplay=0;
 }
 if (fDevs)
 {
  delete fDevs;
  fDevs=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetOwner(Bool_t own)
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
 if (fDevices) fDevices->SetOwner(own);
 fDevCopy=mode;

 NcVertex::SetOwner(own);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetDayTime(TTimeStamp& stamp)
{
// Set the date and time stamp for this event.
// An exact copy of the entered date/time stamp will be saved with an
// accuracy of 1 nanosecond.
//
// Note : Since the introduction of the more versatile class NcTimestamp
//        and the fact that NcEvent has now been derived from it,
//        this memberfunction has become obsolete.
//        It is recommended to use the corresponding NcTimestamp
//        functionality directly for NcEvent instances.
//        This memberfunction is only kept for backward compatibility.

 Set(stamp.GetDate(),stamp.GetTime(),0,kTRUE,0);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetDayTime(TDatime& stamp)
{
// Set the date and time stamp for this event.
// The entered date/time will be interpreted as being the local date/time
// and the accuracy is 1 second.
//
// This function with the TDatime argument is mainly kept for backward
// compatibility reasons.
// It is recommended to use the corresponding NcTimestamp functionality
// directly for NcEvent instances.

 Set(stamp.GetDate(),stamp.GetTime(),0,kFALSE,0);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetRunNumber(Int_t run)
{
// Set the run number for this event
 fRun=run;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetEventNumber(Int_t evt)
{
// Set the event number for this event
 fEvent=evt;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetWeight(Double_t weight)
{
// Set the weight for this event.
// By default the weight is set to 1 in the constructor.
 fWeight=weight;
}
///////////////////////////////////////////////////////////////////////////
TTimeStamp NcEvent::GetDayTime() const
{
// Provide the date and time stamp for this event
//
// Note : Since the introduction of the more versatile class NcTimestamp
//        and the fact that NcEvent has now been derived from it,
//        this memberfunction has become obsolete.
//        It is recommended to use the corresponding NcTimestamp
//        functionality directly for NcEvent instances.
//        This memberfunction is only kept for backward compatibility.

 return (TTimeStamp)(*this);
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetRunNumber() const
{
// Provide the run number for this event
 return fRun;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetEventNumber() const
{
// Provide the event number for this event
 return fEvent;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcEvent::GetWeight() const
{
// Provide the weight for this event
 return fWeight;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetProjectile(Int_t a,Int_t z,Double_t pnuc,Int_t id)
{
// Set the projectile A, Z, momentum per nucleon and user defined particle ID.
// If not explicitly specified by the user, the projectile particle ID is set
// to zero by default and will not be stored in the event structure.
// The projectile specifications will be stored in a device named "Beam"
// which is an instance of NcSignal.
// As such these data are easily retrievable from the event structure.
// However, for backward compatibility reasons the beam data can also be
// retrieved via memberfunctions like GetProjectileA() etc...

 Int_t newdev=0;
 
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 
 if (!beam)
 {
  beam=new NcSignal();
  beam->SetNameTitle("Beam","Beam and target specifications");
  newdev=1;
 }

 if (a || z)
 {
  beam->AddNamedSlot("Aproj");
  beam->SetSignal(a,"Aproj");
  beam->AddNamedSlot("Zproj");
  beam->SetSignal(z,"Zproj");
 }
 beam->AddNamedSlot("Pnucproj");
 beam->SetSignal(pnuc,"Pnucproj");
 if (id)
 {
  beam->AddNamedSlot("Idproj");
  beam->SetSignal(id,"Idproj");
 }

 if (newdev)
 {
  AddDevice(beam);
  if (fDevCopy) delete beam;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetProjectileA() const
{
// Provide the projectile A value.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Aproj"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetProjectileZ() const
{
// Provide the projectile Z value.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Zproj"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcEvent::GetProjectilePnuc() const
{
// Provide the projectile momentum value per nucleon.
 Double_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=beam->GetSignal("Pnucproj");
 return val;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetProjectileId() const
{
// Provide the user defined particle ID of the projectile.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Idproj"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetTarget(Int_t a,Int_t z,Double_t pnuc,Int_t id)
{
// Set the target A, Z, momentum per nucleon and user defined particle ID.
// If not explicitly specified by the user, the target particle ID is set
// to zero by default and will not be stored in the event structure.
// The target specifications will be stored in a device named "Beam"
// which is an instance of NcSignal.
// As such these data are easily retrievable from the event structure.
// However, for backward compatibility reasons the beam data can also be
// retrieved via memberfunctions like GetTargetA() etc...

 Int_t newdev=0;
 
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 
 if (!beam)
 {
  beam=new NcSignal();
  beam->SetNameTitle("Beam","Beam and target specifications");
  newdev=1;
 }

 if (a || z)
 {
  beam->AddNamedSlot("Atarg");
  beam->SetSignal(a,"Atarg");
  beam->AddNamedSlot("Ztarg");
  beam->SetSignal(z,"Ztarg");
 }
 beam->AddNamedSlot("Pnuctarg");
 beam->SetSignal(pnuc,"Pnuctarg");
 if (id)
 {
  beam->AddNamedSlot("Idtarg");
  beam->SetSignal(id,"Idtarg");
 }

 if (newdev)
 {
  AddDevice(beam);
  if (fDevCopy) delete beam;
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetTargetA() const
{
// Provide the target A value.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Atarg"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetTargetZ() const
{
// Provide the target Z value.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Ztarg"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcEvent::GetTargetPnuc() const
{
// Provide the target momentum value per nucleon.
 Double_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=beam->GetSignal("Pnuctarg");
 return val;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetTargetId() const
{
// Provide the user defined particle ID of the target.
 Int_t val=0;
 NcSignal* beam=(NcSignal*)GetDevice("Beam");
 if (beam) val=int(beam->GetSignal("Idtarg"));
 return val;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::HeaderData()
{
// Provide event header information
 const char* name=GetName();
 const char* title=GetTitle();
 cout << " *" << ClassName() << "::Data*";
 if (strlen(name))  cout << " Name : " << GetName();
 if (strlen(title)) cout << " Title : " << GetTitle();
 cout << endl;
 Date(1);
 cout << "  Run : " << fRun << " Event : " << fEvent
      << " Weight : " << fWeight << endl;
 ShowDevices(0);
 ShowTracks(0);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::Data(TString f,TString u)
{
// Provide event information within the coordinate frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 HeaderData();
 NcVertex::Data(f,u);
} 
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetNdevices() const
{
// Provide the number of stored devices
 Int_t ndevs=0;
 if (fDevices) ndevs=fDevices->GetEntries();
 return ndevs;
} 
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetNdevices(const char* classname) const
{
// Provide the number of stored devices of the specified class.

 Int_t ndevs=0;
 for (Int_t idev=1; idev<=GetNdevices(); idev++)
 {
  TObject* obj=GetDevice(idev);
  if (!obj) continue;

  if (obj->InheritsFrom(classname)) ndevs++;
 }
 return ndevs;
} 
///////////////////////////////////////////////////////////////////////////
void NcEvent::AddDevice(TObject& d)
{
// Add a device to the event.
//
// Note :
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All devices (i.e. classes derived from TObject) have the default TObject::Clone() 
// memberfunction.
// However, devices generally contain an internal (signal) data structure
// which may include pointers to other objects. Therefore it is recommended to provide
// for all devices a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from NcCalorimeter.   

 if (!fDevices)
 {
  fDevices=new TObjArray();
  if (fDevCopy) fDevices->SetOwner();
 }
 
 // Add the device to this event
 if (fDevCopy)
 {
  fDevices->Add(d.Clone());
 }
 else
 {
  fDevices->Add(&d);
 }
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::RemoveDevice(TObject* d)
{
// Remove the specified device from the event.

 if (!fDevices || !d) return;

 TObject* obj=fDevices->Remove(d);
 if (obj)
 {
  if (fDevCopy) delete obj;
  fDevices->Compress();
 }
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::SetDevCopy(Int_t j)
{
// (De)activate the creation of private copies of the added devices.
// j=0 ==> No private copies are made; pointers of original devices are stored.
// j=1 ==> Private copies of the devices are made and these pointers are stored.
//
//
// Notes :
//  In case a private copy is made, this is performed via the Clone() memberfunction.
//  All devices (i.e. classes derived from TObject) have the default TObject::Clone() 
//  memberfunction.
//  However, devices generally contain an internal (signal) data structure
//  which may include pointers to other objects. Therefore it is recommended to provide
//  for all devices a specific copy constructor and override the default Clone()
//  memberfunction using this copy constructor.
//  An example for this may be seen from NcCalorimeter.   
//
//  Once the storage contains pointer(s) to device(s) one cannot
//  change the DevCopy mode anymore.
//  To change the DevCopy mode for an existing NcEvent containing
//  devices one first has to invoke Reset().

 if (!fDevices)
 {
  if (j==0 || j==1)
  {
   fDevCopy=j;
  }
  else
  {
   cout << " *" << ClassName() << "::SetDevCopy* Invalid argument : " << j << endl;
  }
 }
 else
 {
  cout << " *" << ClassName() << "::SetDevCopy* Storage already contained devices."
       << "  ==> DevCopy mode not changed." << endl; 
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetDevCopy() const
{
// Provide value of the DevCopy mode.
// 0 ==> No private copies are made; pointers of original devices are stored.
// 1 ==> Private copies of the devices are made and these pointers are stored.
//
// Note :
// In case a private copy is made, this is performed via the Clone() memberfunction.
// All devices (i.e. classes derived from TObject) have the default TObject::Clone() 
// memberfunction.
// However, devices generally contain an internal (signal) data structure
// which may include pointers to other objects. Therefore it is recommended to provide
// for all devices a specific copy constructor and override the default Clone()
// memberfunction using this copy constructor.
// An example for this may be seen from NcCalorimeter.   

 return fDevCopy;
}
///////////////////////////////////////////////////////////////////////////
TObject* NcEvent::GetDevice(Int_t i) const
{
// Return the i-th device of this event.
// The first device corresponds to i=1.

 if (!fDevices)
 {
  return 0;
 }
 else
 {
  Int_t ndevs=GetNdevices();
  if (i<=0 || i>ndevs)
  {
   cout << " *" << ClassName() << "::GetDevice* Invalid argument i : " << i
        << " ndevs = " << ndevs << endl;
   return 0;
  }
  else
  {
   return fDevices->At(i-1);
  }
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* NcEvent::GetDevice(TString name) const
{
// Return the device with name tag "name"
 if (!fDevices)
 {
  return 0;
 }
 else
 {
  TString s;
  Int_t ndevs=GetNdevices();
  for (Int_t i=0; i<ndevs; i++)
  {
   TObject* dev=fDevices->At(i);
   if (dev)
   {
    s=dev->GetName();
    if (s == name) return dev;
   }
  }

  return 0; // No matching name found
 }
}
///////////////////////////////////////////////////////////////////////////
TObject* NcEvent::GetIdDevice(Int_t id,TObjArray* devs) const
{
// Return the device with identifier "id" from the specified array "devs".
// In case devs=0 (which is the default) all devices stored in the event
// structure will be evaluated.
// Note : In case of multiple occurrences of identifier "id", the first
//        encountered matching device will be returned.

 TObjArray* arr=devs;
 if (!arr) arr=fDevices;

 if (!arr || id<0) return 0;

 Int_t idx=0;
 for (Int_t i=0; i<arr->GetSize(); i++)
 {
  TObject* dev=arr->At(i);
  if (dev)
  {
   idx=dev->GetUniqueID();
   if (idx==id) return dev;
  }
 }
 return 0; // No matching id found
}
///////////////////////////////////////////////////////////////////////////
TObject* NcEvent::GetIdDevice(Int_t id,const char* classname) const
{
// Return the device with identifier "id" of the specified class.
// Note : In case of multiple occurrences of identifier "id", the first
//        encountered matching device will be returned.

 if (!fDevices || id<0) return 0;

 Int_t idx=0;
 for (Int_t i=0; i<GetNdevices(); i++)
 {
  TObject* dev=fDevices->At(i);
  if (dev)
  {
   idx=dev->GetUniqueID();
   if (idx==id && dev->InheritsFrom(classname)) return dev;
  }
 }
 return 0; // No matching id found for the specified class
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::ShowDevices(Int_t mode) const
{
// Provide an overview of the available devices.
// The argument mode determines the amount of information as follows :
// mode = 0 ==> Only printout of the number of devices
//        1 ==> Provide a listing with 1 line of info for each device
//
// The default is mode=1.
//
 Int_t ndevs=GetNdevices();
 if (ndevs)
 {
  if (!mode)
  {
   cout << " There are " << ndevs << " devices available." << endl; 
  }
  else
  {
   cout << " The following " << ndevs << " devices are available :" << endl; 
   Int_t nh=0,nw=0;
   for (Int_t i=1; i<=ndevs; i++)
   {
    TObject* dev=GetDevice(i);
    if (dev)
    {
     const char* name=dev->GetName();
     cout << " Device number : " << i;
     cout << " Class : " << dev->ClassName() << " Id : " << dev->GetUniqueID();
     if (strlen(name)) cout << " Name : " << name;
     if (dev->InheritsFrom("NcDevice"))
     {
      nh=((NcDevice*)dev)->GetNhits();
      if (nh) cout << " Nhits : " << nh;
     }
     if (dev->InheritsFrom("NcSignal"))
     {
      nw=((NcSignal*)dev)->GetNwaveforms();
      if (nw) cout << " Nwaveforms : " << nw;
     }
     cout << endl;
    }
   }
  }
 }
 else
 {
  cout << " No devices present for this event." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::ShowDevices(const char* classname,Int_t mode) const
{
// Provide an overview of the available devices of the specified class.
// The argument mode determines the amount of information as follows :
// mode = 0 ==> Only printout of the number of devices
//        1 ==> Provide a listing with 1 line of info for each device
//
// The default is mode=1.
//
 Int_t ndevs=GetNdevices();
 if (ndevs)
 {
  Int_t ndevs2=GetNdevices(classname);
  if (!mode || !ndevs2)
  {
   cout << " There are " << ndevs2 << " selected devices available." << endl; 
  }
  else
  {
   cout << " The following " << ndevs2 << " selected devices are available :" << endl; 
   Int_t nh=0,nw=0;
   for (Int_t i=1; i<=ndevs; i++)
   {
    TObject* dev=GetDevice(i);
    if (dev)
    {
     if (dev->InheritsFrom(classname))
     {
      const char* name=dev->GetName();
      cout << " Device number : " << i;
      cout << " Class : " << dev->ClassName() << " Id : " << dev->GetUniqueID();
      if (strlen(name)) cout << " Name : " << name;
      if (dev->InheritsFrom("NcDevice"))
      {
       nh=((NcDevice*)dev)->GetNhits();
       if (nh) cout << " Nhits : " << nh;
      }
      if (dev->InheritsFrom("NcSignal"))
      {
       nw=((NcSignal*)dev)->GetNwaveforms();
       if (nw) cout << " Nwaveforms : " << nw;
      }
      cout << endl;
     }
    }
   }
  }
 }
 else
 {
  cout << " No devices present for this event." << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::GetDevices(const char* classname)
{
// Provide the references to the various devices derived from the
// specified class.
 if (fDevs) fDevs->Clear();

 Int_t ndev=GetNdevices();
 for (Int_t idev=1; idev<=ndev; idev++)
 {
  TObject* obj=GetDevice(idev);
  if (!obj) continue;

  if (obj->InheritsFrom(classname))
  {
   if (!fDevs) fDevs=new TObjArray();
   fDevs->Add(obj);
  }
 }
 return fDevs;
}
///////////////////////////////////////////////////////////////////////////
Int_t NcEvent::GetNhits(const char* classname)
{
// Provide the number of hits registered to the specified device class.
// The specified device class has to be derived from NcDevice.
// It is possible to indicate with the argument "classname" a specific
// device instead of a whole class of devices. However, in such a case
// it is more efficient to use the GetDevice() memberfunction directly.
 LoadHits(classname);
 Int_t nhits=0;
 if (fHits) nhits=fHits->GetEntries();
 return nhits;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::GetHits(const char* classname)
{
// Provide the references to all the hits registered to the specified
// device class.
// The specified device class has to be derived from NcDevice.
// It is possible to indicate with the argument "classname" a specific
// device instead of a whole class of devices. However, in such a case
// it is more efficient to use the GetDevice() memberfunction directly.
 LoadHits(classname);
 return fHits;
}
///////////////////////////////////////////////////////////////////////////
NcSignal* NcEvent::GetIdHit(Int_t id,const char* classname)
{
// Return the hit with unique identifier "id" for the specified device class.
 if (id<0) return 0;

 Int_t nhits=GetNhits(classname);
 if (!nhits) return 0;

 NcSignal* sx=0;
 Int_t sid=0;
 for (Int_t i=0; i<nhits; i++)
 {
  sx=(NcSignal*)fHits->At(i);
  if (sx)
  {
   sid=sx->GetUniqueID();
   if (id==sid) return sx;
  }
 }
 return 0; // No matching id found
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::LoadHits(const char* classname)
{
// Load the references to the various hits registered to the specified
// device class.
// The specified device class has to be derived from NcDevice.
 if (fHits) fHits->Clear();

 Int_t ndev=GetNdevices();
 for (Int_t idev=1; idev<=ndev; idev++)
 {
  TObject* obj=GetDevice(idev);
  if (!obj) continue;

  if (obj->InheritsFrom(classname) && obj->InheritsFrom("NcDevice"))
  {
   NcDevice* dev=(NcDevice*)GetDevice(idev);
   Int_t nhits=dev->GetNhits();
   if (nhits)
   {
    if (!fHits) fHits=new TObjArray();
    for (Int_t ih=1; ih<=nhits; ih++)
    {
     NcSignal* sx=dev->GetHit(ih);
     if (sx) fHits->Add(sx);
    }
   }
  }
 }
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortHits(const char* classname,Int_t idx,Int_t mode,Int_t mcal)
{
// Order the references to the various hits registered to the specified
// device class. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the index of the signal slot to perform the sorting on.
// By default the slotindex will be 1.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.

 if (idx<=0 || abs(mode)!=1) return 0;

 LoadHits(classname);

 NcDevice dev;
 TObjArray* ordered=dev.SortHits(idx,mode,fHits,mcal);

 if (fHits)
 {
  delete fHits;
  fHits=0;
 } 
 if (ordered) fHits=new TObjArray(*ordered);
 return fHits;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortHits(const char* classname,TString name,Int_t mode,Int_t mcal)
{
// Order the references to the various hits registered to the specified
// device class. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the name of the signal slot to perform the sorting on.
// In case no matching slotname is found, the signal will be skipped.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.
 
 if (abs(mode)!=1) return 0;

 LoadHits(classname);

 NcDevice dev;
 TObjArray* ordered=dev.SortHits(name,mode,fHits,mcal);

 if (fHits)
 {
  delete fHits;
  fHits=0;
 } 
 if (ordered) fHits=new TObjArray(*ordered);
 return fHits;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::GetExtremes(const char* classname,Float_t& vmin,Float_t& vmax,Int_t idx,Int_t mode)
{
// Provide the min. and max. signal values of the various hits registered
// to the specified device class.
// The input argument "idx" denotes the index of the signal slots to be investigated.
// The default is idx=1;
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the process as specified
// by the  "mode" argument. The definition of this "mode" parameter corresponds to
// the description provided in the GetSignal memberfunction of class NcSignal.
// The default is mode=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.

 if (idx<=0) return;

 LoadHits(classname);

 NcDevice dev;
 dev.GetExtremes(vmin,vmax,idx,fHits,mode);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::GetExtremes(const char* classname,Float_t& vmin,Float_t& vmax,TString name,Int_t mode)
{
// Provide the min. and max. signal values of the various hits registered
// to the specified device class.
// The input argument "name" denotes the name of the signal slots to be investigated.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the process as specified
// by the  "mode" argument. The definition of this "mode" parameter corresponds to
// the description provided in the GetSignal memberfunction of class NcSignal.
// The default is mode=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.

 LoadHits(classname);

 NcDevice dev;
 dev.GetExtremes(vmin,vmax,name,fHits,mode);
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::DisplayHits(const char* classname,Int_t idx,Float_t scale,Int_t dp,Int_t mode,Int_t mcol)
{
// 3D color display of the various hits registered to the specified device class.
// The user can specify the index (default=1) of the signal slot to perform the display for.
// The marker size will indicate the absolute value of the signal (specified by the slotindex)
// as a percentage of the input argument "scale".
// In case scale<0 the maximum absolute signal value encountered in the hit array will be used
// to define the 100% scale. The default is scale=-1.
// In case dp=1 the owning device position will be used, otherwise the hit position will
// be used in the display. The default is dp=0.
// Via the "mcol" argument the user can specify the marker color (see TPolyMarker3D).
// The default is mcol=blue.
// Signals which were declared as "Dead" will not be displayed.
// The gain etc... corrected signals will be used to determine the marker size.
// The gain correction is performed according to "mode" argument. The definition of this
// "mode" parameter corresponds to the description provided in the GetSignal
// memberfunction of class NcSignal.
// The default is mode=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.
//
// Note :
// ------
// Before any display activity, a TCanvas and a TView have to be initiated
// first by the user like for instance
// 
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();

 if (idx<=0) return;

 LoadHits(classname);

 NcDevice* dev=new NcDevice();
 dev->DisplayHits(idx,scale,fHits,dp,mode,mcol);

 if (fDisplay)
 {
  delete fDisplay;
  fDisplay=0;
 }
 fDisplay=dev;
}
///////////////////////////////////////////////////////////////////////////
void NcEvent::DisplayHits(const char* classname,TString name,Float_t scale,Int_t dp,Int_t mode,Int_t mcol)
{
// 3D color display of the various hits registered to the specified device class.
// The user can specify the name of the signal slot to perform the display for.
// The marker size will indicate the absolute value of the signal (specified by the slotname)
// as a percentage of the input argument "scale".
// In case scale<0 the maximum absolute signal value encountered in the hit array will be used
// to define the 100% scale. The default is scale=-1.
// In case dp=1 the owning device position will be used, otherwise the hit position will
// be used in the display. The default is dp=0.
// The marker size will indicate the percentage of the maximum encountered value
// of the absolute value of the name-specified input signal slots.
// Via the "mcol" argument the user can specify the marker color (see TPolyMarker3D).
// The default is mcol=blue.
// Signals which were declared as "Dead" will not be displayed.
// The gain etc... corrected signals will be used to determine the marker size.
// The gain correction is performed according to "mode" argument. The definition of this
// "mode" parameter corresponds to the description provided in the GetSignal
// memberfunction of class NcSignal.
// The default is mode=1 (for backward compatibility reasons).
//
// For more extended functionality see class NcDevice.
//
// Note :
// ------
// Before any display activity, a TCanvas and a TView have to be initiated
// first by the user like for instance
// 
// TCanvas* c1=new TCanvas("c1","c1");
// TView* view=new TView(1);
// view->SetRange(-1000,-1000,-1000,1000,1000,1000);
// view->ShowAxis();

 LoadHits(classname);

 NcDevice* dev=new NcDevice();
 dev->DisplayHits(name,scale,fHits,dp,mode,mcol);

 if (fDisplay)
 {
  delete fDisplay;
  fDisplay=0;
 }
 fDisplay=dev;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortDevices(const char* classname,TString name,Int_t mode,Int_t mcal)
{
// Order the references to the various devices based on hit signals registered
// to the specified device class. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the name of the signal slot to perform the sorting on.
// In case no matching slotname is found, the signal will be skipped.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//

 TObjArray* ordered=SortHits(classname,name,mode,mcal);
 
 if (!ordered) return 0;

 TObjArray* devs=SortDevices(ordered,"*",0,mcal);
 return devs;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortDevices(const char* classname,Int_t idx,Int_t mode,Int_t mcal)
{
// Order the references to the various devices based on hit signals registered
// to the specified device class. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the index of the signal slot to perform the sorting on.
// By default the slotindex will be 1.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1) or ordering in increasing order (mode=1).
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//

 TObjArray* ordered=SortHits(classname,idx,mode,mcal);
 
 if (!ordered) return 0;

 TObjArray* devs=SortDevices(ordered,0,0,mcal);
 return devs;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortDevices(TObjArray* hits,TString name,Int_t mode,Int_t mcal)
{
// Order the references to the various devices based on hit signals contained
// in the input array. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the name of the signal slot to perform the sorting on.
// In case no matching slotname is found, the signal will be skipped.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1), ordering in increasing order (mode=1) or no ordering (mode=0).
// The latter option provides a means to quickly obtain an ordered devices list
// when the hits in the array were already ordered by the user. In this case
// the input argument "name" is irrelevant.
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//

 if (!hits) return 0;

 TObjArray* ordered=hits;
 NcDevice dev;
 if (mode) ordered=dev.SortHits(name,mode,hits,mcal);
 
 if (!ordered) return 0;

 if (fOrdered)
 {
  fOrdered->Clear();
 }
 else
 {
  fOrdered=new TObjArray();
 }

 Int_t nhits=ordered->GetEntries();
 Int_t exist=0;
 for (Int_t ih=0; ih<nhits; ih++)
 {
  NcSignal* sx=(NcSignal*)ordered->At(ih);
  if (!sx) continue;
  NcDevice* dx=sx->GetDevice();
  exist=0;
  for (Int_t id=0; id<fOrdered->GetEntries(); id++)
  {
   NcDevice* odx=(NcDevice*)fOrdered->At(id);
   if (dx==odx)
   {
    exist=1;
    break;
   }
  }
  if (!exist) fOrdered->Add(dx);
 }
 return fOrdered;
}
///////////////////////////////////////////////////////////////////////////
TObjArray* NcEvent::SortDevices(TObjArray* hits,Int_t idx,Int_t mode,Int_t mcal)
{
// Order the references to the various devices based on hit signals contained
// in the input array. The ordered array is returned as a TObjArray.
// A "hit" represents an abstract object which is derived from NcSignal.
// The user can specify the index of the signal slot to perform the sorting on.
// By default the slotindex will be 1.
// Via the "mode" argument the user can specify ordering in decreasing
// order (mode=-1), ordering in increasing order (mode=1) or no ordering (mode=0).
// The latter option provides a means to quickly obtain an ordered devices list
// when the hits in the array were already ordered by the user. In this case
// the input argument "idx" is irrelevant.
// The default is mode=-1.
// Signals which were declared as "Dead" will be rejected.
// The gain etc... corrected signals will be used in the ordering process as
// specified by the "mcal" argument. The definition of this "mcal" parameter
// corresponds to the signal correction mode described in the GetSignal
// memberfunction of class NcSignal.
// The default is mcal=1 (for backward compatibility reasons).
//

 if (!hits) return 0;

 TObjArray* ordered=hits;
 NcDevice dev;
 if (mode) ordered=dev.SortHits(idx,mode,hits,mcal);
 
 if (!ordered) return 0;

 if (fOrdered)
 {
  fOrdered->Clear();
 }
 else
 {
  fOrdered=new TObjArray();
 }

 Int_t nhits=ordered->GetEntries();
 Int_t exist=0;
 for (Int_t ih=0; ih<nhits; ih++)
 {
  NcSignal* sx=(NcSignal*)ordered->At(ih);
  if (!sx) continue;
  NcDevice* dx=sx->GetDevice();
  exist=0;
  for (Int_t id=0; id<fOrdered->GetEntries(); id++)
  {
   NcDevice* odx=(NcDevice*)fOrdered->At(id);
   if (dx==odx)
   {
    exist=1;
    break;
   }
  }
  if (!exist) fOrdered->Add(dx);
 }
 return fOrdered;
}
///////////////////////////////////////////////////////////////////////////
TObject* NcEvent::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.
// This feature allows to store either NcEvent objects or objects derived from
// NcEvent via some generic AddEvent memberfunction, provided these derived
// classes also have a proper Clone memberfunction. 

 NcEvent* evt=new NcEvent(*this);
 if (name)
 {
  if (strlen(name)) evt->SetName(name);
 }
 return evt;
}
///////////////////////////////////////////////////////////////////////////

