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
// Class NcPosition
// Handling of positions in various reference frames.
//
// This class is meant to serve as a base class for ALICE objects
// that have a unique position in 3-dimensional space.
//
// Note :
// ------
// Positions (r), errors (e), reference frames (f) and angular units (u)
// are specified via
//
//    SetPosition(Float_t* r,TString f,TString u)
//    SetPositionErrors(Float_t* e,TString f,TString u)
//
// under the following conventions :
//
// f="car" ==> r in Cartesian coordinates   (x,y,z)
// f="sph" ==> r in Spherical coordinates   (r,theta,phi)
// f="cyl" ==> r in Cylindrical coordinates (rho,phi,z)
//
// u="rad" ==> angles in radians
// u="deg" ==> angles in degrees
//
// The "f" and "u" facilities only serve as a convenient user interface.
// Internally the actual storage of the various components is performed
// in a unique way. This allows setting/retrieval of vector components in a
// user selected frame/unit convention at any time. 
//
// The metric unit scale for the coordinates can be defined by the user
// via the SetUnitScale() memberfunction.
// This enables standardised expressions using numerical values of
// physical constants by means of the GetUnitScale() memberfunction.
// By default the unit scale is set to meter, corresponding to invokation
// of SetUnitScale(1).
// The user can specify a certain required metric unit scale in retreiving
// position components and/or distances.
// Please refer to the corresponding member functions for further details.
//   
//
// Example :
// ---------
//
// NcPosition q;
// Float_t pos[3]={-1,25,7};
// Float_t err[3]={0.08,1.85,0.5};
// q.SetPosition(pos,"car");
// q.SetPositionErrors(pos,"car");
// Float_t loc[3],dloc[3];
// q.GetPosition(loc,"sph","deg");
// q.GetPositionErrors(dloc,"sph","deg");
//
//--- Author: Nick van Eijndhoven 06-feb-1999 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcPosition.h"
#include "Riostream.h"
 
ClassImp(NcPosition) // Class implementation to enable ROOT I/O
 
NcPosition::NcPosition()
{
// Creation of an NcPosition object and initialisation of parameters.
// The unit scale for position coordinates is initialised to cm.
 fScale=1;
 fTstamp=0;
}
///////////////////////////////////////////////////////////////////////////
NcPosition::~NcPosition()
{
// Destructor to delete dynamically allocated memory
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
NcPosition::NcPosition(const NcPosition& p) : Nc3Vector(p)
{
// Copy constructor
 fScale=p.fScale;
 fTstamp=0;
 if (p.fTstamp) fTstamp=new NcTimestamp(*(p.fTstamp));
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPosition(Double_t* r,TString f,TString u)
{
// Store position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPosition(Double_t r1,Double_t r2,Double_t r3,TString f,TString u)
{
// Store position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetVector(r1,r2,r3,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::GetPosition(Double_t* r,TString f,TString u,Float_t scale) const
{
// Provide position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinates will be provided in the metric unit scale as
// stored in the NcPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinates in cm.

 Nc3Vector v=(Nc3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPosition(Float_t* r,TString f,TString u)
{
// Store position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::GetPosition(Float_t* r,TString f,TString u,Float_t scale) const
{
// Provide position according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinates will be provided in the metric unit scale as
// stored in the NcPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinates in cm.

 Nc3Vector v=(Nc3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetVector(r,f,u);
}
///////////////////////////////////////////////////////////////////////////
NcPosition& NcPosition::GetPosition()
{
// Provide position
 return (*this);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPosition(Nc3Vector& r)
{
// Set position
 Double_t a[3];
 r.GetVector(a,"sph");
 SetVector(a,"sph");
 r.GetErrors(a,"car");
 SetErrors(a,"car");
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPositionErrors(Double_t* e,TString f,TString u)
{
// Store position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetErrors(e,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPositionErrors(Double_t e1,Double_t e2,Double_t e3,TString f,TString u)
{
// Store position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetErrors(e1,e2,e3,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::GetPositionErrors(Double_t* e,TString f,TString u,Float_t scale) const
{
// Provide position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinate errors will be provided in the metric unit scale as
// stored in the NcPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinate errors in cm.

 Nc3Vector v=(Nc3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetErrors(e,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetPositionErrors(Float_t* e,TString f,TString u)
{
// Store position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".

 SetErrors(e,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::GetPositionErrors(Float_t* e,TString f,TString u,Float_t scale) const
{
// Provide position errors according to reference frame f
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The default is u="rad".
//
// By default the coordinate errors will be provided in the metric unit scale as
// stored in the NcPosition object.
// However, the user can select a different metric unit scale by
// specification of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the position coordinate errors in cm.

 Nc3Vector v=(Nc3Vector)(*this);
 if (scale>0) v*=fScale/scale;
 v.GetErrors(e,f,u);
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::ResetPosition()
{
// Reset the position and corresponding errors.
// When retrieved values of 0 will be returned.

 SetZero();
}
///////////////////////////////////////////////////////////////////////////
Double_t NcPosition::GetDistance(NcPosition& p,Float_t scale)
{
// Provide distance of the current NcPosition to position p.
// The error on the result can be obtained as usual by invoking
// GetResultError() afterwards. 
//
// By default the distance will be provided in the metric unit scale of
// the current NcPosition.
// This implies that the results of r1.GetDistance(r2) and r2.GetDistance(r1)
// may be numerically different in case r1 and r2 have different metric units.
// However, the user can specify a required metric unit scale by specification
// of the scale parameter.
// The convention is that scale=1 corresponds to meter, so specification
// of scale=0.01 will provide the distance in cm.
// As such it is possible to obtain a correctly computed distance even in case
// the position coordinates have a different unit scale.
// However, it is recommended to work always with one single unit scale.
//
 Nc3Vector d=(Nc3Vector)p;
 Float_t pscale=p.GetUnitScale();
 if ((pscale/fScale > 1.1) || (fScale/pscale > 1.1)) d=d*(pscale/fScale);
 Nc3Vector q=(Nc3Vector)(*this);
 d=d-q;
 Double_t dist=d.GetNorm();
 fDresult=d.GetResultError();

 if (scale>0)
 {
  dist*=fScale/scale;
  fDresult*=fScale/scale;
 }
 return dist;
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetUnitScale(Float_t s)
{
// Set the unit scale for the position coordinates.
// The scale is normalised w.r.t. the meter, so setting the unit scale
// to 0.01 means that all position coordinates are in cm.
// By default the unit scale is set to meter in the NcPosition constructor.
// It is recommended to use one single unit scale throughout a complete
// analysis and/or simulation project.
//
// Note : This memberfunction does not modify the numerical values of
//        the position coordinates.
//        It only specifies their numerical meaning.
// 
 if (s>0.)
 {
  fScale=s;
 }
 else
 {
  cout << " *NcPosition::SetUnitScale* Invalid argument s = " << s << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
Float_t NcPosition::GetUnitScale() const
{
// Provide the unit scale for the position coordinates.
// The scale is normalised w.r.t. the meter, so a unit scale of 0.01
// means that all position coordinates are in cm.
 return fScale;
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::SetTimestamp(NcTimestamp& t)
{
// Store the timestamp for this position.
 if (fTstamp) delete fTstamp;
 fTstamp=new NcTimestamp(t);
}
///////////////////////////////////////////////////////////////////////////
NcTimestamp* NcPosition::GetTimestamp()
{
// Provide the timestamp of this position.
 return fTstamp;
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::RemoveTimestamp()
{
// Remove the timestamp from this postion.
 if (fTstamp)
 {
  delete fTstamp;
  fTstamp=0;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcPosition::Data(TString f,TString u) const
{
// Provide all position/time information within the coordinate frame f.
//
// The string argument "u" allows to choose between different angular units
// in case e.g. a spherical frame is selected.
// u = "rad" : angles provided in radians
//     "deg" : angles provided in degrees
//
// The defaults are f="car" and u="rad".

 Nc3Vector::Data(f,u);
 cout << "   Metric unit : " << fScale << " meter" << endl;
 if (fTstamp) fTstamp->Date(1);
} 
///////////////////////////////////////////////////////////////////////////
