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
// Class NcAttribObj
// Generic handling of detector signal (calibration) attributes.
//
// This class is meant to provide an NcAttrib object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// ---------
// NcAttrib a;
// a.SetGain(250.7);
// a.SetGain(1340,3);
// a.SetEdgeOn(3);
// a.SetOffset(-22.5,2);
// a.SetDead(1);
// a.Data();
//
// NcAttribObj b(a);
// b.Data();
//
// NcAttribObj c;
// c.Load(a);
// c.Data();
//
//--- Author: Nick van Eijndhoven 18-sep-2003 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcAttribObj.h"
#include "Riostream.h"
 
ClassImp(NcAttribObj) // Class implementation to enable ROOT I/O
 
NcAttribObj::NcAttribObj() : TObject(),NcAttrib()
{
// Creation of an NcAttrib object and initialisation of parameters.
// Several values of the same type (e.g. gain) can be stored in different slots.
// If needed, the storage for values will be expanded automatically
// when entering values.
}
///////////////////////////////////////////////////////////////////////////
NcAttribObj::NcAttribObj(NcAttrib& a) : TObject(),NcAttrib(a)
{
// Creation of an NcAttrib object and initialisation of parameters.
// All attributes are initialised to the values of the input NcAttrib.
}
///////////////////////////////////////////////////////////////////////////
NcAttribObj::~NcAttribObj()
{
// Destructor to delete dynamically allocated memory
}
///////////////////////////////////////////////////////////////////////////
NcAttribObj::NcAttribObj(const NcAttribObj& a) : TObject(a),NcAttrib(a)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
TObject* NcAttribObj::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 NcAttribObj* att=new NcAttribObj(*this);
 if (name)
 {
  if (strlen(name))
  {
   cout << " *" << ClassName() << "::Clone* No support for SetName." << endl;
  }
 }
 return att;
}
///////////////////////////////////////////////////////////////////////////
