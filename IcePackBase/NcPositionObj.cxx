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
// Class NcPositionObj
// Handling of positions in various reference frames.
//
// This class is meant to provide an NcPosition object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// =========
//
// Float_t a[3]={1,2,3};
// Float_t ea[3]={0.01,0.02,0.03};
// Float_t b[3]={4,5,6};
// Float_t eb[3]={0.04,0.05,0.06};
//
// NcPosition r1,r2;
//
// r1.SetPosition(a,"car");
// r1.SetPositionErrors(ea,"car");
// r2.SetPosition(b,"car");
// r2.SetPositionErrors(eb,"car");
//
// Nc3Vector sum=r1+r2;
// Nc3Vector rel=r1-r2;
//
// NcPositionObj rr1(r1);
// NcPositionObj rr2;
// rr2.Load(r2);
// NcPositionObj ssum(r1+r2);
//
// rr1.Data();
// rr2.Data();
// ssum.Data();
//
//--- Author: Nick van Eijndhoven 18-oct-1999 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcPositionObj.h"
 
ClassImp(NcPositionObj) // Class implementation to enable ROOT I/O
 
NcPositionObj::NcPositionObj() : TObject(),NcPosition()
{
// Creation of an NcPositionObj object and initialisation of parameters.
// All attributes initialised to 0.
}
///////////////////////////////////////////////////////////////////////////
NcPositionObj::NcPositionObj(NcPosition& p) : TObject(),NcPosition(p)
{
// Creation of an NcPositionObj object and initialisation of parameters.
// All attributes are initialised to the values of the input NcPosition.
}
///////////////////////////////////////////////////////////////////////////
NcPositionObj::~NcPositionObj()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
NcPositionObj::NcPositionObj(const NcPositionObj& p) : TObject(p),NcPosition(p)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
