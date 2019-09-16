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
// Class Nc4VectorObj
// Handling of Lorentz 4-vectors in various reference frames.
//
// This class is meant to provide an Nc4Vector object which is derived
// from TObject such that it can be stored in e.g. TObjArray etc...
// and that it can be written out using the ROOT I/O machinery.
//
// Example :
// =========
//
// Float_t a[4]={5,1,2,3};
// Float_t ea[4]={0.05,0.01,0.02,0.03};
// Float_t b[4]={10,4,5,6};
// Float_t eb[4]={0.1,0.04,0.05,0.06};
//
// Nc4Vector v,w;
//
// v.SetVector(a,"car");
// v.SetErrors(ea,"car");
// w.SetVector(b,"car");
// w.SetErrors(eb,"car");
//
// Nc4Vector add=v+w;
//
// Nc4Vector sub=v-w;
//
// Nc4VectorObj vec1(add);
//
// Nc4VectorObj vec2;
// vec2.Load(sub);
//
// vec1.Data();
// vec2.Data();
//
//--- Author: Nick van Eijndhoven 18-oct-1999 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "Nc4VectorObj.h"
 
ClassImp(Nc4VectorObj) // Class implementation to enable ROOT I/O
 
Nc4VectorObj::Nc4VectorObj() : TObject(),Nc4Vector()
{
// Default constructor
// Creation of an Nc4VectorObj object and initialisation of parameters.
// All attributes initialised to 0.
}
///////////////////////////////////////////////////////////////////////////
Nc4VectorObj::Nc4VectorObj(Nc4Vector& q) : TObject(),Nc4Vector(q)
{
// Creation of an Nc3VectorObj object and initialisation of parameters.
// All attributes are initialised to the values of the input Nc4Vector.
}
///////////////////////////////////////////////////////////////////////////
Nc4VectorObj::~Nc4VectorObj()
{
// Destructor to delete dynamically allocated memory.
}
///////////////////////////////////////////////////////////////////////////
Nc4VectorObj::Nc4VectorObj(Nc4VectorObj& v) : TObject(v),Nc4Vector(v)
{
// Copy constructor
}
///////////////////////////////////////////////////////////////////////////
