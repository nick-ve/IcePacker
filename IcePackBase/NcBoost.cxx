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
// Class NcBoost
// Perform various Lorentz transformations.
//
// Example :
// =========
//
// Float_t a[3]={0.1,0.2,0.3};
// Float_t ea[3]={0.01,0.02,0.03};
// Nc3Vector beta;
// beta.SetVector(a,"car");
// beta.SetErrors(ea,"car");
//
// NcBoost b1;
// b1.SetBeta(beta);
// b1.Data();
//
// Float_t b[4]={14,1,2,3};
// Float_t eb[4]={1.4,0.1,0.2,0.3};
// Nc4Vector p;
// p.SetVector(b,"car");
// p.SetErrors(eb,"car");
// Nc4Vector pprim=b1.Boost(p);
// p.Data();
// pprim.Data();
//
// p=b1.Inverse(pprim);
// pprim.Data();
// p.Data();
//
// Float_t c[4]={5,0,0,4};
// Float_t ec[4]={0.5,0,0,0.4};
// Nc4Vector q;
// q.SetVector(c,"car");
// q.SetErrors(ec,"car");
//
// NcBoost b2;
// b2.Set4Momentum(q);
// b2.Data("sph");
//
//--- Author: Nick van Eijndhoven 14-may-1996 Utrecht University
//- Modified: NvE $Date$ NCFS
///////////////////////////////////////////////////////////////////////////

#include "NcBoost.h"
#include "Riostream.h"
 
ClassImp(NcBoost) // Class implementation to enable ROOT I/O
 
NcBoost::NcBoost() : TObject()
{
// Creation of a Lorentz boost object and initialisation of parameters.
// Beta is set to (0,0,0) and consequently Gamma=1. 
// All errors are initialised to 0. 
 Double_t a[3]={0,0,0};
 fBeta.SetVector(a,"sph");
 fGamma=1;
 fDgamma=0;
 fDresult=0;
}
///////////////////////////////////////////////////////////////////////////
NcBoost::~NcBoost()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
NcBoost::NcBoost(const NcBoost& b) : TObject(b)
{
// Copy constructor
 fBeta=b.fBeta;
 fGamma=b.fGamma;
 fDgamma=b.fDgamma;
 fDresult=b.fDresult;
}
///////////////////////////////////////////////////////////////////////////
void NcBoost::SetBeta(Nc3Vector& b)
{
// Setting of boost parameters on basis of beta 3-vector.
// The errors on the beta 3-vector are taken from the input 3-vector.
// The gamma value and its error are calculated accordingly.
 fBeta=b;
 Double_t beta2=fBeta.Dot(fBeta);
 Double_t dbeta2=fBeta.GetResultError();

 if (beta2 > 1.)
 {
  cout << " *NcBoost::SetBeta* beta > 1." << endl;
 }
 fGamma=0;
 fDgamma=0;
 Double_t temp=1.-beta2;
 if (temp > 0.)
 {
  fGamma=sqrt(1./temp);
  fDgamma=fabs(dbeta2/(2.*pow(temp,1.5)));
 }
}
///////////////////////////////////////////////////////////////////////////
void NcBoost::Set4Momentum(Nc4Vector& p)
{
// Setting of boost parameters on basis of momentum 4-vector data.
// The errors of the input 4-vector are used to calculate the
// errors on the beta 3-vector and the gamma factor.
 Double_t E=p.GetScalar();
 Double_t dE=p.GetResultError();
 if (E <= 0.)
 {
  cout << " *NcBoost::Set4Momentum* Unphysical situation." << endl;
  p.Data();
 }
 else
 {
  Nc3Vector b=p.Get3Vector();
  Double_t vb[3],eb[3];
  b.GetVector(vb,"car");
  b.GetErrors(eb,"car");
  b=b/E;
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(eb[i]/E,2)+pow(vb[i]*dE/(E*E),2));
  }
  b.SetErrors(eb,"car");
  SetBeta(b);
 }
}
///////////////////////////////////////////////////////////////////////////
Nc3Vector NcBoost::GetBetaVector() const
{
// Provide the beta 3-vector.
 return fBeta;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcBoost::GetBeta()
{
// Provide the norm of the beta 3-vector.
// The error on the value can be obtained via GetResultError().
 Double_t norm=fBeta.GetNorm();
 fDresult=fBeta.GetResultError();
 return norm;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcBoost::GetGamma()
{
// Provide the gamma factor.
// The error on the value can be obtained via GetResultError().
 fDresult=fDgamma;
 return fGamma;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcBoost::GetResultError() const
{
// Provide the error on the result of an operation yielding a scalar.
// E.g. GetBeta() or GetGamma()
 return fDresult;
}
///////////////////////////////////////////////////////////////////////////
void NcBoost::Data(TString f)
{
// Printing of the boost parameter info in coordinate frame f.
 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();
 cout << " *NcBoost::Data* beta : " << beta << " error : " << dbeta
      << " gamma : " << fGamma << " error : " << fDgamma << endl;
 cout << "  Beta"; 
 fBeta.Data(f);
}
///////////////////////////////////////////////////////////////////////////
Nc4Vector NcBoost::Boost(Nc4Vector& v)
{
// Perform the Lorentz boost on the 4-vector v.
// Error propagation is performed automatically.
// Note : As an approximation Beta and p.Dot(Beta) are considered as
//        independent quantities.

 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();

 Double_t beta2=pow(beta,2);

 if (beta > 1.e-10)
 {
  Double_t E=v.GetScalar();
  Double_t dE=v.GetResultError();

  Nc3Vector p=v.Get3Vector();

  Double_t pdotbeta=p.Dot(fBeta);
  Double_t dpdotbeta=p.GetResultError();

  // Determine the new vector components
  Double_t Eprim=fGamma*(E-pdotbeta);

  Double_t z=((fGamma-1.)*pdotbeta/beta2)-fGamma*E;
  Nc3Vector add=fBeta*z;

  // Determine errors on the new vector components
  Double_t dEprim=sqrt(pow((E-pdotbeta)*fDgamma,2)+pow(fGamma*dE,2)
                      +pow(fGamma*dpdotbeta,2));
  Double_t dz=sqrt( pow(((fGamma-1.)/beta2)*dpdotbeta,2) + pow(fGamma*dE,2)
                   +pow((
    ((2./beta)-(4.*pow(beta,3)-6.*pow(beta,5))/(2.*pow((pow(beta,4)-pow(beta,6)),1.5)))*pdotbeta
    +beta*E/pow(fGamma,3))*dbeta,2) );

  Double_t vb[3],eb[3];
  fBeta.GetVector(vb,"car");
  fBeta.GetErrors(eb,"car");
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(z*eb[i],2)+pow(vb[i]*dz,2));
  }
  add.SetErrors(eb,"car");

  // Calculate the new 3-vector
  Nc3Vector pprim=p+add;

  // Set the components and errors of the new 4-vector 
  Nc4Vector w;
  w.SetVector(Eprim,pprim);
  w.SetScalarError(dEprim);

  return w;
 }
 else
 {
  return v;
 }
}
///////////////////////////////////////////////////////////////////////////
Nc4Vector NcBoost::Inverse(Nc4Vector& vprim)
{
// Perform the inverse Lorentz boost on the 4-vector vprim.
// Error propagation is performed automatically.
// Note : As an approximation Beta and pprim.Dot(Beta) are considered as
//        independent quantities.

 Double_t beta=fBeta.GetNorm();
 Double_t dbeta=fBeta.GetResultError();

 Double_t beta2=pow(beta,2);

 if (beta > 1.e-10)
 {
  Double_t Eprim=vprim.GetScalar();
  Double_t dEprim=vprim.GetResultError();

  Nc3Vector pprim=vprim.Get3Vector();

  Double_t pprimdotbeta=pprim.Dot(fBeta);
  Double_t dpprimdotbeta=pprim.GetResultError();

  // Determine the new vector components
  Double_t E=fGamma*(Eprim+pprimdotbeta);

  Double_t z=((fGamma-1.)*pprimdotbeta/beta2)+fGamma*Eprim;
  Nc3Vector add=fBeta*z;

  // Determine errors on the prime-vector components
  Double_t dE=sqrt(pow((Eprim+pprimdotbeta)*fDgamma,2)+pow(fGamma*dEprim,2)
                      +pow(fGamma*dpprimdotbeta,2));
  Double_t dz=sqrt( pow(((fGamma-1.)/beta2)*dpprimdotbeta,2) + pow(fGamma*dEprim,2)
                   +pow((
    ((2./beta)-(4.*pow(beta,3)-6.*pow(beta,5))/(2.*pow((pow(beta,4)-pow(beta,6)),1.5)))*pprimdotbeta
    -beta*Eprim/pow(fGamma,3))*dbeta,2) );

  Double_t vb[3],eb[3];
  fBeta.GetVector(vb,"car");
  fBeta.GetErrors(eb,"car");
  for (Int_t i=0; i<3; i++)
  {
   eb[i]=sqrt(pow(z*eb[i],2)+pow(vb[i]*dz,2));
  }
  add.SetErrors(eb,"car");

  // Calculate the new 3-vector
  Nc3Vector p=pprim+add;

  // Set the components and errors of the new 4-vector 
  Nc4Vector w;
  w.SetVector(E,p);
  w.SetScalarError(dE);

  return w;
 }
 else
 {
  return vprim;
 }
}
///////////////////////////////////////////////////////////////////////////
