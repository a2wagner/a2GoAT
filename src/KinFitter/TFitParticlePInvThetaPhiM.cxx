// //Classname: TFitParticlePInvThetaPhiM
// Author: Verena Klose, Jan E. Sundermann, Thomas Goepfert (TU Dresden)

//________________________________________________________________
// 
// TFitParticlePInvThetaPhiM::
// --------------------
//
// Particle with 1/p, theta, phi, and d parametrization of the momentum 4vector
// (4 free parameters). The parametrization is chosen as follows:
//
// p_vec = (1/p, theta, phi)
// E(fit) =  Sqrt( |p|^2 + d^2 )
//

using namespace std;

#include <iostream>

#include "KinFitter/TFitParticlePInvThetaPhiM.h"

#include "TLorentzVector.h"
#include "TMath.h"

//ClassImp(TFitParticlePInvThetaPhiM)

//----------------
// Constructor --
//----------------
TFitParticlePInvThetaPhiM::TFitParticlePInvThetaPhiM()
  :TAbsFitParticle()
{ 
  init( 0, 0);
}

TFitParticlePInvThetaPhiM::TFitParticlePInvThetaPhiM( const TFitParticlePInvThetaPhiM& fitParticle )
  :TAbsFitParticle( fitParticle.GetName(), fitParticle.GetTitle() )
{

  _nPar = fitParticle._nPar;
  _u1 = fitParticle._u1;
  _u2 = fitParticle._u2;
  _u3 = fitParticle._u3;
  _covMatrix.ResizeTo(  fitParticle._covMatrix );
  _covMatrix = fitParticle._covMatrix;
  _iniparameters.ResizeTo( fitParticle._iniparameters );
  _iniparameters = fitParticle._iniparameters;
  _parameters.ResizeTo( fitParticle._parameters );
  _parameters = fitParticle._parameters;
  _pini = fitParticle._pini;
  _pcurr = fitParticle._pcurr;

}

TFitParticlePInvThetaPhiM::TFitParticlePInvThetaPhiM(TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle()
{ 
  init(pini, theCovMatrix);
}

TFitParticlePInvThetaPhiM::TFitParticlePInvThetaPhiM(const TString &name, const TString &title,
						     TLorentzVector* pini, const TMatrixD* theCovMatrix)
  :TAbsFitParticle(name, title)
{ 
  init(pini, theCovMatrix);
}

TAbsFitParticle* TFitParticlePInvThetaPhiM::clone( TString newname ) const {
  // Returns a copy of itself
  
  TAbsFitParticle* myclone = new TFitParticlePInvThetaPhiM( *this );
  if ( newname.Length() > 0 ) myclone->SetName(newname);
  return myclone;
  
}

//--------------
// Destructor --
//--------------
TFitParticlePInvThetaPhiM::~TFitParticlePInvThetaPhiM() {

}


//--------------
// Operations --
//--------------
void TFitParticlePInvThetaPhiM::init(TLorentzVector* pini, const TMatrixD* theCovMatrix ) {

  _nPar = 4;
  setIni4Vec(pini);
  setCovMatrix(theCovMatrix);
  
}
 
TLorentzVector* TFitParticlePInvThetaPhiM::calc4Vec( const TMatrixD* params ) {
  // Calculates a 4vector corresponding to the given
  // parameter values

  if (params == 0) {
    return 0;
  }

  if ( params->GetNcols() != 1 || params->GetNrows() !=_nPar ) {
    cout << "Parameter matrix has wrong size." << endl;
    return 0;
  }
  
  Double_t r = (*params)(0,0);
  Double_t theta = (*params)(1,0);
  Double_t phi = (*params)(2,0);
  Double_t d = (*params)(3,0);

  Double_t X = 1./r*TMath::Cos(phi)*TMath::Sin(theta);
  Double_t Y = 1./r*TMath::Sin(phi)*TMath::Sin(theta);
  Double_t Z = 1./r*TMath::Cos(theta);
  Double_t E = TMath::Sqrt( X*X + Y*Y + Z*Z + d*d );

  TLorentzVector* vec = new TLorentzVector( X, Y, Z, E );
  return vec;

}

void TFitParticlePInvThetaPhiM::setIni4Vec(const TLorentzVector* pini) {
  // Set the initial 4vector. Will also set the 
  // inital parameter values 

  if ( pini == 0 ) {

    _u1.SetXYZ(0., 0., 0.);
    _u3.SetXYZ(0., 0., 0.);
    _u2.SetXYZ(0., 0., 0.);
    _pini.SetXYZT(0., 0., 0., 0.);
    _pcurr = _pini;

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = 0.;
    _iniparameters(1,0) = 0.;
    _iniparameters(2,0) = 0.;
    _iniparameters(3,0) = 0.;

    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

  } else {

    _pini = (*pini);
    _pcurr = _pini;

    Double_t r = 1/_pini.P();
    Double_t theta = _pini.Theta();
    Double_t phi = _pini.Phi();
    Double_t m = _pini.M();

    _iniparameters.ResizeTo(_nPar,1);
    _iniparameters(0,0) = r;
    _iniparameters(1,0) = theta;
    _iniparameters(2,0) = phi;
    _iniparameters(3,0) = m;

    _parameters.ResizeTo(_nPar,1);
    _parameters = _iniparameters;

    _u1.SetXYZ( TMath::Cos(phi)*TMath::Sin(theta), TMath::Sin(phi)*TMath::Sin(theta), TMath::Cos(theta) );
    _u2.SetXYZ( TMath::Cos(phi)*TMath::Cos(theta), TMath::Sin(phi)*TMath::Cos(theta), -1.*TMath::Sin(theta) );
    _u3.SetXYZ( -1.*TMath::Sin(phi), TMath::Cos(phi), 0. );
  
  }

}

TMatrixD* TFitParticlePInvThetaPhiM::getDerivative() {
  // returns derivative dP/dy with P=(p,E) and y=(r, theta, phi, d) 
  // the free parameters of the fit. The columns of the matrix contain 
  // (dP/dr, dP/dtheta, ...).

  TMatrixD* DerivativeMatrix = new TMatrixD(4,4);
  (*DerivativeMatrix) *= 0.;

  Double_t r = _parameters(0,0);
  Double_t p = 1./r;
  Double_t theta = _parameters(1,0);
  Double_t phi = _parameters(2,0);
  Double_t m = _parameters(3,0);
  Double_t E = TMath::Sqrt( p*p + m*m );

  //1st column: dP/d(1/p)
  (*DerivativeMatrix)(0,0) = -1.*p*p*TMath::Cos(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(1,0) = -1.*p*p*TMath::Sin(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(2,0) = -1.*p*p*TMath::Cos(theta);
  (*DerivativeMatrix)(3,0) = -1. / (E*r*r*r);

  //2nd column: dP/dtheta
  (*DerivativeMatrix)(0,1) = p*TMath::Cos(phi)*TMath::Cos(theta);
  (*DerivativeMatrix)(1,1) = p*TMath::Sin(phi)*TMath::Cos(theta);
  (*DerivativeMatrix)(2,1) = -1.*p*TMath::Sin(theta);
  (*DerivativeMatrix)(3,1) = 0.;

   //3rd column: dP/dphi
  (*DerivativeMatrix)(0,2) = -1.*p*TMath::Sin(phi)*TMath::Sin(theta);
  (*DerivativeMatrix)(1,2) = p*TMath::Cos(phi)*TMath::Sin(theta);;
  (*DerivativeMatrix)(2,2) = 0.;
  (*DerivativeMatrix)(3,2) = 0.;

   //4rd column: dP/d(d)
  (*DerivativeMatrix)(0,3) = 0.;
  (*DerivativeMatrix)(1,3) = 0.;
  (*DerivativeMatrix)(2,3) = 0.;
  (*DerivativeMatrix)(3,3) = m / E;

  return DerivativeMatrix;  

}

TMatrixD* TFitParticlePInvThetaPhiM::transform(const TLorentzVector& vec) {
  // Returns the parameters corresponding to the given 
  // 4vector

  // retrieve parameters
  TMatrixD* tparams = new TMatrixD( _nPar, 1 );
  (*tparams)(0,0) = 1./vec.P();
  (*tparams)(1,0) = vec.Theta();
  (*tparams)(2,0) = vec.Phi();
  (*tparams)(3,0) = vec.M();

  return tparams;

}
