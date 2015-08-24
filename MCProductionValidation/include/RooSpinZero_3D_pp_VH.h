/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOSPINZERO_3D_PP_VH
#define ROOSPINZERO_3D_PP_VH

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

class RooSpinZero_3D_pp_VH : public RooAbsPdf {
public:
    RooSpinZero_3D_pp_VH() {} ; 
    RooSpinZero_3D_pp_VH(const char *name, const char *title,
                      RooAbsReal& _h1,
                      RooAbsReal& _h2,
                      RooAbsReal& _Phi,
                      RooAbsReal& _m,
                      RooAbsReal& _Y,
                      RooAbsReal& _sqrts,
                      RooAbsReal& _mX,
                      RooAbsReal& _mZ,
                      RooAbsReal& _R1val,
                      RooAbsReal& _R2Val,
                      int _parameterizatiion, 
                      RooAbsReal& _a1Val,
                      RooAbsReal& _phi1Val,
                      RooAbsReal& _a2Val,
                      RooAbsReal& _phi2Val,
                      RooAbsReal& _a3Val,
                      RooAbsReal& _phi3Val,
                      RooAbsReal& _g1Val,
                      RooAbsReal& _g2Val,
                      RooAbsReal& _g3Val,
                      RooAbsReal& _g4Val,
                      RooAbsReal& _g1ValIm,
                      RooAbsReal& _g2ValIm,
                      RooAbsReal& _g3ValIm,
                      RooAbsReal& _g4ValIm,			
                      RooAbsReal& _fa2,
                      RooAbsReal& _fa3,
                      RooAbsReal& _phia2,
                      RooAbsReal& _phia3,
                      bool  _withAcc);
    
    RooSpinZero_3D_pp_VH(const RooSpinZero_3D_pp_VH& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new RooSpinZero_3D_pp_VH(*this,newname); }
    inline virtual ~RooSpinZero_3D_pp_VH() { }
    
    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
    
    
protected:
    
    RooRealProxy h1 ;
    RooRealProxy h2 ;
    RooRealProxy Phi ;
    RooRealProxy m ;    
    RooRealProxy Y ;        
    RooRealProxy sqrts ;
    RooRealProxy mX ;
    RooRealProxy mZ ;
    RooRealProxy R1Val ;
    RooRealProxy R2Val ;
    int parameterization ;
    RooRealProxy a1Val ;
    RooRealProxy phi1Val ;
    RooRealProxy a2Val ;
    RooRealProxy phi2Val ;
    RooRealProxy a3Val ;
    RooRealProxy phi3Val ;
    RooRealProxy g1Val ;
    RooRealProxy g2Val ;
    RooRealProxy g3Val ;
    RooRealProxy g4Val ;
    RooRealProxy g1ValIm ;
    RooRealProxy g2ValIm ;
    RooRealProxy g3ValIm ;
    RooRealProxy g4ValIm ;
    RooRealProxy fa2;
    RooRealProxy fa3;
    RooRealProxy phia2;
    RooRealProxy phia3;
    bool withAcc;
    
    
    Double_t evaluate() const ;
    
private:
    
    vector<TLorentzVector> Calculate4Momentum(float Mx,float M1,float M2,float theta,float theta1,float theta2,float _Phi1_,float _Phi_) const ;
    float partonicLuminosity(float mVal, float YVal, float sqrts) const ;
    
    ClassDef(RooSpinZero_3D_pp_VH,1) // Your description goes here...
};

#endif
