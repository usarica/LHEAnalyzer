/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOSPINZERO_7DCOMPLEX_WITHACCEP_GGH
#define ROOSPINZERO_7DCOMPLEX_WITHACCEP_GGH

#include "RooSpinZero.h"

class RooSpinZero_7DComplex_withAccep_ggH : public RooAbsPdf {

 public:

  struct measurables{
    RooRealVar* m1;
    RooRealVar* m2;
    RooRealVar* h1;
    RooRealVar* h2;
    RooRealVar* hs;
    RooRealVar* Phi;
    RooRealVar* Phi1;
  };
  
  struct modelParameters{
    RooRealVar* a1Val;
    RooRealVar* phi1Val;
    RooRealVar* a2Val;
    RooRealVar* phi2Val;
    RooRealVar* a3Val;
    RooRealVar* phi3Val;
    int parameterization;
    RooRealVar* g1Val;
    RooRealVar* g2Val;
    RooRealVar* g3Val;
    RooRealVar* g4Val;
    RooRealVar* g1ValIm;
    RooRealVar* g2ValIm;
    RooRealVar* g3ValIm;
    RooRealVar* g4ValIm;			
    RooRealVar* fL1;
    RooRealVar* fa2;
    RooRealVar* fa3;
    RooRealVar* phia2;
    RooRealVar* phia3;
    RooRealVar* mZ;
    RooRealVar* gamZ;
    RooRealVar* mX;
    RooRealVar* R1Val;
    RooRealVar* R2Val;
  };  	      
	      
  struct accepParameters{
    RooRealVar* aPhi;
    RooRealVar* bPhi;
    RooRealVar* cPhi;
    RooRealVar* dPhi;
    RooRealVar* ePhi;
    RooRealVar* aPhi1;
    RooRealVar* bPhi1;
    RooRealVar* cPhi1;
    RooRealVar* dPhi1;
    RooRealVar* ePhi1;
    RooRealVar* aH1;
    RooRealVar* bH1;
    RooRealVar* cH1;
    RooRealVar* dH1;
    RooRealVar* eH1;
    RooRealVar* aH2;
    RooRealVar* bH2;
    RooRealVar* cH2;
    RooRealVar* dH2;
    RooRealVar* eH2;
    RooRealVar* aHs;
    RooRealVar* bHs;
    RooRealVar* cHs;
    RooRealVar* dHs;
    RooRealVar* eHs;
    RooRealVar* aM1;
    RooRealVar* bM1;
    RooRealVar* cM1;
    RooRealVar* dM1;
    RooRealVar* aM2;
    RooRealVar* bM2;
    RooRealVar* cM2;
    RooRealVar* dM2;
  };

  RooSpinZero_7DComplex_withAccep_ggH() {} ; 
  RooSpinZero_7DComplex_withAccep_ggH(const char *name, const char *title,
				  measurables _measurables,
				  modelParameters _modelParams,
				  accepParameters _accepParams);

    RooSpinZero_7DComplex_withAccep_ggH(const RooSpinZero_7DComplex_withAccep_ggH& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const { return new RooSpinZero_7DComplex_withAccep_ggH(*this,newname); }
    inline virtual ~RooSpinZero_7DComplex_withAccep_ggH() { }
    
    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
    
protected:
    
    RooRealProxy m1 ;
    RooRealProxy m2 ;
    RooRealProxy h1 ;
    RooRealProxy h2 ;
    RooRealProxy hs ;
    RooRealProxy Phi ;
    RooRealProxy Phi1 ;

    RooRealProxy a1Val ;
    RooRealProxy phi1Val ;
    RooRealProxy a2Val ;
    RooRealProxy phi2Val ;
    RooRealProxy a3Val ;
    RooRealProxy phi3Val ;
    int parameterization ;
    RooRealProxy g1Val ;
    RooRealProxy g2Val ;
    RooRealProxy g3Val ;
    RooRealProxy g4Val ;
    RooRealProxy g1ValIm ;
    RooRealProxy g2ValIm ;
    RooRealProxy g3ValIm ;
    RooRealProxy g4ValIm ;
    RooRealProxy fL1;
    RooRealProxy fa2;
    RooRealProxy fa3;
    RooRealProxy phia2;
    RooRealProxy phia3;
    RooRealProxy mZ ;
    RooRealProxy gamZ ;
    RooRealProxy mX ;
    RooRealProxy R1Val ;
    RooRealProxy R2Val ;
    // acceptance parameters
    RooRealProxy aPhi ;
    RooRealProxy bPhi ;
    RooRealProxy cPhi ;
    RooRealProxy dPhi ;
    RooRealProxy ePhi ;
    RooRealProxy aPhi1 ;
    RooRealProxy bPhi1 ;
    RooRealProxy cPhi1 ;
    RooRealProxy dPhi1 ;
    RooRealProxy ePhi1 ;
    RooRealProxy aH1 ;
    RooRealProxy bH1 ;
    RooRealProxy cH1 ;
    RooRealProxy dH1 ;
    RooRealProxy eH1 ;
    RooRealProxy aH2 ;
    RooRealProxy bH2 ;
    RooRealProxy cH2 ;
    RooRealProxy dH2 ;
    RooRealProxy eH2 ;
    RooRealProxy aHs ;
    RooRealProxy bHs ;
    RooRealProxy cHs ;
    RooRealProxy dHs ;
    RooRealProxy eHs ;

    RooRealProxy aM1 ;
    RooRealProxy bM1 ;
    RooRealProxy cM1 ;
    RooRealProxy dM1 ;

    RooRealProxy aM2 ;
    RooRealProxy bM2 ;
    RooRealProxy cM2 ;
    RooRealProxy dM2 ;

    Double_t evaluate() const ;
    
private:
    
};

#endif