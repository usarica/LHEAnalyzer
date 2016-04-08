#ifdef _def_melatools_
#include <ZZMatrixElement/MELA/interface/TensorPdfFactory_HVV.h>
#else
#include "../include/TensorPdfFactory_HVV.h"
#endif


TensorPdfFactory_HVV::TensorPdfFactory_HVV(RooSpinTwo::modelMeasurables measurables_, RooSpin::VdecayType V1decay_, RooSpin::VdecayType V2decay_) :
TensorPdfFactory(measurables_, V1decay_, V2decay_)
{
  measurables.Y=0;
  makeParamsConst(true);
  initPDF();
}

TensorPdfFactory_HVV::~TensorPdfFactory_HVV(){
  destroyPDF();
}

void TensorPdfFactory_HVV::makeParamsConst(bool yesNo){
  couplings.Lambda->setConstant(true);

  parameters.mX->setConstant(yesNo);
  parameters.gamX->setConstant(yesNo);
  parameters.mV->setConstant(yesNo);
  parameters.gamV->setConstant(yesNo);
  parameters.R1Val->setConstant(yesNo);
  parameters.R2Val->setConstant(yesNo);
}

void TensorPdfFactory_HVV::initPDF(){
  PDF = new RooSpinTwo_7DComplex_HVV(
    "PDF", "PDF",
    measurables,
    parameters,
    couplings,
    V1decay, V2decay
    );
  PDF_base = (RooSpinTwo*)PDF;
}




