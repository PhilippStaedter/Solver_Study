%module model0_Morris2002_CellCycle_CDK2Cyclin
%import amici.i
// Add necessary symbols to generated header

%{
#include "wrapfunctions.h"
#include "amici/model_ode.h"
#include "amici/model_dae.h"
using namespace amici;
%}


// Process symbols in header
%include "wrapfunctions.h"
