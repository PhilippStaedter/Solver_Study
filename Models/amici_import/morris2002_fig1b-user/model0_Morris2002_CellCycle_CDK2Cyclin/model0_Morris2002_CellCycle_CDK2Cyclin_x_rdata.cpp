#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = CDK2cycA;
    x_rdata[1] = CDK2cycA_star_;
    x_rdata[2] = Cdk2;
    x_rdata[3] = CyclinA;
}