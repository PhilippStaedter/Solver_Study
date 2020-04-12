#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 9.9999999999999998e-13*Activation_kf;
    dwdx[1] = -9.9999999999999998e-13*Binding_kb;
    dwdx[2] = -9.9999999999999998e-13*Activation_kb;
    dwdx[3] = 9.9999999999999998e-13*Binding_kf*CyclinA;
    dwdx[4] = 9.9999999999999998e-13*Binding_kf*Cdk2;
}