#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = -9.9999999999999998e-13*Activation_kb*CDK2cycA_star_ + 9.9999999999999998e-13*Activation_kf*CDK2cycA;
    w[1] = -9.9999999999999998e-13*Binding_kb*CDK2cycA + 9.9999999999999998e-13*Binding_kf*Cdk2*CyclinA;
}