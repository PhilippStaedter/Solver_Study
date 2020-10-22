#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_Morris2002_CellCycle_CDK2Cyclin(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 1:
            dwdp[0] = -9.9999999999999998e-13*CDK2cycA_star_;
            break;
        case 2:
            dwdp[0] = 9.9999999999999998e-13*CDK2cycA;
            break;
        case 3:
            dwdp[1] = -9.9999999999999998e-13*CDK2cycA;
            break;
        case 4:
            dwdp[1] = 9.9999999999999998e-13*Cdk2*CyclinA;
            break;
    }
}