#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kowald1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 0.01*k10;
    dwdx[1] = 1.0*k2*species_0000002;
    dwdx[2] = 1.0*k3*(-species_0000002 + species_0000016);
    dwdx[3] = 1.0*k4*species_0000007;
    dwdx[4] = 1.0*k5*species_0000006;
    dwdx[5] = 0.01*k10;
    dwdx[6] = 1.0*k2*species_0000001;
    dwdx[7] = -1.0*k3*species_0000001;
    dwdx[8] = -1.0*k13a;
    dwdx[9] = 1.0*k13b;
    dwdx[10] = 1.0*k6*species_0000006;
    dwdx[11] = 1.0*k5*species_0000001;
    dwdx[12] = 1.0*k6*species_0000002;
    dwdx[13] = 1.0*k7*species_0000017;
    dwdx[14] = 1.0*k18;
    dwdx[15] = 2.0*k19*species_0000007;
    dwdx[16] = 1.0*k4*species_0000001;
    dwdx[17] = 1.0*k9;
    dwdx[18] = 1.0*k11;
    dwdx[19] = 1.0*k12;
    dwdx[20] = 1.0*k17;
    dwdx[21] = 1.0*k3*species_0000001;
    dwdx[22] = 1.0*k13a;
    dwdx[23] = 1.0*k7*species_0000006;
}