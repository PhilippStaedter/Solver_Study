#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_piedrafita1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*STU*k1;
    dwdx[1] = 1.0*k4;
    dwdx[2] = -1.0*STU*k3r;
    dwdx[3] = 1.0*SU*k5;
    dwdx[4] = 1.0*S*k1;
    dwdx[5] = -1.0*SU*k10r;
    dwdx[6] = -1.0*ST*k3r;
    dwdx[7] = 1.0*k4;
    dwdx[8] = -1.0*SU*k7r;
    dwdx[9] = -1.0*k1r;
    dwdx[10] = 1.0*T*k2;
    dwdx[11] = 1.0*U*k9;
    dwdx[12] = -1.0*k2r;
    dwdx[13] = 1.0*k3;
    dwdx[14] = 1.0*k10;
    dwdx[15] = -1.0*k9r;
    dwdx[16] = -1.0*STU*k10r;
    dwdx[17] = 1.0*ST*k5;
    dwdx[18] = -1.0*STU*k7r;
    dwdx[19] = 1.0*k4;
    dwdx[20] = -1.0*k5r;
    dwdx[21] = 1.0*U*k6;
    dwdx[22] = -1.0*k6r;
    dwdx[23] = 1.0*k7;
    dwdx[24] = 1.0*STUS*k2;
    dwdx[25] = 1.0*SUST*k6;
    dwdx[26] = 1.0*STUS*k9;
}