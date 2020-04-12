#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_brands1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = _J11_K11;
    dwdx[1] = _J8_K8;
    dwdx[2] = _J9_K9;
    dwdx[3] = _J10_K10*lys_R;
    dwdx[4] = _J2_K2;
    dwdx[5] = _J4_K4;
    dwdx[6] = _J5_K5;
    dwdx[7] = _J1_K1;
    dwdx[8] = _J3_K3;
    dwdx[9] = _J7_K7*lys_R;
    dwdx[10] = _J6_K6;
    dwdx[11] = Fru*_J10_K10;
    dwdx[12] = Glu*_J7_K7;
}