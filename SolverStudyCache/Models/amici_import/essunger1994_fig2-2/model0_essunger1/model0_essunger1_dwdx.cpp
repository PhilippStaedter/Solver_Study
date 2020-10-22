#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_essunger1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = V*ki;
    dwdx[1] = mua;
    dwdx[2] = V*ki*phi;
    dwdx[3] = 2*Ttot*r/Tmax;
    dwdx[4] = r;
    dwdx[5] = Ttot*r/Tmax;
    dwdx[6] = 2*r;
    dwdx[7] = NN*muastarstar;
    dwdx[8] = muastarstar;
    dwdx[9] = rstar;
    dwdx[10] = Ttot*rstar/Tmax;
    dwdx[11] = 2*rstar;
    dwdx[12] = 2*Ttot*rstar/Tmax;
    dwdx[13] = kb;
    dwdx[14] = mum;
    dwdx[15] = kb;
    dwdx[16] = mum;
    dwdx[17] = 2*Ta*r/Tmax;
    dwdx[18] = Tastarstar*rstar/Tmax;
    dwdx[19] = 2*Tastarstar*rstar/Tmax;
    dwdx[20] = Ta*r/Tmax;
    dwdx[21] = ka;
    dwdx[22] = muv;
    dwdx[23] = Ta*ki;
    dwdx[24] = muV;
    dwdx[25] = Ta*ki*phi;
}