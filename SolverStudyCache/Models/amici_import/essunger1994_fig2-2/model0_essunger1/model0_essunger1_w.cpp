#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_essunger1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = NN*Tastarstar*muastarstar;
    w[1] = Ta*V*ki;
    w[2] = Tmstarstar*kb;
    w[3] = Tm*kb;
    w[4] = Tv*ka;
    w[5] = V*muV;
    w[6] = Tastarstar*muastarstar;
    w[7] = Ta*mua;
    w[8] = Tmstarstar*mum;
    w[9] = Tm*mum;
    w[10] = Tv*muv;
    w[11] = Ta*V*ki*phi;
    w[12] = s;
    w[13] = 2*Ta*Ttot*r/Tmax;
    w[14] = Tastarstar*rstar;
    w[15] = Tastarstar*Ttot*rstar/Tmax;
    w[16] = 2*Tastarstar*rstar;
    w[17] = 2*Tastarstar*Ttot*rstar/Tmax;
    w[18] = Ta*r;
    w[19] = Ta*Ttot*r/Tmax;
    w[20] = 2*Ta*r;
}