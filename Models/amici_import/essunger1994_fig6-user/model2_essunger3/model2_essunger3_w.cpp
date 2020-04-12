#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_essunger3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = N0*Tastarstar*muAA*(gamma*pow(t, 2)/(pow(Tc, 2) + pow(t, 2)) + 1);
    w[1] = Ta*V*kI;
    w[2] = Tmstarstar*kB;
    w[3] = Tm*kB;
    w[4] = Tv*kA;
    w[5] = V*muVV;
    w[6] = Tastarstar*muAA;
    w[7] = Ta*muA;
    w[8] = Tmstarstar*muM;
    w[9] = Tm*muM;
    w[10] = Tv*muV;
    w[11] = Ta*V*kI*phi;
    w[12] = Ta*r*(1 - Ttot/Tmax);
    w[13] = s0*s1/(cbrt(V) + s1);
    w[14] = 2*Ta*r*(1 - Ttot/Tmax);
    w[15] = 2*Tastarstar*rStar*(1 - Ttot/Tmax);
    w[16] = Tastarstar*rStar*(1 - Ttot/Tmax);
}