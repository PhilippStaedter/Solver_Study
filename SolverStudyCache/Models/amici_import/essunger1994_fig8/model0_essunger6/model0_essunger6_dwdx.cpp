#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_essunger6(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = V*kI;
    dwdx[1] = muA;
    dwdx[2] = V*kI*phi;
    dwdx[3] = amici_p;
    dwdx[4] = 2*amici_p;
    dwdx[5] = N0*muAA*(gamma*pow(t, 2)/(pow(Tc, 2) + pow(t, 2)) + 1);
    dwdx[6] = muAA;
    dwdx[7] = 2*rStar*(1 - Ttot/Tmax);
    dwdx[8] = rStar*(1 - Ttot/Tmax);
    dwdx[9] = kB;
    dwdx[10] = muM;
    dwdx[11] = kB;
    dwdx[12] = muM;
    dwdx[13] = -2*Tastarstar*rStar/Tmax;
    dwdx[14] = -Tastarstar*rStar/Tmax;
    dwdx[15] = kA;
    dwdx[16] = muV;
    dwdx[17] = Ta*kI;
    dwdx[18] = muVV;
    dwdx[19] = Ta*kI*phi;
    dwdx[20] = -1.0/3.0*s0*s1/(pow(V, 2.0/3.0)*pow(cbrt(V) + s1, 2));
}