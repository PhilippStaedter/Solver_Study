#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_essunger3(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = Tastarstar*muAA*(gamma*pow(t, 2)/(pow(Tc, 2) + pow(t, 2)) + 1);
            break;
        case 1:
            dwdp[0] = -2*N0*Tastarstar*Tc*gamma*muAA*pow(t, 2)/pow(pow(Tc, 2) + pow(t, 2), 2);
            break;
        case 2:
            dwdp[12] = Ta*Ttot*r/pow(Tmax, 2);
            dwdp[14] = 2*Ta*Ttot*r/pow(Tmax, 2);
            dwdp[15] = 2*Tastarstar*Ttot*rStar/pow(Tmax, 2);
            dwdp[16] = Tastarstar*Ttot*rStar/pow(Tmax, 2);
            break;
        case 3:
            dwdp[0] = N0*Tastarstar*muAA*pow(t, 2)/(pow(Tc, 2) + pow(t, 2));
            break;
        case 4:
            dwdp[4] = Tv;
            break;
        case 5:
            dwdp[2] = Tmstarstar;
            dwdp[3] = Tm;
            break;
        case 6:
            dwdp[1] = Ta*V;
            dwdp[11] = Ta*V*phi;
            break;
        case 7:
            dwdp[7] = Ta;
            break;
        case 8:
            dwdp[0] = N0*Tastarstar*(gamma*pow(t, 2)/(pow(Tc, 2) + pow(t, 2)) + 1);
            dwdp[6] = Tastarstar;
            break;
        case 9:
            dwdp[8] = Tmstarstar;
            dwdp[9] = Tm;
            break;
        case 10:
            dwdp[10] = Tv;
            break;
        case 11:
            dwdp[5] = V;
            break;
        case 12:
            dwdp[11] = Ta*V*kI;
            break;
        case 13:
            dwdp[12] = Ta*(1 - Ttot/Tmax);
            dwdp[14] = 2*Ta*(1 - Ttot/Tmax);
            break;
        case 14:
            dwdp[15] = 2*Tastarstar*(1 - Ttot/Tmax);
            dwdp[16] = Tastarstar*(1 - Ttot/Tmax);
            break;
        case 15:
            dwdp[13] = s1/(cbrt(V) + s1);
            break;
        case 16:
            dwdp[13] = -s0*s1/pow(cbrt(V) + s1, 2) + s0/(cbrt(V) + s1);
            break;
    }
}