#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_essunger1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = Tastarstar*muastarstar;
            break;
        case 1:
            dwdp[13] = -2*Ta*Ttot*r/pow(Tmax, 2);
            dwdp[15] = -Tastarstar*Ttot*rstar/pow(Tmax, 2);
            dwdp[17] = -2*Tastarstar*Ttot*rstar/pow(Tmax, 2);
            dwdp[19] = -Ta*Ttot*r/pow(Tmax, 2);
            break;
        case 2:
            dwdp[4] = Tv;
            break;
        case 3:
            dwdp[2] = Tmstarstar;
            dwdp[3] = Tm;
            break;
        case 4:
            dwdp[1] = Ta*V;
            dwdp[11] = Ta*V*phi;
            break;
        case 5:
            dwdp[5] = V;
            break;
        case 6:
            dwdp[7] = Ta;
            break;
        case 7:
            dwdp[0] = NN*Tastarstar;
            dwdp[6] = Tastarstar;
            break;
        case 8:
            dwdp[8] = Tmstarstar;
            dwdp[9] = Tm;
            break;
        case 9:
            dwdp[10] = Tv;
            break;
        case 10:
            dwdp[11] = Ta*V*ki;
            break;
        case 11:
            dwdp[13] = 2*Ta*Ttot/Tmax;
            dwdp[18] = Ta;
            dwdp[19] = Ta*Ttot/Tmax;
            dwdp[20] = 2*Ta;
            break;
        case 12:
            dwdp[14] = Tastarstar;
            dwdp[15] = Tastarstar*Ttot/Tmax;
            dwdp[16] = 2*Tastarstar;
            dwdp[17] = 2*Tastarstar*Ttot/Tmax;
            break;
        case 13:
            dwdp[12] = 1;
            break;
    }
}