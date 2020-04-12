#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_perelson2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = Tstarstar*muB;
            break;
        case 1:
            dwdp[7] = T*Ttot*r/pow(Tmax, 2);
            break;
        case 2:
            dwdp[1] = T*V;
            break;
        case 3:
            dwdp[2] = Tstar;
            break;
        case 4:
            dwdp[0] = N0*Tstarstar;
            dwdp[4] = Tstarstar;
            break;
        case 5:
            dwdp[5] = T;
            dwdp[6] = Tstar;
            break;
        case 6:
            dwdp[3] = V;
            break;
        case 7:
            dwdp[7] = T*(1 - Ttot/Tmax);
            break;
        case 8:
            dwdp[8] = 1;
            break;
    }
}