#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_li1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = TStar*beta;
            break;
        case 1:
            dwdp[5] = -T*r1*(T + TStar)/pow(Tmax, 2);
            dwdp[10] = -TStar*r2*(T + TStar)/pow(Tmax, 2);
            break;
        case 2:
            dwdp[3] = T;
            break;
        case 3:
            dwdp[1] = NN*TStar;
            dwdp[8] = TStar;
            break;
        case 4:
            dwdp[2] = V;
            break;
        case 5:
            dwdp[6] = T*V;
            dwdp[7] = T*V*sigma;
            break;
        case 6:
            dwdp[4] = T;
            dwdp[5] = T*(T + TStar)/Tmax;
            break;
        case 7:
            dwdp[9] = TStar;
            dwdp[10] = TStar*(T + TStar)/Tmax;
            break;
        case 8:
            dwdp[0] = 1;
            break;
        case 9:
            dwdp[7] = T*V*amici_k;
            break;
    }
}