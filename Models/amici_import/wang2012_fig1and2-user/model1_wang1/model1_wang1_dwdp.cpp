#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_wang1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[7] = TStar*beta;
            break;
        case 1:
            dwdp[3] = -T*r*(T + TStar)/pow(Tmax, 2);
            break;
        case 2:
            dwdp[1] = T;
            break;
        case 3:
            dwdp[6] = TStar;
            dwdp[7] = NN*TStar;
            break;
        case 4:
            dwdp[8] = V;
            break;
        case 5:
            dwdp[4] = T*V;
            dwdp[5] = T*V*sigma1;
            break;
        case 6:
            dwdp[2] = T;
            dwdp[3] = T*(T + TStar)/Tmax;
            break;
        case 7:
            dwdp[0] = 1;
            break;
        case 8:
            dwdp[5] = T*V*amici_k;
            break;
    }
}