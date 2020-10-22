#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_bier2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[2] = 1;
            break;
        case 1:
            dwdp[5] = -T1 + T2;
            break;
        case 2:
            dwdp[0] = G1*T1;
            dwdp[1] = G2*T2;
            break;
        case 3:
            dwdp[3] = -T1*kp/pow(T1 + km, 2);
            dwdp[4] = -T2*kp/pow(T2 + km, 2);
            break;
        case 4:
            dwdp[3] = T1/(T1 + km);
            dwdp[4] = T2/(T2 + km);
            break;
    }
}