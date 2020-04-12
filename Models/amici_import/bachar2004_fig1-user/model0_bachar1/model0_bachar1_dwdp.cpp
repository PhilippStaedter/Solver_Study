#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_bachar1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[5] = I2;
            break;
        case 1:
            dwdp[1] = S*c*(I1 + I2*a)/(I1 + I2 + S);
            break;
        case 2:
            dwdp[1] = Beta*I2*S*c/(I1 + I2 + S);
            break;
        case 3:
            dwdp[0] = I1 + I2 + S;
            break;
        case 4:
            dwdp[1] = Beta*S*(I1 + I2*a)/(I1 + I2 + S);
            break;
        case 5:
            dwdp[2] = S;
            dwdp[3] = I1;
            dwdp[6] = I2;
            break;
        case 6:
            dwdp[4] = I1;
            break;
        case 7:
            dwdp[7] = I2;
            break;
    }
}