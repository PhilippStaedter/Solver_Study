#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model2_naresh2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[3] = A;
            break;
        case 1:
            dwdp[5] = I1*S/(A + I1 + I2 + S);
            break;
        case 2:
            dwdp[6] = Epsilon*I2*S/(A + I1 + I2 + S);
            dwdp[8] = I2*S*(-Epsilon + 1)/(A + I1 + I2 + S);
            break;
        case 3:
            dwdp[10] = I1;
            break;
        case 4:
            dwdp[1] = I2;
            break;
        case 5:
            dwdp[6] = Beta_2*I2*S/(A + I1 + I2 + S);
            dwdp[8] = -Beta_2*I2*S/(A + I1 + I2 + S);
            break;
        case 6:
            dwdp[2] = I2;
            dwdp[4] = A;
            dwdp[7] = S;
            dwdp[11] = I1;
            break;
        case 7:
            dwdp[0] = 1;
            break;
        case 8:
            dwdp[9] = I1;
            break;
        case 9:
            dwdp[12] = I1*I2;
            break;
    }
}