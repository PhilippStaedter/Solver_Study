#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model2_tripathi1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[9] = A;
            break;
        case 1:
            dwdp[2] = I1*S/(A + I1 + I2 + S);
            break;
        case 2:
            dwdp[2] = I2*S/(A + I1 + I2 + S);
            break;
        case 3:
            dwdp[5] = I1;
            dwdp[7] = I2;
            break;
        case 4:
            dwdp[0] = 1;
            break;
        case 5:
            dwdp[4] = I1;
            break;
        case 6:
            dwdp[1] = A;
            dwdp[3] = S;
            dwdp[6] = I1;
            dwdp[8] = I2;
            break;
    }
}