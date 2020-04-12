#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model2_band2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = VENUS*p2;
            break;
        case 1:
            dwdp[0] = -pow(VENUS, 2)*p2/pow(VENUS*p1_star + qj_star, 2);
            break;
        case 2:
            dwdp[0] = VENUS/(VENUS*p1_star + qj_star);
            dwdp[1] = VENUS*lambda_star;
            dwdp[2] = 1;
            break;
        case 3:
            dwdp[0] = -VENUS*p2/pow(VENUS*p1_star + qj_star, 2);
            break;
    }
}