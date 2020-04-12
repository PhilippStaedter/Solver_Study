#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_ribeiro1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = sigma;
            break;
        case 1:
            dwdp[1] = Nplus;
            dwdp[4] = NN;
            break;
        case 2:
            dwdp[3] = NN;
            break;
        case 3:
            dwdp[0] = alpha;
            dwdp[2] = 1;
            break;
    }
}