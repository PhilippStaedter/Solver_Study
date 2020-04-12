#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_karin5(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = -0.051000000000000004*beta*pow(alpha/g, 1.7)/(alpha*pow(pow(alpha/g, 1.7) + 1, 2)) - 0.051000000000000004*mbeta*pow(alpha/(amici_k*g), 1.7)/(alpha*pow(pow(alpha/(amici_k*g), 1.7) + 1, 2));
            break;
        case 1:
            dwdp[1] = 0.051000000000000004*mbeta*pow(alpha/(amici_k*g), 1.7)/(amici_k*pow(pow(alpha/(amici_k*g), 1.7) + 1, 2));
            dwdp[3] = (1.0/1440.0)*mbeta*(2.4414062500000001e-5*pow(amici_k, 7)*pow(g, 8)/pow((1.0/65536.0)*pow(amici_k, 8)*pow(g, 8) + 1, 2) + 6.3346658010681942*pow(1/(amici_k*g), 1.7)/(amici_k*pow(37.262740006283494*pow(1/(amici_k*g), 1.7) + 1, 2)) - 13668750.0/(pow(amici_k, 7)*pow(g, 6)*pow(1 + 11390625/(pow(amici_k, 6)*pow(g, 6)), 2)));
            break;
    }
}