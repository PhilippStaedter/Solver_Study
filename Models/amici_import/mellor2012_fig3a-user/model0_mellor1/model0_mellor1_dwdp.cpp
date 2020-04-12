#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_mellor1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1000.0*parameter_2*species_1/pow(parameter_1 + species_1, 2);
            break;
        case 1:
            dwdp[0] = 1000.0*species_1/(parameter_1 + species_1);
            break;
        case 2:
            dwdp[1] = -1000.0*parameter_4*species_1/pow(parameter_3 + species_1, 2);
            break;
        case 3:
            dwdp[1] = 1000.0*species_1/(parameter_3 + species_1);
            break;
        case 4:
            dwdp[2] = -1000.0*parameter_6*species_1/pow(parameter_5 + species_1, 2);
            break;
        case 5:
            dwdp[2] = 1000.0*species_1/(parameter_5 + species_1);
            break;
        case 6:
            dwdp[3] = -1000.0*parameter_8*species_7/pow(parameter_7 + species_7, 2);
            dwdp[4] = -135.0*parameter_8*species_8/pow(parameter_7 + species_8, 2);
            break;
        case 7:
            dwdp[3] = 1000.0*species_7/(parameter_7 + species_7);
            dwdp[4] = 135.0*species_8/(parameter_7 + species_8);
            break;
    }
}