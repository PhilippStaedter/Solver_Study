#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kim5(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[4] = 1.0*species_3;
            break;
        case 1:
            dwdp[5] = 1.0*species_4;
            break;
        case 3:
            dwdp[6] = 1.0*reaction_7_V*(-pow(reaction_7_Shalve, parameter_4)*log(reaction_7_Shalve) - pow(species_2, parameter_4)*log(species_2))/pow(pow(reaction_7_Shalve, parameter_4) + pow(species_2, parameter_4), 2);
            break;
        case 4:
            dwdp[7] = 1.0*reaction_8_V*pow(species_1, parameter_5)*log(species_1)/(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5)) + 1.0*reaction_8_V*pow(species_1, parameter_5)*(-pow(reaction_8_Shalve, parameter_5)*log(reaction_8_Shalve) - pow(species_1, parameter_5)*log(species_1))/pow(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5), 2);
            break;
        case 5:
            dwdp[8] = 1.0*species_4;
            break;
        case 6:
            dwdp[0] = 1.0*species_1;
            break;
        case 7:
            dwdp[1] = 1.0*species_2;
            break;
        case 8:
            dwdp[2] = 1.0*species_3;
            break;
        case 9:
            dwdp[3] = 1.0*species_4;
            break;
        case 10:
            dwdp[6] = -1.0*parameter_4*pow(reaction_7_Shalve, parameter_4)*reaction_7_V/(reaction_7_Shalve*pow(pow(reaction_7_Shalve, parameter_4) + pow(species_2, parameter_4), 2));
            break;
        case 11:
            dwdp[6] = 1.0/(pow(reaction_7_Shalve, parameter_4) + pow(species_2, parameter_4));
            break;
        case 12:
            dwdp[7] = 1.0*pow(species_1, parameter_5)/(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5));
            break;
        case 13:
            dwdp[7] = -1.0*parameter_5*pow(reaction_8_Shalve, parameter_5)*reaction_8_V*pow(species_1, parameter_5)/(reaction_8_Shalve*pow(pow(reaction_8_Shalve, parameter_5) + pow(species_1, parameter_5), 2));
            break;
    }
}