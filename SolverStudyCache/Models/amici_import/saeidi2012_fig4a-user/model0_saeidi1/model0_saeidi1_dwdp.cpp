#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_saeidi1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*s1;
            break;
        case 1:
            dwdp[0] = -1.0*s2;
            break;
        case 2:
            dwdp[1] = re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45;
            break;
        case 3:
            dwdp[1] = pow(s17, re14_n1)*(re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1));
            break;
        case 4:
            dwdp[1] = -re14_K8*pow(re14_K9, re14_n1)*re14_n1*pow(s17, re14_n1)*(re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45)/(re14_K9*pow(pow(re14_K9, re14_n1) + pow(s17, re14_n1), 2));
            break;
        case 5:
            dwdp[1] = re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1));
            break;
        case 6:
            dwdp[1] = pow(s17, re14_n2)*(re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)))/(pow(re14_K12, re14_n2) + pow(s17, re14_n2));
            break;
        case 7:
            dwdp[1] = -re14_K11*pow(re14_K12, re14_n2)*re14_n2*pow(s17, re14_n2)*(re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)))/(re14_K12*pow(pow(re14_K12, re14_n2) + pow(s17, re14_n2), 2));
            break;
        case 8:
            dwdp[1] = (re14_K8*pow(s17, re14_n1)*log(s17)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)) + re14_K8*pow(s17, re14_n1)*(-pow(re14_K9, re14_n1)*log(re14_K9) - pow(s17, re14_n1)*log(s17))/pow(pow(re14_K9, re14_n1) + pow(s17, re14_n1), 2))*(re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45);
            break;
        case 9:
            dwdp[1] = (re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)))*(re14_K11*pow(s17, re14_n2)*log(s17)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_K11*pow(s17, re14_n2)*(-pow(re14_K12, re14_n2)*log(re14_K12) - pow(s17, re14_n2)*log(s17))/pow(pow(re14_K12, re14_n2) + pow(s17, re14_n2), 2));
            break;
        case 10:
            dwdp[2] = 1.0*s2;
            break;
        case 11:
            dwdp[3] = 1.0*s19;
            break;
        case 12:
            dwdp[4] = 1.0*s19*s4;
            break;
        case 13:
            dwdp[4] = -s42;
            break;
        case 14:
            dwdp[5] = 1.0*s16*s42;
            break;
        case 15:
            dwdp[5] = -s17;
            break;
        case 16:
            dwdp[6] = s4;
            break;
    }
}