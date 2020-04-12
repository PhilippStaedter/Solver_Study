#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_saeidi1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*re1_k1;
    dwdx[1] = 1.0*re5_k5*s42;
    dwdx[2] = (re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)))*(-re14_K11*re14_n2*pow(s17, 2*re14_n2)/(s17*pow(pow(re14_K12, re14_n2) + pow(s17, re14_n2), 2)) + re14_K11*re14_n2*pow(s17, re14_n2)/(s17*(pow(re14_K12, re14_n2) + pow(s17, re14_n2)))) + (-re14_K8*re14_n1*pow(s17, 2*re14_n1)/(s17*pow(pow(re14_K9, re14_n1) + pow(s17, re14_n1), 2)) + re14_K8*re14_n1*pow(s17, re14_n1)/(s17*(pow(re14_K9, re14_n1) + pow(s17, re14_n1))))*(re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45);
    dwdx[3] = -re5_k6;
    dwdx[4] = 1.0*re3_Y2;
    dwdx[5] = 1.0*re4_k3*s4;
    dwdx[6] = -1.0*re1_Y1;
    dwdx[7] = 1.0*re2_k2;
    dwdx[8] = 1.0*re4_k3*s19;
    dwdx[9] = re8_Y3;
    dwdx[10] = -re4_k4;
    dwdx[11] = 1.0*re5_k5*s16;
    dwdx[12] = -re14_K7 - re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1));
}