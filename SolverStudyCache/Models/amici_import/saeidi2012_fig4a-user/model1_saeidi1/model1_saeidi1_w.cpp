#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_saeidi1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = -1.0*re1_Y1*s2 + 1.0*re1_k1*s1;
    w[1] = (re14_K7 + re14_K8*pow(s17, re14_n1)/(pow(re14_K9, re14_n1) + pow(s17, re14_n1)))*(re14_K11*pow(s17, re14_n2)/(pow(re14_K12, re14_n2) + pow(s17, re14_n2)) + re14_k10 - s45);
    w[2] = 1.0*re2_k2*s2;
    w[3] = 1.0*re3_Y2*s19;
    w[4] = 1.0*re4_k3*s19*s4 - re4_k4*s42;
    w[5] = 1.0*re5_k5*s16*s42 - re5_k6*s17;
    w[6] = re8_Y3*s4;
}