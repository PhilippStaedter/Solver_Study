#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bindschadler1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*k3*(-h1 + 1)/(R5 + c1);
    w[1] = 1.0*k3*(-h2 + 1)/(R5 + c2);
    w[2] = 1.0*Jleak_Cell1_Jleak;
    w[3] = 1.0*Jleak_Cell2_Jleak;
    w[4] = 1.0*Jpump_Cell1_Vp*pow(c1, 2)/(pow(Jpump_Cell1_Kp, 2) + pow(c1, 2));
    w[5] = 1.0*Jpump_Cell2_Vp*pow(c2, 2)/(pow(Jpump_Cell2_Kp, 2) + pow(c2, 2));
    w[6] = 1.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
    w[7] = 1.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
    w[8] = 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
    w[9] = 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
    w[10] = 1.0*diffusion_D*(-c1 + c2);
}