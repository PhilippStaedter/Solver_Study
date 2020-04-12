#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_bindschadler1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*k3*(-h1 + 1)/pow(R5 + c1, 2);
    dwdx[1] = -2.0*Jpump_Cell1_Vp*pow(c1, 3)/pow(pow(Jpump_Cell1_Kp, 2) + pow(c1, 2), 2) + 2.0*Jpump_Cell1_Vp*c1/(pow(Jpump_Cell1_Kp, 2) + pow(c1, 2));
    dwdx[2] = 1.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)*(4*Jreceptor_Cell1_p*c1*r2/pow(R1 + c1, 2) - 4*Jreceptor_Cell1_p*r2/(R1 + c1) + 4*k1/pow(R3 + c1, 2))/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5)) - 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 5)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4)) + 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 3)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
    dwdx[3] = 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*r4/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1))) + 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)*(Open_to_Inactivated_Cell1_p*c1*r2/pow(R1 + c1, 2) - Open_to_Inactivated_Cell1_p*r2/(R1 + c1) + k1/pow(R3 + c1, 2))/((R1 + c1)*(R3 + c1)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2)) - 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/((R1 + c1)*pow(R3 + c1, 2)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1))) - 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/(pow(R1 + c1, 2)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1))) + 1.0*Open_to_Inactivated_Cell1_p*h1*r2*(c1*r4 + k2)/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
    dwdx[4] = -1.0*diffusion_D;
    dwdx[5] = -1.0*k3*(-h2 + 1)/pow(R5 + c2, 2);
    dwdx[6] = -2.0*Jpump_Cell2_Vp*pow(c2, 3)/pow(pow(Jpump_Cell2_Kp, 2) + pow(c2, 2), 2) + 2.0*Jpump_Cell2_Vp*c2/(pow(Jpump_Cell2_Kp, 2) + pow(c2, 2));
    dwdx[7] = 1.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)*(4*Jreceptor_Cell2_p*c2*r2/pow(R1 + c2, 2) - 4*Jreceptor_Cell2_p*r2/(R1 + c2) + 4*k1/pow(R3 + c2, 2))/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5)) - 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 5)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4)) + 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 3)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
    dwdx[8] = 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*r4/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2))) + 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)*(Open_to_Inactivated_Cell2_p*c2*r2/pow(R1 + c2, 2) - Open_to_Inactivated_Cell2_p*r2/(R1 + c2) + k1/pow(R3 + c2, 2))/((R1 + c2)*(R3 + c2)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2)) - 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/((R1 + c2)*pow(R3 + c2, 2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2))) - 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/(pow(R1 + c2, 2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2))) + 1.0*Open_to_Inactivated_Cell2_p*h2*r2*(c2*r4 + k2)/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
    dwdx[9] = 1.0*diffusion_D;
    dwdx[10] = -1.0*k3/(R5 + c1);
    dwdx[11] = 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 3)*pow(r2, 4)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
    dwdx[12] = 1.0*Open_to_Inactivated_Cell1_p*c1*r2*(c1*r4 + k2)/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
    dwdx[13] = -1.0*k3/(R5 + c2);
    dwdx[14] = 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 3)*pow(r2, 4)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
    dwdx[15] = 1.0*Open_to_Inactivated_Cell2_p*c2*r2*(c2*r4 + k2)/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
}