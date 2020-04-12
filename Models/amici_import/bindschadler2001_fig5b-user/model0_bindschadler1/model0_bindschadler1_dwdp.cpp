#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_bindschadler1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[6] = 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 5)*pow(c1, 5)*pow(h1, 4)*pow(r2, 5)/(pow(R1 + c1, 6)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5)) - 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 5)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
            dwdp[7] = 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 5)*pow(c2, 5)*pow(h2, 4)*pow(r2, 5)/(pow(R1 + c2, 6)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5)) - 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 5)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
            dwdp[8] = 1.0*pow(Open_to_Inactivated_Cell1_p, 2)*pow(c1, 2)*h1*pow(r2, 2)*(c1*r4 + k2)/(pow(R1 + c1, 3)*(R3 + c1)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2)) - 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/(pow(R1 + c1, 2)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            dwdp[9] = 1.0*pow(Open_to_Inactivated_Cell2_p, 2)*pow(c2, 2)*h2*pow(r2, 2)*(c2*r4 + k2)/(pow(R1 + c2, 3)*(R3 + c2)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2)) - 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/(pow(R1 + c2, 2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 1:
            dwdp[6] = 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*k1*pow(r2, 4)/(pow(R1 + c1, 4)*pow(R3 + c1, 2)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5));
            dwdp[7] = 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*k1*pow(r2, 4)/(pow(R1 + c2, 4)*pow(R3 + c2, 2)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5));
            dwdp[8] = 1.0*Open_to_Inactivated_Cell1_p*c1*h1*k1*r2*(c1*r4 + k2)/((R1 + c1)*pow(R3 + c1, 3)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2)) - 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/((R1 + c1)*pow(R3 + c1, 2)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            dwdp[9] = 1.0*Open_to_Inactivated_Cell2_p*c2*h2*k1*r2*(c2*r4 + k2)/((R1 + c2)*pow(R3 + c2, 3)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2)) - 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/((R1 + c2)*pow(R3 + c2, 2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 2:
            dwdp[0] = -1.0*k3*(-h1 + 1)/pow(R5 + c1, 2);
            dwdp[1] = -1.0*k3*(-h2 + 1)/pow(R5 + c2, 2);
            break;
        case 3:
            dwdp[6] = -4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 4)*(R3 + c1)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5));
            dwdp[7] = -4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 4)*(R3 + c2)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5));
            dwdp[8] = -1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2*(c1*r4 + k2)/((R1 + c1)*pow(R3 + c1, 2)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2));
            dwdp[9] = -1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2*(c2*r4 + k2)/((R1 + c2)*pow(R3 + c2, 2)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2));
            break;
        case 4:
            dwdp[8] = 1.0*Open_to_Inactivated_Cell1_p*c1*h1*r2/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            dwdp[9] = 1.0*Open_to_Inactivated_Cell2_p*c2*h2*r2/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 5:
            dwdp[0] = 1.0*(-h1 + 1)/(R5 + c1);
            dwdp[1] = 1.0*(-h2 + 1)/(R5 + c2);
            break;
        case 6:
            dwdp[6] = -4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 5)*pow(c1, 5)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 5)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5)) + 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 3)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
            dwdp[7] = -4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 5)*pow(c2, 5)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 5)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5)) + 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 3)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
            dwdp[8] = -1.0*pow(Open_to_Inactivated_Cell1_p, 2)*pow(c1, 2)*h1*r2*(c1*r4 + k2)/(pow(R1 + c1, 2)*(R3 + c1)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2)) + 1.0*Open_to_Inactivated_Cell1_p*c1*h1*(c1*r4 + k2)/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            dwdp[9] = -1.0*pow(Open_to_Inactivated_Cell2_p, 2)*pow(c2, 2)*h2*r2*(c2*r4 + k2)/(pow(R1 + c2, 2)*(R3 + c2)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2)) + 1.0*Open_to_Inactivated_Cell2_p*c2*h2*(c2*r4 + k2)/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 7:
            dwdp[8] = 1.0*Open_to_Inactivated_Cell1_p*pow(c1, 2)*h1*r2/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            dwdp[9] = 1.0*Open_to_Inactivated_Cell2_p*pow(c2, 2)*h2*r2/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 8:
            dwdp[2] = 1.0;
            break;
        case 9:
            dwdp[3] = 1.0;
            break;
        case 10:
            dwdp[4] = -2.0*Jpump_Cell1_Kp*Jpump_Cell1_Vp*pow(c1, 2)/pow(pow(Jpump_Cell1_Kp, 2) + pow(c1, 2), 2);
            break;
        case 11:
            dwdp[4] = 1.0*pow(c1, 2)/(pow(Jpump_Cell1_Kp, 2) + pow(c1, 2));
            break;
        case 12:
            dwdp[5] = -2.0*Jpump_Cell2_Kp*Jpump_Cell2_Vp*pow(c2, 2)/pow(pow(Jpump_Cell2_Kp, 2) + pow(c2, 2), 2);
            break;
        case 13:
            dwdp[5] = 1.0*pow(c2, 2)/(pow(Jpump_Cell2_Kp, 2) + pow(c2, 2));
            break;
        case 14:
            dwdp[6] = -4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 4)*pow(c1, 5)*pow(h1, 4)*pow(r2, 5)/(pow(R1 + c1, 5)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 5)) + 4.0*Jreceptor_Cell1_kf*pow(Jreceptor_Cell1_p, 3)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
            break;
        case 15:
            dwdp[6] = 1.0*pow(Jreceptor_Cell1_p, 4)*pow(c1, 4)*pow(h1, 4)*pow(r2, 4)/(pow(R1 + c1, 4)*pow(Jreceptor_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 4));
            break;
        case 16:
            dwdp[7] = -4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 4)*pow(c2, 5)*pow(h2, 4)*pow(r2, 5)/(pow(R1 + c2, 5)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 5)) + 4.0*Jreceptor_Cell2_kf*pow(Jreceptor_Cell2_p, 3)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
            break;
        case 17:
            dwdp[7] = 1.0*pow(Jreceptor_Cell2_p, 4)*pow(c2, 4)*pow(h2, 4)*pow(r2, 4)/(pow(R1 + c2, 4)*pow(Jreceptor_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 4));
            break;
        case 18:
            dwdp[8] = -1.0*Open_to_Inactivated_Cell1_p*pow(c1, 2)*h1*pow(r2, 2)*(c1*r4 + k2)/(pow(R1 + c1, 2)*(R3 + c1)*pow(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1), 2)) + 1.0*c1*h1*r2*(c1*r4 + k2)/((R1 + c1)*(R3 + c1)*(Open_to_Inactivated_Cell1_p*c1*r2/(R1 + c1) + k1/(R3 + c1)));
            break;
        case 19:
            dwdp[9] = -1.0*Open_to_Inactivated_Cell2_p*pow(c2, 2)*h2*pow(r2, 2)*(c2*r4 + k2)/(pow(R1 + c2, 2)*(R3 + c2)*pow(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2), 2)) + 1.0*c2*h2*r2*(c2*r4 + k2)/((R1 + c2)*(R3 + c2)*(Open_to_Inactivated_Cell2_p*c2*r2/(R1 + c2) + k1/(R3 + c2)));
            break;
        case 20:
            dwdp[10] = -1.0*c1 + 1.0*c2;
            break;
    }
}