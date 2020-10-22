#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_fung1_Fig3A_Vgly_0_5(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[6] = 1.0*H*OAc - 1.0*HOAc*Keq;
            break;
        case 1:
            dwdp[6] = 1.0*C*OAc;
            break;
        case 2:
            dwdp[9] = -1.0*AcCoA*Pta*k1/pow(AcCoA + KM1, 2);
            break;
        case 3:
            dwdp[8] = -1.0*Acs*OAc*k2/pow(KM2 + OAc, 2);
            break;
        case 4:
            dwdp[6] = -1.0*C*HOAc;
            break;
        case 5:
            dwdp[1] = 1.0*alpha1*n*pow(AcP/Kg1, 2*n)/(Kg1*pow(pow(AcP/Kg1, n) + 1, 2)) - 1.0*alpha1*n*pow(AcP/Kg1, n)/(Kg1*(pow(AcP/Kg1, n) + 1));
            break;
        case 6:
            dwdp[0] = 1.0*alpha2*n*pow(AcP/Kg2, 2*n)/(Kg2*pow(pow(AcP/Kg2, n) + 1, 2)) - 1.0*alpha2*n*pow(AcP/Kg2, n)/(Kg2*(pow(AcP/Kg2, n) + 1));
            break;
        case 7:
            dwdp[2] = alpha3*n*pow(LacI/Kg3, n)/(Kg3*pow(pow(LacI/Kg3, n) + 1, 2));
            break;
        case 8:
            dwdp[11] = 1.0;
            break;
        case 9:
            dwdp[0] = 1.0;
            dwdp[1] = 1.0;
            dwdp[2] = 1;
            break;
        case 10:
            dwdp[1] = 1.0*pow(AcP/Kg1, n)/(pow(AcP/Kg1, n) + 1);
            break;
        case 11:
            dwdp[0] = 1.0*pow(AcP/Kg2, n)/(pow(AcP/Kg2, n) + 1);
            break;
        case 12:
            dwdp[2] = 1.0/(pow(LacI/Kg3, n) + 1);
            break;
        case 13:
            dwdp[9] = 1.0*AcCoA*Pta/(AcCoA + KM1);
            break;
        case 14:
            dwdp[8] = 1.0*Acs*OAc/(KM2 + OAc);
            break;
        case 15:
            dwdp[12] = 1.0*HOAc - 1.0*HOAc_E;
            break;
        case 16:
            dwdp[7] = 1.0*AcP;
            break;
        case 17:
            dwdp[7] = -1.0*OAc;
            break;
        case 18:
            dwdp[10] = 1.0*AcCoA;
            break;
        case 19:
            dwdp[3] = 1.0*Acs;
            dwdp[4] = 1.0*LacI;
            dwdp[5] = 1.0*Pta;
            break;
        case 20:
            dwdp[0] = -1.0*alpha2*pow(AcP/Kg2, 2*n)*log(AcP/Kg2)/pow(pow(AcP/Kg2, n) + 1, 2) + 1.0*alpha2*pow(AcP/Kg2, n)*log(AcP/Kg2)/(pow(AcP/Kg2, n) + 1);
            dwdp[1] = -1.0*alpha1*pow(AcP/Kg1, 2*n)*log(AcP/Kg1)/pow(pow(AcP/Kg1, n) + 1, 2) + 1.0*alpha1*pow(AcP/Kg1, n)*log(AcP/Kg1)/(pow(AcP/Kg1, n) + 1);
            dwdp[2] = -alpha3*pow(LacI/Kg3, n)*log(LacI/Kg3)/pow(pow(LacI/Kg3, n) + 1, 2);
            break;
    }
}