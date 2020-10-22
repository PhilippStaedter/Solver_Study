#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_fung1_Fig3A_Vgly_0_5(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*AcCoA*Pta*k1/pow(AcCoA + KM1, 2) + 1.0*Pta*k1/(AcCoA + KM1);
    dwdx[1] = 1.0*kTCA;
    dwdx[2] = -1.0*alpha2*n*pow(AcP/Kg2, 2*n)/(AcP*pow(pow(AcP/Kg2, n) + 1, 2)) + 1.0*alpha2*n*pow(AcP/Kg2, n)/(AcP*(pow(AcP/Kg2, n) + 1));
    dwdx[3] = -1.0*alpha1*n*pow(AcP/Kg1, 2*n)/(AcP*pow(pow(AcP/Kg1, n) + 1, 2)) + 1.0*alpha1*n*pow(AcP/Kg1, n)/(AcP*(pow(AcP/Kg1, n) + 1));
    dwdx[4] = 1.0*kAck_f;
    dwdx[5] = 1.0*kd;
    dwdx[6] = 1.0*OAc*k2/(KM2 + OAc);
    dwdx[7] = -1.0*C*Keq;
    dwdx[8] = 1.0*k3;
    dwdx[9] = -1.0*k3;
    dwdx[10] = -alpha3*n*pow(LacI/Kg3, n)/(LacI*pow(pow(LacI/Kg3, n) + 1, 2));
    dwdx[11] = 1.0*kd;
    dwdx[12] = 1.0*C*H;
    dwdx[13] = -1.0*kAck_r;
    dwdx[14] = -1.0*Acs*OAc*k2/pow(KM2 + OAc, 2) + 1.0*Acs*k2/(KM2 + OAc);
    dwdx[15] = 1.0*kd;
    dwdx[16] = 1.0*AcCoA*k1/(AcCoA + KM1);
}