#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_fung1_Fig3A_Vgly_0_5(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*alpha0 + 1.0*alpha2*pow(AcP/Kg2, n)/(pow(AcP/Kg2, n) + 1);
    w[1] = 1.0*alpha0 + 1.0*alpha1*pow(AcP/Kg1, n)/(pow(AcP/Kg1, n) + 1);
    w[2] = alpha0 + alpha3/(pow(LacI/Kg3, n) + 1);
    w[3] = 1.0*Acs*kd;
    w[4] = 1.0*LacI*kd;
    w[5] = 1.0*Pta*kd;
    w[6] = 1.0*C*(H*OAc - HOAc*Keq);
    w[7] = 1.0*AcP*kAck_f - 1.0*OAc*kAck_r;
    w[8] = 1.0*Acs*OAc*k2/(KM2 + OAc);
    w[9] = 1.0*AcCoA*Pta*k1/(AcCoA + KM1);
    w[10] = 1.0*AcCoA*kTCA;
    w[11] = 1.0*S0;
    w[12] = 1.0*k3*(HOAc - HOAc_E);
}