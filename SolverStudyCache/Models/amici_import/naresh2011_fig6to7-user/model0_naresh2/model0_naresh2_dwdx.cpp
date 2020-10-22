#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_naresh2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = Alpha;
    dwdx[1] = Mu;
    dwdx[2] = -Beta_1*I1*S/pow(A + I1 + I2 + S, 2);
    dwdx[3] = -Beta_2*Epsilon*I2*S/pow(A + I1 + I2 + S, 2);
    dwdx[4] = -Beta_2*I2*S*(-Epsilon + 1)/pow(A + I1 + I2 + S, 2);
    dwdx[5] = -Beta_1*I1*S/pow(A + I1 + I2 + S, 2) + Beta_1*S/(A + I1 + I2 + S);
    dwdx[6] = -Beta_2*Epsilon*I2*S/pow(A + I1 + I2 + S, 2);
    dwdx[7] = -Beta_2*I2*S*(-Epsilon + 1)/pow(A + I1 + I2 + S, 2);
    dwdx[8] = Theta;
    dwdx[9] = Delta_1;
    dwdx[10] = Mu;
    dwdx[11] = I2*amici_k;
    dwdx[12] = Delta_2;
    dwdx[13] = Mu;
    dwdx[14] = -Beta_1*I1*S/pow(A + I1 + I2 + S, 2);
    dwdx[15] = -Beta_2*Epsilon*I2*S/pow(A + I1 + I2 + S, 2) + Beta_2*Epsilon*S/(A + I1 + I2 + S);
    dwdx[16] = -Beta_2*I2*S*(-Epsilon + 1)/pow(A + I1 + I2 + S, 2) + Beta_2*S*(-Epsilon + 1)/(A + I1 + I2 + S);
    dwdx[17] = I1*amici_k;
    dwdx[18] = -Beta_1*I1*S/pow(A + I1 + I2 + S, 2) + Beta_1*I1/(A + I1 + I2 + S);
    dwdx[19] = -Beta_2*Epsilon*I2*S/pow(A + I1 + I2 + S, 2) + Beta_2*Epsilon*I2/(A + I1 + I2 + S);
    dwdx[20] = Mu;
    dwdx[21] = -Beta_2*I2*S*(-Epsilon + 1)/pow(A + I1 + I2 + S, 2) + Beta_2*I2*(-Epsilon + 1)/(A + I1 + I2 + S);
}