#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model3_naresh2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Q_0;
    w[1] = Delta_2*I2;
    w[2] = I2*Mu;
    w[3] = A*Alpha;
    w[4] = A*Mu;
    w[5] = Beta_1*I1*S/(A + I1 + I2 + S);
    w[6] = Beta_2*Epsilon*I2*S/(A + I1 + I2 + S);
    w[7] = Mu*S;
    w[8] = Beta_2*I2*S*(-Epsilon + 1)/(A + I1 + I2 + S);
    w[9] = I1*Theta;
    w[10] = Delta_1*I1;
    w[11] = I1*Mu;
    w[12] = I1*I2*amici_k;
}