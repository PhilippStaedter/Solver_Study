#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_alsheihk1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Q_0;
    w[1] = A*d;
    w[2] = S*(Beta_1*I1 + Beta_2*I2);
    w[3] = Mu*S;
    w[4] = I1*Theta;
    w[5] = I1*Mu;
    w[6] = Delta*I1;
    w[7] = I2*Mu;
    w[8] = Delta*I2;
    w[9] = A*Mu;
}