#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_tripathi1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Q_0;
    w[1] = A*d;
    w[2] = S*(Beta_1*I1 + Beta_2*I2)/(A + I1 + I2 + S);
    w[3] = S*d;
    w[4] = I1*Theta;
    w[5] = Delta*I1;
    w[6] = I1*d;
    w[7] = Delta*I2;
    w[8] = I2*d;
    w[9] = A*Alpha;
}