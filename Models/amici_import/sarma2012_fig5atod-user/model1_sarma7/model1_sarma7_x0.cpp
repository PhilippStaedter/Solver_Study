#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model1_sarma7(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 999.99990368875297;
    x0[2] = 499.99995184437699;
    x0[3] = 99.999990368875203;
    x0[5] = 3999.99961475501;
    x0[8] = 999.99990368875297;
    x0[11] = 99.999990368875203;
    x0[12] = 499.99995184437699;
}