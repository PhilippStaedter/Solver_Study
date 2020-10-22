#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_model3_panteleev3(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = 0.5*(-2*my0 + 2*y0)/pow(sigmay0, 2);
            break;
        case 1:
            dJydy[0] = 0.5*(-2*my1 + 2*y1)/pow(sigmay1, 2);
            break;
        case 2:
            dJydy[0] = 0.5*(-2*my2 + 2*y2)/pow(sigmay2, 2);
            break;
        case 3:
            dJydy[0] = 0.5*(-2*my3 + 2*y3)/pow(sigmay3, 2);
            break;
        case 4:
            dJydy[0] = 0.5*(-2*my4 + 2*y4)/pow(sigmay4, 2);
            break;
        case 5:
            dJydy[0] = 0.5*(-2*my5 + 2*y5)/pow(sigmay5, 2);
            break;
        case 6:
            dJydy[0] = 0.5*(-2*my6 + 2*y6)/pow(sigmay6, 2);
            break;
        case 7:
            dJydy[0] = 0.5*(-2*my7 + 2*y7)/pow(sigmay7, 2);
            break;
    }
}