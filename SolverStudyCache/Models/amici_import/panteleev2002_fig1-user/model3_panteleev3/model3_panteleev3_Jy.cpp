#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_model3_panteleev3(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay0, 2)) + 0.5*pow(-my0 + y0, 2)/pow(sigmay0, 2);
            break;
        case 1:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay1, 2)) + 0.5*pow(-my1 + y1, 2)/pow(sigmay1, 2);
            break;
        case 2:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay2, 2)) + 0.5*pow(-my2 + y2, 2)/pow(sigmay2, 2);
            break;
        case 3:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay3, 2)) + 0.5*pow(-my3 + y3, 2)/pow(sigmay3, 2);
            break;
        case 4:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay4, 2)) + 0.5*pow(-my4 + y4, 2)/pow(sigmay4, 2);
            break;
        case 5:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay5, 2)) + 0.5*pow(-my5 + y5, 2)/pow(sigmay5, 2);
            break;
        case 6:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay6, 2)) + 0.5*pow(-my6 + y6, 2)/pow(sigmay6, 2);
            break;
        case 7:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay7, 2)) + 0.5*pow(-my7 + y7, 2)/pow(sigmay7, 2);
            break;
    }
}