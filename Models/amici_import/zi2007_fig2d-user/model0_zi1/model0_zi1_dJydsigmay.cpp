#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_model0_zi1(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmay0 - 1.0*pow(-my0 + y0, 2)/pow(sigmay0, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmay1 - 1.0*pow(-my1 + y1, 2)/pow(sigmay1, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmay2 - 1.0*pow(-my2 + y2, 2)/pow(sigmay2, 3);
            break;
        case 3:
            dJydsigmay[3] = 1.0/sigmay3 - 1.0*pow(-my3 + y3, 2)/pow(sigmay3, 3);
            break;
        case 4:
            dJydsigmay[4] = 1.0/sigmay4 - 1.0*pow(-my4 + y4, 2)/pow(sigmay4, 3);
            break;
        case 5:
            dJydsigmay[5] = 1.0/sigmay5 - 1.0*pow(-my5 + y5, 2)/pow(sigmay5, 3);
            break;
        case 6:
            dJydsigmay[6] = 1.0/sigmay6 - 1.0*pow(-my6 + y6, 2)/pow(sigmay6, 3);
            break;
        case 7:
            dJydsigmay[7] = 1.0/sigmay7 - 1.0*pow(-my7 + y7, 2)/pow(sigmay7, 3);
            break;
        case 8:
            dJydsigmay[8] = 1.0/sigmay8 - 1.0*pow(-my8 + y8, 2)/pow(sigmay8, 3);
            break;
        case 9:
            dJydsigmay[9] = 1.0/sigmay9 - 1.0*pow(-my9 + y9, 2)/pow(sigmay9, 3);
            break;
        case 10:
            dJydsigmay[10] = 1.0/sigmay10 - 1.0*pow(-my10 + y10, 2)/pow(sigmay10, 3);
            break;
        case 11:
            dJydsigmay[11] = 1.0/sigmay11 - 1.0*pow(-my11 + y11, 2)/pow(sigmay11, 3);
            break;
        case 12:
            dJydsigmay[12] = 1.0/sigmay12 - 1.0*pow(-my12 + y12, 2)/pow(sigmay12, 3);
            break;
        case 13:
            dJydsigmay[13] = 1.0/sigmay13 - 1.0*pow(-my13 + y13, 2)/pow(sigmay13, 3);
            break;
        case 14:
            dJydsigmay[14] = 1.0/sigmay14 - 1.0*pow(-my14 + y14, 2)/pow(sigmay14, 3);
            break;
        case 15:
            dJydsigmay[15] = 1.0/sigmay15 - 1.0*pow(-my15 + y15, 2)/pow(sigmay15, 3);
            break;
    }
}