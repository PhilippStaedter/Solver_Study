#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydsigmay_model2_sarma1(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 16:
            dJydsigmay[16] = 1.0/sigmay16 - 1.0*pow(-my16 + y16, 2)/pow(sigmay16, 3);
            break;
        case 17:
            dJydsigmay[17] = 1.0/sigmay17 - 1.0*pow(-my17 + y17, 2)/pow(sigmay17, 3);
            break;
        case 18:
            dJydsigmay[18] = 1.0/sigmay18 - 1.0*pow(-my18 + y18, 2)/pow(sigmay18, 3);
            break;
        case 19:
            dJydsigmay[19] = 1.0/sigmay19 - 1.0*pow(-my19 + y19, 2)/pow(sigmay19, 3);
            break;
        case 20:
            dJydsigmay[20] = 1.0/sigmay20 - 1.0*pow(-my20 + y20, 2)/pow(sigmay20, 3);
            break;
        case 21:
            dJydsigmay[21] = 1.0/sigmay21 - 1.0*pow(-my21 + y21, 2)/pow(sigmay21, 3);
            break;
        case 22:
            dJydsigmay[22] = 1.0/sigmay22 - 1.0*pow(-my22 + y22, 2)/pow(sigmay22, 3);
            break;
        case 23:
            dJydsigmay[23] = 1.0/sigmay23 - 1.0*pow(-my23 + y23, 2)/pow(sigmay23, 3);
            break;
        case 24:
            dJydsigmay[24] = 1.0/sigmay24 - 1.0*pow(-my24 + y24, 2)/pow(sigmay24, 3);
            break;
        case 25:
            dJydsigmay[25] = 1.0/sigmay25 - 1.0*pow(-my25 + y25, 2)/pow(sigmay25, 3);
            break;
        case 26:
            dJydsigmay[26] = 1.0/sigmay26 - 1.0*pow(-my26 + y26, 2)/pow(sigmay26, 3);
            break;
    }
}