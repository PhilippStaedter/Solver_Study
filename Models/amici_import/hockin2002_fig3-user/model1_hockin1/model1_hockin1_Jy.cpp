#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void Jy_model1_hockin1(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 8:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay8, 2)) + 0.5*pow(-my8 + y8, 2)/pow(sigmay8, 2);
            break;
        case 9:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay9, 2)) + 0.5*pow(-my9 + y9, 2)/pow(sigmay9, 2);
            break;
        case 10:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay10, 2)) + 0.5*pow(-my10 + y10, 2)/pow(sigmay10, 2);
            break;
        case 11:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay11, 2)) + 0.5*pow(-my11 + y11, 2)/pow(sigmay11, 2);
            break;
        case 12:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay12, 2)) + 0.5*pow(-my12 + y12, 2)/pow(sigmay12, 2);
            break;
        case 13:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay13, 2)) + 0.5*pow(-my13 + y13, 2)/pow(sigmay13, 2);
            break;
        case 14:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay14, 2)) + 0.5*pow(-my14 + y14, 2)/pow(sigmay14, 2);
            break;
        case 15:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay15, 2)) + 0.5*pow(-my15 + y15, 2)/pow(sigmay15, 2);
            break;
        case 16:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay16, 2)) + 0.5*pow(-my16 + y16, 2)/pow(sigmay16, 2);
            break;
        case 17:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay17, 2)) + 0.5*pow(-my17 + y17, 2)/pow(sigmay17, 2);
            break;
        case 18:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay18, 2)) + 0.5*pow(-my18 + y18, 2)/pow(sigmay18, 2);
            break;
        case 19:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay19, 2)) + 0.5*pow(-my19 + y19, 2)/pow(sigmay19, 2);
            break;
        case 20:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay20, 2)) + 0.5*pow(-my20 + y20, 2)/pow(sigmay20, 2);
            break;
        case 21:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay21, 2)) + 0.5*pow(-my21 + y21, 2)/pow(sigmay21, 2);
            break;
        case 22:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay22, 2)) + 0.5*pow(-my22 + y22, 2)/pow(sigmay22, 2);
            break;
        case 23:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay23, 2)) + 0.5*pow(-my23 + y23, 2)/pow(sigmay23, 2);
            break;
        case 24:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay24, 2)) + 0.5*pow(-my24 + y24, 2)/pow(sigmay24, 2);
            break;
        case 25:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay25, 2)) + 0.5*pow(-my25 + y25, 2)/pow(sigmay25, 2);
            break;
        case 26:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay26, 2)) + 0.5*pow(-my26 + y26, 2)/pow(sigmay26, 2);
            break;
        case 27:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay27, 2)) + 0.5*pow(-my27 + y27, 2)/pow(sigmay27, 2);
            break;
        case 28:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay28, 2)) + 0.5*pow(-my28 + y28, 2)/pow(sigmay28, 2);
            break;
        case 29:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay29, 2)) + 0.5*pow(-my29 + y29, 2)/pow(sigmay29, 2);
            break;
        case 30:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay30, 2)) + 0.5*pow(-my30 + y30, 2)/pow(sigmay30, 2);
            break;
        case 31:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay31, 2)) + 0.5*pow(-my31 + y31, 2)/pow(sigmay31, 2);
            break;
        case 32:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay32, 2)) + 0.5*pow(-my32 + y32, 2)/pow(sigmay32, 2);
            break;
        case 33:
            Jy[0] = 0.5*log(2*amici::pi*pow(sigmay33, 2)) + 0.5*pow(-my33 + y33, 2)/pow(sigmay33, 2);
            break;
    }
}