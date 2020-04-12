#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

void dJydy_Pathak2013a(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
        case 8:
            dJydy[0] = 0.5*(-2*my8 + 2*y8)/pow(sigmay8, 2);
            break;
        case 9:
            dJydy[0] = 0.5*(-2*my9 + 2*y9)/pow(sigmay9, 2);
            break;
        case 10:
            dJydy[0] = 0.5*(-2*my10 + 2*y10)/pow(sigmay10, 2);
            break;
        case 11:
            dJydy[0] = 0.5*(-2*my11 + 2*y11)/pow(sigmay11, 2);
            break;
        case 12:
            dJydy[0] = 0.5*(-2*my12 + 2*y12)/pow(sigmay12, 2);
            break;
        case 13:
            dJydy[0] = 0.5*(-2*my13 + 2*y13)/pow(sigmay13, 2);
            break;
        case 14:
            dJydy[0] = 0.5*(-2*my14 + 2*y14)/pow(sigmay14, 2);
            break;
        case 15:
            dJydy[0] = 0.5*(-2*my15 + 2*y15)/pow(sigmay15, 2);
            break;
        case 16:
            dJydy[0] = 0.5*(-2*my16 + 2*y16)/pow(sigmay16, 2);
            break;
        case 17:
            dJydy[0] = 0.5*(-2*my17 + 2*y17)/pow(sigmay17, 2);
            break;
        case 18:
            dJydy[0] = 0.5*(-2*my18 + 2*y18)/pow(sigmay18, 2);
            break;
        case 19:
            dJydy[0] = 0.5*(-2*my19 + 2*y19)/pow(sigmay19, 2);
            break;
        case 20:
            dJydy[0] = 0.5*(-2*my20 + 2*y20)/pow(sigmay20, 2);
            break;
        case 21:
            dJydy[0] = 0.5*(-2*my21 + 2*y21)/pow(sigmay21, 2);
            break;
        case 22:
            dJydy[0] = 0.5*(-2*my22 + 2*y22)/pow(sigmay22, 2);
            break;
        case 23:
            dJydy[0] = 0.5*(-2*my23 + 2*y23)/pow(sigmay23, 2);
            break;
        case 24:
            dJydy[0] = 0.5*(-2*my24 + 2*y24)/pow(sigmay24, 2);
            break;
        case 25:
            dJydy[0] = 0.5*(-2*my25 + 2*y25)/pow(sigmay25, 2);
            break;
        case 26:
            dJydy[0] = 0.5*(-2*my26 + 2*y26)/pow(sigmay26, 2);
            break;
        case 27:
            dJydy[0] = 0.5*(-2*my27 + 2*y27)/pow(sigmay27, 2);
            break;
        case 28:
            dJydy[0] = 0.5*(-2*my28 + 2*y28)/pow(sigmay28, 2);
            break;
        case 29:
            dJydy[0] = 0.5*(-2*my29 + 2*y29)/pow(sigmay29, 2);
            break;
        case 30:
            dJydy[0] = 0.5*(-2*my30 + 2*y30)/pow(sigmay30, 2);
            break;
        case 31:
            dJydy[0] = 0.5*(-2*my31 + 2*y31)/pow(sigmay31, 2);
            break;
        case 32:
            dJydy[0] = 0.5*(-2*my32 + 2*y32)/pow(sigmay32, 2);
            break;
        case 33:
            dJydy[0] = 0.5*(-2*my33 + 2*y33)/pow(sigmay33, 2);
            break;
        case 34:
            dJydy[0] = 0.5*(-2*my34 + 2*y34)/pow(sigmay34, 2);
            break;
        case 35:
            dJydy[0] = 0.5*(-2*my35 + 2*y35)/pow(sigmay35, 2);
            break;
        case 36:
            dJydy[0] = 0.5*(-2*my36 + 2*y36)/pow(sigmay36, 2);
            break;
        case 37:
            dJydy[0] = 0.5*(-2*my37 + 2*y37)/pow(sigmay37, 2);
            break;
        case 38:
            dJydy[0] = 0.5*(-2*my38 + 2*y38)/pow(sigmay38, 2);
            break;
        case 39:
            dJydy[0] = 0.5*(-2*my39 + 2*y39)/pow(sigmay39, 2);
            break;
        case 40:
            dJydy[0] = 0.5*(-2*my40 + 2*y40)/pow(sigmay40, 2);
            break;
        case 41:
            dJydy[0] = 0.5*(-2*my41 + 2*y41)/pow(sigmay41, 2);
            break;
        case 42:
            dJydy[0] = 0.5*(-2*my42 + 2*y42)/pow(sigmay42, 2);
            break;
        case 43:
            dJydy[0] = 0.5*(-2*my43 + 2*y43)/pow(sigmay43, 2);
            break;
        case 44:
            dJydy[0] = 0.5*(-2*my44 + 2*y44)/pow(sigmay44, 2);
            break;
        case 45:
            dJydy[0] = 0.5*(-2*my45 + 2*y45)/pow(sigmay45, 2);
            break;
        case 46:
            dJydy[0] = 0.5*(-2*my46 + 2*y46)/pow(sigmay46, 2);
            break;
        case 47:
            dJydy[0] = 0.5*(-2*my47 + 2*y47)/pow(sigmay47, 2);
            break;
        case 48:
            dJydy[0] = 0.5*(-2*my48 + 2*y48)/pow(sigmay48, 2);
            break;
        case 49:
            dJydy[0] = 0.5*(-2*my49 + 2*y49)/pow(sigmay49, 2);
            break;
        case 50:
            dJydy[0] = 0.5*(-2*my50 + 2*y50)/pow(sigmay50, 2);
            break;
        case 51:
            dJydy[0] = 0.5*(-2*my51 + 2*y51)/pow(sigmay51, 2);
            break;
    }
}