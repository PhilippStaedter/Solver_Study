#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_moriya1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[6] = kscig;
            dwdp[20] = -kac10*(Cdc10T - s71)/pow(Cdc10T + Jac10 - s71, 2) + kac10/(Cdc10T + Jac10 - s71);
            dwdp[21] = Vamik;
            dwdp[35] = ksc18;
            break;
        case 1:
            dwdp[5] = k25_dash*s137;
            dwdp[8] = k255*k25_dash*s63;
            dwdp[13] = k255*k25_dash*s153;
            dwdp[33] = -Va25_dash2*s56*(Cdc25T - s83)/pow(Cdc25T + Ja25 - s83, 2) + Va25_dash2*s56/(Cdc25T + Ja25 - s83);
            dwdp[40] = k25_dash*s60;
            break;
        case 2:
            dwdp[33] = -Va25_dash2*s56*(Cdc25T - s83)/pow(Cdc25T + Ja25 - s83, 2);
            break;
        case 3:
            dwdp[20] = -kac10*(Cdc10T - s71)/pow(Cdc10T + Jac10 - s71, 2);
            break;
        case 4:
            dwdp[37] = -(-s50 + 1)*(kaie*s56 + kaie_dash*s75)/pow(Jaie - s50 + 1, 2);
            break;
        case 5:
            dwdp[49] = -kaslp*s50*(Slp1T - s48)/pow(Jaslp + Slp1T - s48, 2);
            break;
        case 6:
            dwdp[16] = -(Srw1T - s47)*(kasrw + kasrw_dash*s48)/pow(Jasrw + Srw1T - s47, 2);
            break;
        case 7:
            dwdp[26] = -(Vawee_dash + Vawee_dash2*s81)*(Wee1T - s80)/pow(Jawee + Wee1T - s80, 2);
            break;
        case 8:
            dwdp[34] = -s83*(Rad3*Vi25*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + Vi25_dash + Vi25_dash2*s81)/pow(Ji25 + s83, 2);
            break;
        case 9:
            dwdp[19] = -s71*(kic10 + kic10_dash*s67)/pow(Jic10 + s71, 2);
            break;
        case 10:
            dwdp[43] = -kiie*s50/pow(Jiie + s50, 2);
            break;
        case 11:
            dwdp[47] = kori*n*s91*pow((kipre*s56 + kipre_dash*s67)/Jipre, n)/(Jipre*pow(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1, 2));
            break;
        case 12:
            dwdp[50] = -kislp*s48/pow(Jislp + s48, 2);
            break;
        case 13:
            dwdp[27] = -s47*(Puc1*kisrw_dash3 + kisrw + kisrw_dash*s67 + kisrw_dash2*s56 + kisrw_dash4*s75)/pow(Jisrw + s47, 2);
            break;
        case 14:
            dwdp[28] = -s80*(Viwee_dash + Viwee_dash2*s56)/pow(Jiwee + s80, 2);
            break;
        case 15:
            dwdp[1] = kdrumpuc*s137;
            dwdp[3] = kdrumpuc*s161;
            dwdp[15] = kdrumpuc*s153;
            dwdp[17] = kdrumpuc*s149;
            dwdp[27] = kisrw_dash3*s47/(Jisrw + s47);
            dwdp[44] = kdrumpuc*s166;
            break;
        case 16:
            dwdp[8] = -beta*kpyp*s63*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -beta*kpyp*s153*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Vi25*s83*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(Ji25 + s83);
            break;
        case 17:
            dwdp[2] = kdcycslp_dash*s161;
            dwdp[4] = kdcycslp_dash*s137;
            dwdp[29] = kdcycslp_dash*s60;
            dwdp[30] = kdcycslp_dash*s56;
            dwdp[49] = -kaslp*s50*(Slp1T - s48)/pow(Jaslp + Slp1T - s48, 2) + kaslp*s50/(Jaslp + Slp1T - s48);
            break;
        case 18:
            dwdp[2] = kdcycsrw_dash*s161;
            dwdp[4] = kdcycsrw_dash*s137;
            dwdp[16] = -(Srw1T - s47)*(kasrw + kasrw_dash*s48)/pow(Jasrw + Srw1T - s47, 2) + (kasrw + kasrw_dash*s48)/(Jasrw + Srw1T - s47);
            dwdp[29] = kdcycsrw_dash*s60;
            dwdp[30] = kdcycsrw_dash*s56;
            break;
        case 19:
            dwdp[33] = s56*(Cdc25T - s83)/(Cdc25T + Ja25 - s83);
            break;
        case 20:
            dwdp[21] = Cdc10T;
            break;
        case 21:
            dwdp[21] = s71;
            break;
        case 22:
            dwdp[26] = (Wee1T - s80)/(Jawee + Wee1T - s80);
            break;
        case 23:
            dwdp[26] = s81*(Wee1T - s80)/(Jawee + Wee1T - s80);
            break;
        case 24:
            dwdp[34] = Rad3*s83*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(Ji25 + s83);
            break;
        case 25:
            dwdp[34] = s83/(Ji25 + s83);
            break;
        case 26:
            dwdp[34] = s81*s83/(Ji25 + s83);
            break;
        case 27:
            dwdp[22] = s72;
            break;
        case 28:
            dwdp[22] = s67*s72;
            break;
        case 29:
            dwdp[22] = s56*s72;
            break;
        case 30:
            dwdp[22] = s60*s72;
            break;
        case 31:
            dwdp[28] = s80/(Jiwee + s80);
            break;
        case 32:
            dwdp[28] = s56*s80/(Jiwee + s80);
            break;
        case 33:
            dwdp[26] = -(Vawee_dash + Vawee_dash2*s81)*(Wee1T - s80)/pow(Jawee + Wee1T - s80, 2) + (Vawee_dash + Vawee_dash2*s81)/(Jawee + Wee1T - s80);
            dwdp[31] = kwee_dash*s161;
            dwdp[41] = kwee_dash*s56;
            break;
        case 34:
            dwdp[8] = -Rad3*kpyp*s63*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*kpyp*s153*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            break;
        case 35:
            dwdp[8] = -Rad3*beta*kpyp*s63*s84/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*s84/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*s84/(Ji25 + s83);
            break;
        case 36:
            dwdp[8] = s63*(Cdc25T*k25_dash + s83*(-k25_dash + k25_dash2));
            dwdp[13] = s153*(Cdc25T*k25_dash + s83*(-k25_dash + k25_dash2));
            break;
        case 37:
            dwdp[5] = s137*(Cdc25T - s83);
            dwdp[8] = k255*s63*(Cdc25T - s83);
            dwdp[13] = k255*s153*(Cdc25T - s83);
            dwdp[40] = s60*(Cdc25T - s83);
            break;
        case 38:
            dwdp[5] = s137*s83;
            dwdp[8] = k255*s63*s83;
            dwdp[13] = k255*s153*s83;
            dwdp[40] = s60*s83;
            break;
        case 39:
            dwdp[8] = -Rad3*beta*kpyp*s63*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2)*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18));
            dwdp[13] = -Rad3*beta*kpyp*s153*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2)*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18));
            dwdp[34] = Rad3*Vi25*s83*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/((Ji25 + s83)*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18));
            break;
        case 40:
            dwdp[20] = (Cdc10T - s71)/(Cdc10T + Jac10 - s71);
            break;
        case 41:
            dwdp[37] = s56*(-s50 + 1)/(Jaie - s50 + 1);
            break;
        case 42:
            dwdp[37] = s75*(-s50 + 1)/(Jaie - s50 + 1);
            break;
        case 43:
            dwdp[49] = s50*(Slp1T - s48)/(Jaslp + Slp1T - s48);
            break;
        case 44:
            dwdp[16] = (Srw1T - s47)/(Jasrw + Srw1T - s47);
            break;
        case 45:
            dwdp[16] = s48*(Srw1T - s47)/(Jasrw + Srw1T - s47);
            break;
        case 46:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[36] = s84;
            dwdp[46] = -2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 47:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s56/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s56/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s56/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[36] = s56*s84;
            dwdp[46] = -2*oriT*s84*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s56/ko18 - s56*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s56/ko18)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 48:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s67/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s67/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s67/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[36] = s67*s84;
            dwdp[46] = -2*oriT*s84*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s67/ko18 - s67*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s67/ko18)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 49:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s63/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s63/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s63/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[36] = s63*s84;
            dwdp[46] = -2*oriT*s84*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-s63/ko18 - s63*(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - s63/ko18)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 50:
            dwdp[24] = s75;
            break;
        case 51:
            dwdp[24] = s48*s75;
            break;
        case 52:
            dwdp[24] = s47*s75;
            break;
        case 53:
            dwdp[10] = s149;
            dwdp[12] = s153;
            dwdp[14] = s67;
            dwdp[18] = s63;
            break;
        case 54:
            dwdp[10] = s149*s48;
            dwdp[12] = s153*s48;
            dwdp[14] = s48*s67;
            dwdp[18] = s48*s63;
            break;
        case 55:
            dwdp[2] = s161;
            dwdp[4] = s137;
            dwdp[29] = s60;
            dwdp[30] = s56;
            break;
        case 56:
            dwdp[2] = s161*s48;
            dwdp[4] = s137*s48;
            dwdp[29] = s48*s60;
            dwdp[30] = s48*s56;
            break;
        case 57:
            dwdp[2] = Slp1T*s161;
            dwdp[4] = Slp1T*s137;
            dwdp[29] = Slp1T*s60;
            dwdp[30] = Slp1T*s56;
            break;
        case 58:
            dwdp[2] = s161*s47;
            dwdp[4] = s137*s47;
            dwdp[29] = s47*s60;
            dwdp[30] = s47*s56;
            break;
        case 59:
            dwdp[2] = Srw1T*s161;
            dwdp[4] = Srw1T*s137;
            dwdp[29] = Srw1T*s60;
            dwdp[30] = Srw1T*s56;
            break;
        case 60:
            dwdp[38] = s81;
            break;
        case 61:
            dwdp[1] = s137;
            dwdp[3] = s161;
            dwdp[15] = s153;
            dwdp[17] = s149;
            dwdp[44] = s166;
            break;
        case 62:
            dwdp[1] = s137*s56;
            dwdp[3] = s161*s56;
            dwdp[15] = s153*s56;
            dwdp[17] = s149*s56;
            dwdp[44] = s166*s56;
            break;
        case 63:
            dwdp[1] = s137*s60;
            dwdp[3] = s161*s60;
            dwdp[15] = s153*s60;
            dwdp[17] = s149*s60;
            dwdp[44] = s166*s60;
            break;
        case 64:
            dwdp[1] = s137*s75;
            dwdp[3] = s161*s75;
            dwdp[15] = s153*s75;
            dwdp[17] = s149*s75;
            dwdp[44] = s166*s75;
            break;
        case 65:
            dwdp[1] = s137*s67;
            dwdp[3] = s161*s67;
            dwdp[15] = s153*s67;
            dwdp[17] = s149*s67;
            dwdp[44] = s166*s67;
            break;
        case 66:
            dwdp[1] = s137*s63;
            dwdp[3] = s161*s63;
            dwdp[15] = s153*s63;
            dwdp[17] = s149*s63;
            dwdp[44] = s166*s63;
            break;
        case 67:
            dwdp[1] = Puc1*s137;
            dwdp[3] = Puc1*s161;
            dwdp[15] = Puc1*s153;
            dwdp[17] = Puc1*s149;
            dwdp[44] = Puc1*s166;
            break;
        case 68:
            dwdp[19] = s71/(Jic10 + s71);
            break;
        case 69:
            dwdp[19] = s67*s71/(Jic10 + s71);
            break;
        case 70:
            dwdp[43] = s50/(Jiie + s50);
            break;
        case 71:
            dwdp[46] = s56*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(oriT - s90 - s91)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18);
            break;
        case 72:
            dwdp[46] = s67*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(oriT - s90 - s91)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18);
            break;
        case 73:
            dwdp[46] = s63*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(oriT - s90 - s91)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18);
            break;
        case 74:
            dwdp[47] = -kori*n*s56*s91*pow((kipre*s56 + kipre_dash*s67)/Jipre, n)/(pow(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1, 2)*(kipre*s56 + kipre_dash*s67));
            break;
        case 75:
            dwdp[47] = -kori*n*s67*s91*pow((kipre*s56 + kipre_dash*s67)/Jipre, n)/(pow(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1, 2)*(kipre*s56 + kipre_dash*s67));
            break;
        case 76:
            dwdp[50] = s48/(Jislp + s48);
            break;
        case 77:
            dwdp[27] = s47/(Jisrw + s47);
            break;
        case 78:
            dwdp[27] = s47*s67/(Jisrw + s47);
            break;
        case 79:
            dwdp[27] = s47*s56/(Jisrw + s47);
            break;
        case 80:
            dwdp[27] = Puc1*s47/(Jisrw + s47);
            break;
        case 81:
            dwdp[27] = s47*s75/(Jisrw + s47);
            break;
        case 83:
            dwdp[31] = s161*s72;
            dwdp[41] = s56*s72;
            break;
        case 84:
            dwdp[7] = s67*s72;
            dwdp[32] = s149*s72;
            break;
        case 85:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[46] = -2*oriT*s84*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*((kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2) + (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/(pow(ko18, 2)*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/pow(ko18, 2))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 86:
            dwdp[8] = -Rad3*beta*kpyp*s63*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(-2*k_dash*oriT*s84*s90*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))/(Ji25 + s83);
            dwdp[46] = -2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/((-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)*pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1/ko18 - (oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/(ko18*sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2))))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 1/ko18)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2);
            break;
        case 87:
            dwdp[47] = s91/(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1);
            break;
        case 88:
            dwdp[8] = s63/(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1);
            dwdp[13] = s153/(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1);
            break;
        case 89:
            dwdp[5] = s137;
            dwdp[40] = s60;
            break;
        case 90:
            dwdp[48] = s90;
            break;
        case 91:
            dwdp[35] = Cdc10T;
            break;
        case 92:
            dwdp[35] = s71;
            break;
        case 93:
            dwdp[23] = 1;
            break;
        case 94:
            dwdp[6] = Cdc10T;
            break;
        case 95:
            dwdp[6] = s71;
            break;
        case 96:
            dwdp[0] = 1;
            break;
        case 97:
            dwdp[25] = 1;
            break;
        case 98:
            dwdp[25] = s48;
            break;
        case 99:
            dwdp[45] = 1;
            break;
        case 100:
            dwdp[31] = s161*(Wee1T - s80);
            dwdp[41] = s56*(Wee1T - s80);
            break;
        case 101:
            dwdp[31] = s161*s80;
            dwdp[41] = s56*s80;
            break;
        case 102:
            dwdp[9] = -s149;
            dwdp[11] = -s153;
            break;
        case 103:
            dwdp[9] = s166*s67;
            dwdp[11] = s166*s63;
            break;
        case 104:
            dwdp[39] = -s161;
            dwdp[42] = -s137;
            break;
        case 105:
            dwdp[39] = s166*s56;
            dwdp[42] = s166*s60;
            break;
        case 106:
            dwdp[47] = -kori*s91*pow((kipre*s56 + kipre_dash*s67)/Jipre, n)*log((kipre*s56 + kipre_dash*s67)/Jipre)/pow(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1, 2);
            break;
        case 107:
            dwdp[8] = -Rad3*beta*kpyp*s63*(k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + k_dash*s90*(-2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[13] = -Rad3*beta*kpyp*s153*(k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + k_dash*s90*(-2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1, 2);
            dwdp[34] = Rad3*Vi25*s83*(k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + k_dash*s90*(-2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))/(Ji25 + s83);
            dwdp[46] = (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/pow(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) + (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + (-2*oriT*s84*(-1 - (oriT - s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)/sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)))/pow(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2) - 2*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18))*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18);
            break;
    }
}