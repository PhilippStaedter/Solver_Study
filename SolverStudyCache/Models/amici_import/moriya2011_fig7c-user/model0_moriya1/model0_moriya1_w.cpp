#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_moriya1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = kscyc;
    w[1] = s137*(Puc1*kdrumpuc + kdrum + kdrumc13*s56 + kdrumc13_dash*s60 + kdrumci1*s75 + kdrumcig*s67 + kdrumcig_dash*s63);
    w[2] = s161*(Slp1T*kdcycslp_dash + Srw1T*kdcycsrw_dash + kdcyc + kdcycslp*s48 + kdcycsrw*s47);
    w[3] = s161*(Puc1*kdrumpuc + kdrum + kdrumc13*s56 + kdrumc13_dash*s60 + kdrumci1*s75 + kdrumcig*s67 + kdrumcig_dash*s63);
    w[4] = s137*(Slp1T*kdcycslp_dash + Srw1T*kdcycsrw_dash + kdcyc + kdcycslp*s48 + kdcycsrw*s47);
    w[5] = s137*(Cdc25T*k25_dash + kpyp2 + s83*(-k25_dash + k25_dash2));
    w[6] = Cdc10T*kscig + kscig_dash*s71;
    w[7] = kmik_dash2*s67*s72;
    w[8] = k255*s63*(Cdc25T*k25_dash + s83*(-k25_dash + k25_dash2)) + kpyp*s63/(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1);
    w[9] = -lcm*s149 + lcp*s166*s67;
    w[10] = s149*(kdcig + kdcig_dash*s48);
    w[11] = -lcm*s153 + lcp*s166*s63;
    w[12] = s153*(kdcig + kdcig_dash*s48);
    w[13] = k255*s153*(Cdc25T*k25_dash + s83*(-k25_dash + k25_dash2)) + kpyp*s153/(Rad3*beta*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + 1);
    w[14] = s67*(kdcig + kdcig_dash*s48);
    w[15] = s153*(Puc1*kdrumpuc + kdrum + kdrumc13*s56 + kdrumc13_dash*s60 + kdrumci1*s75 + kdrumcig*s67 + kdrumcig_dash*s63);
    w[16] = (Srw1T - s47)*(kasrw + kasrw_dash*s48)/(Jasrw + Srw1T - s47);
    w[17] = s149*(Puc1*kdrumpuc + kdrum + kdrumc13*s56 + kdrumc13_dash*s60 + kdrumci1*s75 + kdrumcig*s67 + kdrumcig_dash*s63);
    w[18] = s63*(kdcig + kdcig_dash*s48);
    w[19] = s71*(kic10 + kic10_dash*s67)/(Jic10 + s71);
    w[20] = kac10*(Cdc10T - s71)/(Cdc10T + Jac10 - s71);
    w[21] = Cdc10T*Vamik + Vamik_dash*s71;
    w[22] = s72*(Vimik + Vimik_dash*s67 + Vimik_dash2*s56 + Vimik_dash3*s60);
    w[23] = ksci1;
    w[24] = s75*(kdci1 + kdci1_dash*s48 + kdci1_dash2*s47);
    w[25] = ksflp + ksflp_dash*s48;
    w[26] = (Vawee_dash + Vawee_dash2*s81)*(Wee1T - s80)/(Jawee + Wee1T - s80);
    w[27] = s47*(Puc1*kisrw_dash3 + kisrw + kisrw_dash*s67 + kisrw_dash2*s56 + kisrw_dash4*s75)/(Jisrw + s47);
    w[28] = s80*(Viwee_dash + Viwee_dash2*s56)/(Jiwee + s80);
    w[29] = s60*(Slp1T*kdcycslp_dash + Srw1T*kdcycsrw_dash + kdcyc + kdcycslp*s48 + kdcycsrw*s47);
    w[30] = s56*(Slp1T*kdcycslp_dash + Srw1T*kdcycsrw_dash + kdcyc + kdcycslp*s48 + kdcycsrw*s47);
    w[31] = s161*(Wee1T*kwee_dash + kmik_dash*s72 + s80*(-kwee_dash + kwee_dash2));
    w[32] = kmik_dash2*s149*s72;
    w[33] = Va25_dash2*s56*(Cdc25T - s83)/(Cdc25T + Ja25 - s83);
    w[34] = s83*(Rad3*Vi25*(amici_k*s84 + k_dash*s90*(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18)) + Vi25_dash + Vi25_dash2*s81)/(Ji25 + s83);
    w[35] = Cdc10T*ksc18 + ksc18_dash*s71;
    w[36] = s84*(kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63);
    w[37] = (-s50 + 1)*(kaie*s56 + kaie_dash*s75)/(Jaie - s50 + 1);
    w[38] = kdflp*s81;
    w[39] = -lm*s161 + lp*s166*s56;
    w[40] = s60*(Cdc25T*k25_dash + kpyp2 + s83*(-k25_dash + k25_dash2));
    w[41] = s56*(Wee1T*kwee_dash + kmik_dash*s72 + s80*(-kwee_dash + kwee_dash2));
    w[42] = -lm*s137 + lp*s166*s60;
    w[43] = kiie*s50/(Jiie + s50);
    w[44] = s166*(Puc1*kdrumpuc + kdrum + kdrumc13*s56 + kdrumc13_dash*s60 + kdrumci1*s75 + kdrumcig*s67 + kdrumcig_dash*s63);
    w[45] = ksrum;
    w[46] = (-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84)*(oriT - s90 - s91)*(kini_dash*s56 + kini_dash2*s67 + kini_dash3*s63)/(-2*oriT*s84/(oriT + s84 + sqrt(-4*oriT*s84 + pow(oriT + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18, 2)) + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18) + s84 + (kdc18 + kdc18c13*s56 + kdc18cig*s67 + kdc18cig_dash*s63 + ko18r)/ko18);
    w[47] = kori*s91/(pow((kipre*s56 + kipre_dash*s67)/Jipre, n) + 1);
    w[48] = krepl*s90;
    w[49] = kaslp*s50*(Slp1T - s48)/(Jaslp + Slp1T - s48);
    w[50] = kislp*s48/(Jislp + s48);
}