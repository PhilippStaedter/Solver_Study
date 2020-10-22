#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Ueda2001(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Reaction1_c1 + 1.0*Reaction1_s1*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1);
    w[1] = 1.0*Perm*Reaction2_D0;
    w[2] = 1.0*Reaction3_c2 + 1.0*Reaction3_s3*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1);
    w[3] = 1.0*Reaction4_D0*Timm;
    w[4] = 1.0*Reaction5_c3 + 1.0*Reaction5_s5*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1);
    w[5] = 1.0*Clkm*Reaction6_D0;
    w[6] = 1.0*CCc*Reaction7_T3/(CCc + Reaction7_k3);
    w[7] = 1.0*CCn*Reaction8_T4/(CCn + Reaction8_k4);
    w[8] = 1.0*PTn*Reaction9_T2/(PTn + Reaction9_k2);
    w[9] = 1.0*PTc*Reaction10_T1/(PTc + Reaction10_k1);
    w[10] = -1.0*CCc*Reaction11_parameter_0000073 + 1.0*Clkc*Reaction11_v3*species_0000012;
    w[11] = -1.0*PTc*Reaction12_parameter_0000072 + 1.0*Perc*Reaction12_v1*Timc;
    w[12] = 1.0*Reaction16_s4*Timm;
    w[13] = 1.0*Clkm*Reaction18_s6;
    w[14] = 1.0*Perm*Reaction19_s2;
    w[15] = 1.0*Perc*Reaction20_D0;
    w[16] = 1.0*PTc*Reaction21_D0;
    w[17] = 1.0*PTn*Reaction23_D0;
    w[18] = 1.0*CCc*Reaction24_D0;
    w[19] = 1.0*Clkc*Reaction25_D0;
    w[20] = 1.0*CCn*Reaction26_D0;
    w[21] = 1.0*Reaction27_D0*Timc;
    w[22] = 1.0*Perm*Reaction28_D1/(Perm + Reaction28_L1);
    w[23] = 1.0*Perc*Reaction29_D2*species_0000013/(Perc + Reaction29_L2);
    w[24] = 1.0*Reaction30_D3*Timm/(Reaction30_L3 + Timm);
    w[25] = 1.0*Reaction31_D4*Timc/(Reaction31_L4 + Timc);
    w[26] = 1.0*PTc*Reaction32_D5/(PTc + Reaction32_L5);
    w[27] = 1.0*PTn*Reaction33_D6/(PTn + Reaction33_L6);
    w[28] = 1.0*Clkm*Reaction34_D7/(Clkm + Reaction34_L7);
    w[29] = 1.0*Clkc*Reaction35_D8/(Clkc + Reaction35_L8);
    w[30] = 1.0*CCc*Reaction36_D9/(CCc + Reaction36_L9);
    w[31] = 1.0*CCn*Reaction37_D10/(CCn + Reaction37_L10);
}