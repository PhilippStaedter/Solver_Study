#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Ueda2001(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*CCc*Reaction7_T3/pow(CCc + Reaction7_k3, 2) + 1.0*Reaction7_T3/(CCc + Reaction7_k3);
    dwdx[1] = -1.0*Reaction11_parameter_0000073;
    dwdx[2] = 1.0*Reaction24_D0;
    dwdx[3] = -1.0*CCc*Reaction36_D9/pow(CCc + Reaction36_L9, 2) + 1.0*Reaction36_D9/(CCc + Reaction36_L9);
    dwdx[4] = -1.0*Reaction1_a*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(CCn*pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2)) + 1.0*Reaction1_a*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)/(CCn*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1));
    dwdx[5] = -1.0*Reaction3_a*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(CCn*pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2)) + 1.0*Reaction3_a*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)/(CCn*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1));
    dwdx[6] = -1.0*Reaction5_r*Reaction5_s5*pow(CCn/Reaction5_r3, Reaction5_r)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(CCn*pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2));
    dwdx[7] = -1.0*CCn*Reaction8_T4/pow(CCn + Reaction8_k4, 2) + 1.0*Reaction8_T4/(CCn + Reaction8_k4);
    dwdx[8] = 1.0*Reaction26_D0;
    dwdx[9] = -1.0*CCn*Reaction37_D10/pow(CCn + Reaction37_L10, 2) + 1.0*Reaction37_D10/(CCn + Reaction37_L10);
    dwdx[10] = 1.0*Reaction11_v3*species_0000012;
    dwdx[11] = 1.0*Reaction25_D0;
    dwdx[12] = -1.0*Clkc*Reaction35_D8/pow(Clkc + Reaction35_L8, 2) + 1.0*Reaction35_D8/(Clkc + Reaction35_L8);
    dwdx[13] = 1.0*Reaction6_D0;
    dwdx[14] = 1.0*Reaction18_s6;
    dwdx[15] = -1.0*Clkm*Reaction34_D7/pow(Clkm + Reaction34_L7, 2) + 1.0*Reaction34_D7/(Clkm + Reaction34_L7);
    dwdx[16] = 1.0*Reaction12_v1*Timc;
    dwdx[17] = 1.0*Reaction20_D0;
    dwdx[18] = -1.0*Perc*Reaction29_D2*species_0000013/pow(Perc + Reaction29_L2, 2) + 1.0*Reaction29_D2*species_0000013/(Perc + Reaction29_L2);
    dwdx[19] = 1.0*Reaction2_D0;
    dwdx[20] = 1.0*Reaction19_s2;
    dwdx[21] = -1.0*Perm*Reaction28_D1/pow(Perm + Reaction28_L1, 2) + 1.0*Reaction28_D1/(Perm + Reaction28_L1);
    dwdx[22] = -1.0*PTc*Reaction10_T1/pow(PTc + Reaction10_k1, 2) + 1.0*Reaction10_T1/(PTc + Reaction10_k1);
    dwdx[23] = -1.0*Reaction12_parameter_0000072;
    dwdx[24] = 1.0*Reaction21_D0;
    dwdx[25] = -1.0*PTc*Reaction32_D5/pow(PTc + Reaction32_L5, 2) + 1.0*Reaction32_D5/(PTc + Reaction32_L5);
    dwdx[26] = -1.0*Reaction1_r*Reaction1_s1*pow(PTn/Reaction1_r1, Reaction1_r)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(PTn*pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2));
    dwdx[27] = -1.0*Reaction3_r*Reaction3_s3*pow(PTn/Reaction3_r2, Reaction3_r)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(PTn*pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2));
    dwdx[28] = -1.0*Reaction5_a*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(PTn*pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2)) + 1.0*Reaction5_a*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)/(PTn*(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1));
    dwdx[29] = -1.0*PTn*Reaction9_T2/pow(PTn + Reaction9_k2, 2) + 1.0*Reaction9_T2/(PTn + Reaction9_k2);
    dwdx[30] = 1.0*Reaction23_D0;
    dwdx[31] = -1.0*PTn*Reaction33_D6/pow(PTn + Reaction33_L6, 2) + 1.0*Reaction33_D6/(PTn + Reaction33_L6);
    dwdx[32] = 1.0*Perc*Reaction12_v1;
    dwdx[33] = 1.0*Reaction27_D0;
    dwdx[34] = -1.0*Reaction31_D4*Timc/pow(Reaction31_L4 + Timc, 2) + 1.0*Reaction31_D4/(Reaction31_L4 + Timc);
    dwdx[35] = 1.0*Reaction4_D0;
    dwdx[36] = 1.0*Reaction16_s4;
    dwdx[37] = -1.0*Reaction30_D3*Timm/pow(Reaction30_L3 + Timm, 2) + 1.0*Reaction30_D3/(Reaction30_L3 + Timm);
    dwdx[38] = 1.0*Clkc*Reaction11_v3;
    dwdx[39] = 1.0*Perc*Reaction29_D2/(Perc + Reaction29_L2);
}