#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Ueda2001(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*Reaction1_s1*pow(PTn/Reaction1_r1, Reaction1_r)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))*log(PTn/Reaction1_r1)/pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2);
            break;
        case 1:
            dwdp[0] = 1.0*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1);
            break;
        case 2:
            dwdp[0] = 1.0*Reaction1_r*Reaction1_s1*pow(PTn/Reaction1_r1, Reaction1_r)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(Reaction1_r1*pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2));
            break;
        case 3:
            dwdp[0] = 1.0;
            break;
        case 4:
            dwdp[0] = -1.0*Reaction1_s1*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2) + 1.0*Reaction1_s1/(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1);
            break;
        case 5:
            dwdp[0] = 1.0*Reaction1_a*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))/(Reaction1_A1*pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2)) - 1.0*Reaction1_a*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)/(Reaction1_A1*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1));
            break;
        case 6:
            dwdp[0] = -1.0*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)*(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a))*log(CCn/Reaction1_A1)/pow(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1, 2) + 1.0*Reaction1_s1*pow(CCn/Reaction1_A1, Reaction1_a)*log(CCn/Reaction1_A1)/(Reaction1_B1 + pow(CCn/Reaction1_A1, Reaction1_a) + pow(PTn/Reaction1_r1, Reaction1_r) + 1);
            break;
        case 7:
            dwdp[1] = 1.0*Perm;
            break;
        case 8:
            dwdp[2] = -1.0*Reaction3_s3*pow(PTn/Reaction3_r2, Reaction3_r)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))*log(PTn/Reaction3_r2)/pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2);
            break;
        case 9:
            dwdp[2] = 1.0*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1);
            break;
        case 10:
            dwdp[2] = 1.0*Reaction3_r*Reaction3_s3*pow(PTn/Reaction3_r2, Reaction3_r)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(Reaction3_r2*pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2));
            break;
        case 11:
            dwdp[2] = 1.0;
            break;
        case 12:
            dwdp[2] = -1.0*Reaction3_s3*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2) + 1.0*Reaction3_s3/(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1);
            break;
        case 13:
            dwdp[2] = 1.0*Reaction3_a*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))/(Reaction3_A2*pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2)) - 1.0*Reaction3_a*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)/(Reaction3_A2*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1));
            break;
        case 14:
            dwdp[2] = -1.0*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)*(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a))*log(CCn/Reaction3_A2)/pow(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1, 2) + 1.0*Reaction3_s3*pow(CCn/Reaction3_A2, Reaction3_a)*log(CCn/Reaction3_A2)/(Reaction3_B2 + pow(CCn/Reaction3_A2, Reaction3_a) + pow(PTn/Reaction3_r2, Reaction3_r) + 1);
            break;
        case 15:
            dwdp[3] = 1.0*Timm;
            break;
        case 16:
            dwdp[4] = -1.0*Reaction5_s5*pow(CCn/Reaction5_r3, Reaction5_r)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))*log(CCn/Reaction5_r3)/pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2);
            break;
        case 17:
            dwdp[4] = 1.0*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1);
            break;
        case 18:
            dwdp[4] = 1.0*Reaction5_r*Reaction5_s5*pow(CCn/Reaction5_r3, Reaction5_r)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(Reaction5_r3*pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2));
            break;
        case 19:
            dwdp[4] = 1.0;
            break;
        case 20:
            dwdp[4] = -1.0*Reaction5_s5*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2) + 1.0*Reaction5_s5/(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1);
            break;
        case 21:
            dwdp[4] = 1.0*Reaction5_a*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))/(Reaction5_A3*pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2)) - 1.0*Reaction5_a*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)/(Reaction5_A3*(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1));
            break;
        case 22:
            dwdp[4] = -1.0*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)*(Reaction5_B3 + pow(PTn/Reaction5_A3, Reaction5_a))*log(PTn/Reaction5_A3)/pow(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1, 2) + 1.0*Reaction5_s5*pow(PTn/Reaction5_A3, Reaction5_a)*log(PTn/Reaction5_A3)/(Reaction5_B3 + pow(CCn/Reaction5_r3, Reaction5_r) + pow(PTn/Reaction5_A3, Reaction5_a) + 1);
            break;
        case 23:
            dwdp[5] = 1.0*Clkm;
            break;
        case 24:
            dwdp[6] = 1.0*CCc/(CCc + Reaction7_k3);
            break;
        case 25:
            dwdp[6] = -1.0*CCc*Reaction7_T3/pow(CCc + Reaction7_k3, 2);
            break;
        case 26:
            dwdp[7] = 1.0*CCn/(CCn + Reaction8_k4);
            break;
        case 27:
            dwdp[7] = -1.0*CCn*Reaction8_T4/pow(CCn + Reaction8_k4, 2);
            break;
        case 28:
            dwdp[8] = 1.0*PTn/(PTn + Reaction9_k2);
            break;
        case 29:
            dwdp[8] = -1.0*PTn*Reaction9_T2/pow(PTn + Reaction9_k2, 2);
            break;
        case 30:
            dwdp[9] = 1.0*PTc/(PTc + Reaction10_k1);
            break;
        case 31:
            dwdp[9] = -1.0*PTc*Reaction10_T1/pow(PTc + Reaction10_k1, 2);
            break;
        case 32:
            dwdp[10] = -1.0*CCc;
            break;
        case 33:
            dwdp[10] = 1.0*Clkc*species_0000012;
            break;
        case 34:
            dwdp[11] = -1.0*PTc;
            break;
        case 35:
            dwdp[11] = 1.0*Perc*Timc;
            break;
        case 36:
            dwdp[12] = 1.0*Timm;
            break;
        case 37:
            dwdp[13] = 1.0*Clkm;
            break;
        case 38:
            dwdp[14] = 1.0*Perm;
            break;
        case 39:
            dwdp[15] = 1.0*Perc;
            break;
        case 40:
            dwdp[16] = 1.0*PTc;
            break;
        case 41:
            dwdp[17] = 1.0*PTn;
            break;
        case 42:
            dwdp[18] = 1.0*CCc;
            break;
        case 43:
            dwdp[19] = 1.0*Clkc;
            break;
        case 44:
            dwdp[20] = 1.0*CCn;
            break;
        case 45:
            dwdp[21] = 1.0*Timc;
            break;
        case 46:
            dwdp[22] = -1.0*Perm*Reaction28_D1/pow(Perm + Reaction28_L1, 2);
            break;
        case 47:
            dwdp[22] = 1.0*Perm/(Perm + Reaction28_L1);
            break;
        case 48:
            dwdp[23] = -1.0*Perc*Reaction29_D2*species_0000013/pow(Perc + Reaction29_L2, 2);
            break;
        case 49:
            dwdp[23] = 1.0*Perc*species_0000013/(Perc + Reaction29_L2);
            break;
        case 50:
            dwdp[24] = -1.0*Reaction30_D3*Timm/pow(Reaction30_L3 + Timm, 2);
            break;
        case 51:
            dwdp[24] = 1.0*Timm/(Reaction30_L3 + Timm);
            break;
        case 52:
            dwdp[25] = -1.0*Reaction31_D4*Timc/pow(Reaction31_L4 + Timc, 2);
            break;
        case 53:
            dwdp[25] = 1.0*Timc/(Reaction31_L4 + Timc);
            break;
        case 54:
            dwdp[26] = -1.0*PTc*Reaction32_D5/pow(PTc + Reaction32_L5, 2);
            break;
        case 55:
            dwdp[26] = 1.0*PTc/(PTc + Reaction32_L5);
            break;
        case 56:
            dwdp[27] = -1.0*PTn*Reaction33_D6/pow(PTn + Reaction33_L6, 2);
            break;
        case 57:
            dwdp[27] = 1.0*PTn/(PTn + Reaction33_L6);
            break;
        case 58:
            dwdp[28] = -1.0*Clkm*Reaction34_D7/pow(Clkm + Reaction34_L7, 2);
            break;
        case 59:
            dwdp[28] = 1.0*Clkm/(Clkm + Reaction34_L7);
            break;
        case 60:
            dwdp[29] = -1.0*Clkc*Reaction35_D8/pow(Clkc + Reaction35_L8, 2);
            break;
        case 61:
            dwdp[29] = 1.0*Clkc/(Clkc + Reaction35_L8);
            break;
        case 62:
            dwdp[30] = -1.0*CCc*Reaction36_D9/pow(CCc + Reaction36_L9, 2);
            break;
        case 63:
            dwdp[30] = 1.0*CCc/(CCc + Reaction36_L9);
            break;
        case 64:
            dwdp[31] = -1.0*CCn*Reaction37_D10/pow(CCn + Reaction37_L10, 2);
            break;
        case 65:
            dwdp[31] = 1.0*CCn/(CCn + Reaction37_L10);
            break;
    }
}