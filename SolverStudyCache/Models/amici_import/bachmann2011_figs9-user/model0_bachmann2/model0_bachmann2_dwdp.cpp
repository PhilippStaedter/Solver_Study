#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_bachmann2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[8] = -0.27500000000000002*CISRNAEqc*CISRNATurn*npSTAT5/init_STAT5;
            dwdp[19] = -0.27500000000000002*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5/init_STAT5;
            break;
        case 1:
            dwdp[5] = 0.40000000000000002*CIS*CISInh*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(pow(CISEqc, 2)*pow(init_EpoRJAK2, 2)*pow(CIS*CISInh/CISEqc + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            dwdp[16] = 0.40000000000000002*CISRNA*CISTurn/CISRNAEqc;
            dwdp[18] = CISEqcOE*CISTurn*CISoe;
            break;
        case 2:
            dwdp[18] = CISEqc*CISTurn*CISoe;
            break;
        case 3:
            dwdp[5] = -0.40000000000000002*CIS*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(CISEqc*pow(init_EpoRJAK2, 2)*pow(CIS*CISInh/CISEqc + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            break;
        case 4:
            dwdp[9] = 0.27500000000000002*CISnRNA1;
            dwdp[10] = 0.27500000000000002*CISnRNA2;
            dwdp[12] = 0.27500000000000002*CISnRNA3;
            dwdp[13] = 0.27500000000000002*CISnRNA4;
            dwdp[14] = 0.27500000000000002*CISnRNA5;
            break;
        case 5:
            dwdp[8] = -0.27500000000000002*CISRNATurn*npSTAT5*(ActD - 1)/init_STAT5;
            dwdp[16] = -0.40000000000000002*CISEqc*CISRNA*CISTurn/pow(CISRNAEqc, 2);
            break;
        case 6:
            dwdp[8] = -0.27500000000000002*CISRNAEqc*npSTAT5*(ActD - 1)/init_STAT5;
            dwdp[15] = 0.40000000000000002*CISRNA;
            break;
        case 7:
            dwdp[16] = 0.40000000000000002*CISEqc*CISRNA/CISRNAEqc;
            dwdp[17] = 0.40000000000000002*CIS;
            dwdp[18] = CISEqc*CISEqcOE*CISoe;
            break;
        case 8:
            dwdp[18] = CISEqc*CISEqcOE*CISTurn;
            break;
        case 9:
            dwdp[22] = 0.40000000000000002*EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
            dwdp[30] = 1.2000000000000002*EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            dwdp[31] = 1.2000000000000002*p1EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            dwdp[32] = 0.40000000000000002*p2EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
            break;
        case 10:
            dwdp[30] = -1.2000000000000002*EpoRActJAK2*EpoRJAK2_CIS*EpoRpJAK2/(pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            dwdp[31] = -1.2000000000000002*EpoRActJAK2*EpoRJAK2_CIS*p1EpoRpJAK2/(pow(EpoRCISInh*EpoRJAK2_CIS + 1, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            break;
        case 11:
            dwdp[1] = 0.40000000000000002*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/init_EpoRJAK2;
            break;
        case 12:
            dwdp[0] = 0.40000000000000002*Epo*EpoRJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
            break;
        case 13:
            dwdp[11] = 0.40000000000000002*EpoRpJAK2*SHP1Act/init_SHP1;
            dwdp[33] = 0.40000000000000002*SHP1Act*p1EpoRpJAK2/init_SHP1;
            dwdp[34] = 0.40000000000000002*SHP1Act*p2EpoRpJAK2/init_SHP1;
            dwdp[35] = 0.40000000000000002*SHP1Act*p12EpoRpJAK2/init_SHP1;
            break;
        case 14:
            dwdp[2] = 0.40000000000000002*SHP1*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/init_EpoRJAK2;
            break;
        case 15:
            dwdp[3] = 0.40000000000000002*SHP1Act;
            break;
        case 16:
            dwdp[0] = 0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo*SOCS3*SOCS3Inh/(pow(SOCS3Eqc, 2)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[4] = 0.40000000000000002*SOCS3*SOCS3Inh*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(pow(SOCS3Eqc, 2)*init_EpoRJAK2*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[5] = 0.40000000000000002*SOCS3*SOCS3Inh*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(pow(SOCS3Eqc, 2)*pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[22] = 0.40000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3*SOCS3Inh/(pow(SOCS3Eqc, 2)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[27] = 0.40000000000000002*SOCS3RNA*SOCS3Turn/SOCS3RNAEqc;
            dwdp[29] = SOCS3EqcOE*SOCS3Turn*SOCS3oe;
            dwdp[30] = 1.2000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3*SOCS3Inh/(pow(SOCS3Eqc, 2)*(EpoRCISInh*EpoRJAK2_CIS + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[31] = 1.2000000000000002*EpoRActJAK2*SOCS3*SOCS3Inh*p1EpoRpJAK2/(pow(SOCS3Eqc, 2)*(EpoRCISInh*EpoRJAK2_CIS + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[32] = 0.40000000000000002*EpoRActJAK2*SOCS3*SOCS3Inh*p2EpoRpJAK2/(pow(SOCS3Eqc, 2)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            break;
        case 17:
            dwdp[29] = SOCS3Eqc*SOCS3Turn*SOCS3oe;
            break;
        case 18:
            dwdp[0] = -0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo*SOCS3/(SOCS3Eqc*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[4] = -0.40000000000000002*SOCS3*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(SOCS3Eqc*init_EpoRJAK2*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[5] = -0.40000000000000002*SOCS3*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(SOCS3Eqc*pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[22] = -0.40000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3/(SOCS3Eqc*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[30] = -1.2000000000000002*EpoRActJAK2*EpoRpJAK2*SOCS3/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[31] = -1.2000000000000002*EpoRActJAK2*SOCS3*p1EpoRpJAK2/(SOCS3Eqc*(EpoRCISInh*EpoRJAK2_CIS + 1)*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            dwdp[32] = -0.40000000000000002*EpoRActJAK2*SOCS3*p2EpoRpJAK2/(SOCS3Eqc*pow(SOCS3*SOCS3Inh/SOCS3Eqc + 1, 2));
            break;
        case 19:
            dwdp[20] = 0.27500000000000002*SOCS3nRNA1;
            dwdp[21] = 0.27500000000000002*SOCS3nRNA2;
            dwdp[23] = 0.27500000000000002*SOCS3nRNA3;
            dwdp[24] = 0.27500000000000002*SOCS3nRNA4;
            dwdp[25] = 0.27500000000000002*SOCS3nRNA5;
            break;
        case 20:
            dwdp[19] = -0.27500000000000002*SOCS3RNATurn*npSTAT5*(ActD - 1)/init_STAT5;
            dwdp[27] = -0.40000000000000002*SOCS3Eqc*SOCS3RNA*SOCS3Turn/pow(SOCS3RNAEqc, 2);
            break;
        case 21:
            dwdp[19] = -0.27500000000000002*SOCS3RNAEqc*npSTAT5*(ActD - 1)/init_STAT5;
            dwdp[26] = 0.40000000000000002*SOCS3RNA;
            break;
        case 22:
            dwdp[27] = 0.40000000000000002*SOCS3Eqc*SOCS3RNA/SOCS3RNAEqc;
            dwdp[28] = 0.40000000000000002*SOCS3;
            dwdp[29] = SOCS3Eqc*SOCS3EqcOE*SOCS3oe;
            break;
        case 23:
            dwdp[29] = SOCS3Eqc*SOCS3EqcOE*SOCS3Turn;
            break;
        case 24:
            dwdp[5] = 0.40000000000000002*STAT5*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            break;
        case 25:
            dwdp[4] = 0.40000000000000002*STAT5*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            break;
        case 26:
            dwdp[7] = 0.27500000000000002*npSTAT5;
            break;
        case 27:
            dwdp[6] = 0.40000000000000002*pSTAT5;
            break;
        case 29:
            dwdp[1] = -0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/pow(init_EpoRJAK2, 2);
            dwdp[2] = -0.40000000000000002*SHP1*SHP1ActEpoR*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/pow(init_EpoRJAK2, 2);
            dwdp[4] = -0.40000000000000002*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(pow(init_EpoRJAK2, 2)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            dwdp[5] = -0.80000000000000004*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(pow(init_EpoRJAK2, 3)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
            break;
        case 30:
            dwdp[11] = -0.40000000000000002*EpoRpJAK2*JAK2EpoRDeaSHP1*SHP1Act/pow(init_SHP1, 2);
            dwdp[33] = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2/pow(init_SHP1, 2);
            dwdp[34] = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2/pow(init_SHP1, 2);
            dwdp[35] = -0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2/pow(init_SHP1, 2);
            break;
        case 31:
            dwdp[8] = 0.27500000000000002*CISRNAEqc*CISRNATurn*npSTAT5*(ActD - 1)/pow(init_STAT5, 2);
            dwdp[19] = 0.27500000000000002*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5*(ActD - 1)/pow(init_STAT5, 2);
            break;
    }
}