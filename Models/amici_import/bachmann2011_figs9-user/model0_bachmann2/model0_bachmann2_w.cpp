#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bachmann2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 0.40000000000000002*Epo*EpoRJAK2*JAK2ActEpo/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
    w[1] = 0.40000000000000002*EpoRCISRemove*EpoRJAK2_CIS*(p12EpoRpJAK2 + p1EpoRpJAK2)/init_EpoRJAK2;
    w[2] = 0.40000000000000002*SHP1*SHP1ActEpoR*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/init_EpoRJAK2;
    w[3] = 0.40000000000000002*SHP1Act*SHP1Dea;
    w[4] = 0.40000000000000002*STAT5*STAT5ActJAK2*(EpoRpJAK2 + p12EpoRpJAK2 + p1EpoRpJAK2 + p2EpoRpJAK2)/(init_EpoRJAK2*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
    w[5] = 0.40000000000000002*STAT5*STAT5ActEpoR*pow(p12EpoRpJAK2 + p1EpoRpJAK2, 2)/(pow(init_EpoRJAK2, 2)*(CIS*CISInh/CISEqc + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
    w[6] = 0.40000000000000002*STAT5Imp*pSTAT5;
    w[7] = 0.27500000000000002*STAT5Exp*npSTAT5;
    w[8] = -0.27500000000000002*CISRNAEqc*CISRNATurn*npSTAT5*(ActD - 1)/init_STAT5;
    w[9] = 0.27500000000000002*CISRNADelay*CISnRNA1;
    w[10] = 0.27500000000000002*CISRNADelay*CISnRNA2;
    w[11] = 0.40000000000000002*EpoRpJAK2*JAK2EpoRDeaSHP1*SHP1Act/init_SHP1;
    w[12] = 0.27500000000000002*CISRNADelay*CISnRNA3;
    w[13] = 0.27500000000000002*CISRNADelay*CISnRNA4;
    w[14] = 0.27500000000000002*CISRNADelay*CISnRNA5;
    w[15] = 0.40000000000000002*CISRNA*CISRNATurn;
    w[16] = 0.40000000000000002*CISEqc*CISRNA*CISTurn/CISRNAEqc;
    w[17] = 0.40000000000000002*CIS*CISTurn;
    w[18] = CISEqc*CISEqcOE*CISTurn*CISoe;
    w[19] = -0.27500000000000002*SOCS3RNAEqc*SOCS3RNATurn*npSTAT5*(ActD - 1)/init_STAT5;
    w[20] = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA1;
    w[21] = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA2;
    w[22] = 0.40000000000000002*EpoRActJAK2*EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
    w[23] = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA3;
    w[24] = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA4;
    w[25] = 0.27500000000000002*SOCS3RNADelay*SOCS3nRNA5;
    w[26] = 0.40000000000000002*SOCS3RNA*SOCS3RNATurn;
    w[27] = 0.40000000000000002*SOCS3Eqc*SOCS3RNA*SOCS3Turn/SOCS3RNAEqc;
    w[28] = 0.40000000000000002*SOCS3*SOCS3Turn;
    w[29] = SOCS3Eqc*SOCS3EqcOE*SOCS3Turn*SOCS3oe;
    w[30] = 1.2000000000000002*EpoRActJAK2*EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
    w[31] = 1.2000000000000002*EpoRActJAK2*p1EpoRpJAK2/((EpoRCISInh*EpoRJAK2_CIS + 1)*(SOCS3*SOCS3Inh/SOCS3Eqc + 1));
    w[32] = 0.40000000000000002*EpoRActJAK2*p2EpoRpJAK2/(SOCS3*SOCS3Inh/SOCS3Eqc + 1);
    w[33] = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p1EpoRpJAK2/init_SHP1;
    w[34] = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p2EpoRpJAK2/init_SHP1;
    w[35] = 0.40000000000000002*JAK2EpoRDeaSHP1*SHP1Act*p12EpoRpJAK2/init_SHP1;
}