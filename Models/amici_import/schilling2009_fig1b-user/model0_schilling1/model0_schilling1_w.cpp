#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_schilling1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Epo*JAK2*JAK2_phosphorylation_by_Epo;
    w[1] = 1.0*Delay06_mSHP1*SHP1_delay;
    w[2] = 1.0*Delay07_mSHP1*SHP1_delay;
    w[3] = 1.0*Delay08_mSHP1*SHP1_delay;
    w[4] = 1.0*actSHP1*actSHP1_deactivation;
    w[5] = 1.0*actSHP1*pEpoR*pEpoR_dephosphorylation_by_actSHP1;
    w[6] = 1.0*actSHP1*pJAK2*pJAK2_dephosphorylation_by_actSHP1;
    w[7] = 1.0*SOS*SOS_recruitment_by_pEpoR*pEpoR;
    w[8] = 1.0*mSOS*mSOS_release_from_membrane;
    w[9] = 1.0*Raf*mSOS*mSOS_induced_Raf_phosphorylation;
    w[10] = 1.0*pRaf*pRaf_dephosphorylation;
    w[11] = 1.0*EpoR*EpoR_phosphorylation_by_pJAK2*pJAK2;
    w[12] = 1.0*First_MEK2_phosphorylation_by_pRaf*MEK2*pRaf;
    w[13] = 1.0*First_MEK1_phosphorylation_by_pRaf*MEK1*pRaf;
    w[14] = 1.0*Second_MEK2_phosphorylation_by_pRaf*pMEK2*pRaf;
    w[15] = 1.0*Second_MEK1_phosphorylation_by_pRaf*pMEK1*pRaf;
    w[16] = 1.0*First_MEK_dephosphorylation*ppMEK2;
    w[17] = 1.0*First_MEK_dephosphorylation*ppMEK1;
    w[18] = 1.0*Second_MEK_dephosphorylation*pMEK2;
    w[19] = 1.0*Second_MEK_dephosphorylation*pMEK1;
    w[20] = 1.0*ERK1*First_ERK1_phosphorylation_by_ppMEK*ppMEK2;
    w[21] = 1.0*ERK2*First_ERK2_phosphorylation_by_ppMEK*ppMEK2;
    w[22] = 1.0*SHP1*SHP1_activation_by_pEpoR*pEpoR;
    w[23] = 1.0*ERK1*First_ERK1_phosphorylation_by_ppMEK*ppMEK1;
    w[24] = 1.0*ERK2*First_ERK2_phosphorylation_by_ppMEK*ppMEK1;
    w[25] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*pERK1*ppMEK2;
    w[26] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*pERK2*ppMEK2;
    w[27] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*pERK1*ppMEK1;
    w[28] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*pERK2*ppMEK1;
    w[29] = 1.0*First_ERK_dephosphorylation*ppERK1;
    w[30] = 1.0*First_ERK_dephosphorylation*ppERK2;
    w[31] = 1.0*Second_ERK_dephosphorylation*pERK1;
    w[32] = 1.0*Second_ERK_dephosphorylation*pERK2;
    w[33] = 1.0*SHP1_delay*mSHP1;
    w[34] = 1.0*mSOS*ppERK1*ppERK_neg_feedback_on_mSOS;
    w[35] = 1.0*mSOS*ppERK2*ppERK_neg_feedback_on_mSOS;
    w[36] = 1.0*pSOS*pSOS_dephosphorylation;
    w[37] = 1.0*Delay01_mSHP1*SHP1_delay;
    w[38] = 1.0*Delay02_mSHP1*SHP1_delay;
    w[39] = 1.0*Delay03_mSHP1*SHP1_delay;
    w[40] = 1.0*Delay04_mSHP1*SHP1_delay;
    w[41] = 1.0*Delay05_mSHP1*SHP1_delay;
}