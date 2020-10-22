#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_schilling1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*SHP1_delay;
    dwdx[1] = 1.0*SHP1_delay;
    dwdx[2] = 1.0*SHP1_delay;
    dwdx[3] = 1.0*SHP1_delay;
    dwdx[4] = 1.0*SHP1_delay;
    dwdx[5] = 1.0*SHP1_delay;
    dwdx[6] = 1.0*SHP1_delay;
    dwdx[7] = 1.0*SHP1_delay;
    dwdx[8] = 1.0*First_ERK1_phosphorylation_by_ppMEK*ppMEK2;
    dwdx[9] = 1.0*First_ERK1_phosphorylation_by_ppMEK*ppMEK1;
    dwdx[10] = 1.0*First_ERK2_phosphorylation_by_ppMEK*ppMEK2;
    dwdx[11] = 1.0*First_ERK2_phosphorylation_by_ppMEK*ppMEK1;
    dwdx[12] = 1.0*JAK2*JAK2_phosphorylation_by_Epo;
    dwdx[13] = 1.0*EpoR_phosphorylation_by_pJAK2*pJAK2;
    dwdx[14] = 1.0*Epo*JAK2_phosphorylation_by_Epo;
    dwdx[15] = 1.0*First_MEK1_phosphorylation_by_pRaf*pRaf;
    dwdx[16] = 1.0*First_MEK2_phosphorylation_by_pRaf*pRaf;
    dwdx[17] = 1.0*mSOS*mSOS_induced_Raf_phosphorylation;
    dwdx[18] = 1.0*SHP1_activation_by_pEpoR*pEpoR;
    dwdx[19] = 1.0*SOS_recruitment_by_pEpoR*pEpoR;
    dwdx[20] = 1.0*actSHP1_deactivation;
    dwdx[21] = 1.0*pEpoR*pEpoR_dephosphorylation_by_actSHP1;
    dwdx[22] = 1.0*pJAK2*pJAK2_dephosphorylation_by_actSHP1;
    dwdx[23] = 1.0*SHP1_delay;
    dwdx[24] = 1.0*mSOS_release_from_membrane;
    dwdx[25] = 1.0*Raf*mSOS_induced_Raf_phosphorylation;
    dwdx[26] = 1.0*ppERK1*ppERK_neg_feedback_on_mSOS;
    dwdx[27] = 1.0*ppERK2*ppERK_neg_feedback_on_mSOS;
    dwdx[28] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*ppMEK2;
    dwdx[29] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*ppMEK1;
    dwdx[30] = 1.0*Second_ERK_dephosphorylation;
    dwdx[31] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*ppMEK2;
    dwdx[32] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*ppMEK1;
    dwdx[33] = 1.0*Second_ERK_dephosphorylation;
    dwdx[34] = 1.0*actSHP1*pEpoR_dephosphorylation_by_actSHP1;
    dwdx[35] = 1.0*SOS*SOS_recruitment_by_pEpoR;
    dwdx[36] = 1.0*SHP1*SHP1_activation_by_pEpoR;
    dwdx[37] = 1.0*actSHP1*pJAK2_dephosphorylation_by_actSHP1;
    dwdx[38] = 1.0*EpoR*EpoR_phosphorylation_by_pJAK2;
    dwdx[39] = 1.0*Second_MEK1_phosphorylation_by_pRaf*pRaf;
    dwdx[40] = 1.0*Second_MEK_dephosphorylation;
    dwdx[41] = 1.0*Second_MEK2_phosphorylation_by_pRaf*pRaf;
    dwdx[42] = 1.0*Second_MEK_dephosphorylation;
    dwdx[43] = 1.0*pRaf_dephosphorylation;
    dwdx[44] = 1.0*First_MEK2_phosphorylation_by_pRaf*MEK2;
    dwdx[45] = 1.0*First_MEK1_phosphorylation_by_pRaf*MEK1;
    dwdx[46] = 1.0*Second_MEK2_phosphorylation_by_pRaf*pMEK2;
    dwdx[47] = 1.0*Second_MEK1_phosphorylation_by_pRaf*pMEK1;
    dwdx[48] = 1.0*pSOS_dephosphorylation;
    dwdx[49] = 1.0*First_ERK_dephosphorylation;
    dwdx[50] = 1.0*mSOS*ppERK_neg_feedback_on_mSOS;
    dwdx[51] = 1.0*First_ERK_dephosphorylation;
    dwdx[52] = 1.0*mSOS*ppERK_neg_feedback_on_mSOS;
    dwdx[53] = 1.0*First_MEK_dephosphorylation;
    dwdx[54] = 1.0*ERK1*First_ERK1_phosphorylation_by_ppMEK;
    dwdx[55] = 1.0*ERK2*First_ERK2_phosphorylation_by_ppMEK;
    dwdx[56] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*pERK1;
    dwdx[57] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*pERK2;
    dwdx[58] = 1.0*First_MEK_dephosphorylation;
    dwdx[59] = 1.0*ERK1*First_ERK1_phosphorylation_by_ppMEK;
    dwdx[60] = 1.0*ERK2*First_ERK2_phosphorylation_by_ppMEK;
    dwdx[61] = 1.0*Second_ERK1_phosphorylation_by_ppMEK*pERK1;
    dwdx[62] = 1.0*Second_ERK2_phosphorylation_by_ppMEK*pERK2;
}