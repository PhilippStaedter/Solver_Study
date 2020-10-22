#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_chassagnole2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*vALDO_rmaxALDO*(-cdhap*cgap/vALDO_kALDOeq + cfdp)/(cdhap*cgap/(vALDO_VALDOblf*vALDO_kALDOeq) + cdhap*vALDO_kALDOgap/(vALDO_VALDOblf*vALDO_kALDOeq) + cfdp*cgap/vALDO_kALDOgapinh + cfdp + cgap*vALDO_kALDOdhap/(vALDO_VALDOblf*vALDO_kALDOeq) + vALDO_kALDOfdp);
    w[1] = 1.0*pow(ce4p, vDAHPS_nDAHPSe4p)*pow(cpep, vDAHPS_nDAHPSpep)*vDAHPS_rmaxDAHPS/((pow(ce4p, vDAHPS_nDAHPSe4p) + vDAHPS_KDAHPSe4p)*(pow(cpep, vDAHPS_nDAHPSpep) + vDAHPS_KDAHPSpep));
    w[2] = 1.0*cdhap*vDHAP_mu;
    w[3] = 1.0*ce4p*vE4P_mu;
    w[4] = 1.0*vENO_rmaxENO*(-cpep/vENO_KENOeq + cpg2)/(cpg2 + vENO_KENOpg2*(cpep/vENO_KENOpep + 1));
    w[5] = 1.0*vEXTER_Dil*(-cglcex + vEXTER_cfeed);
    w[6] = 1.0*cg1p*vG1PAT_rmaxG1PAT*(pow(cfdp/vG1PAT_KG1PATfdp, vG1PAT_nG1PATfdp) + 1)*(-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/((cg1p + vG1PAT_KG1PATg1p)*(-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + vG1PAT_KG1PATatp + 4.2699999999999996));
    w[7] = 1.0*cdhap*vG3PDH_rmaxG3PDH/(cdhap + vG3PDH_KG3PDHdhap);
    w[8] = 1.0*cg6p*vG6P_mu;
    w[9] = 1.0*cg6p*vG6PDH_rmaxG6PDH*(-0.0055399999999999998*t/(0.01*pow(t, 2) - 0.27100000000000002*t + 2.7999999999999998) + 0.159 + 0.182/(0.52600000000000002*t + 4.8200000000000003))/((1 + (0.062 + 0.33200000000000002*pow(2.718, -0.46400000000000002*t)*(1.3619999999999999e-13*pow(t, 11) + 0.0166*pow(t, 1.5800000000000001) + 0.000166*pow(t, 4.7300000000000004) + 1.3120000000000002e-10*pow(t, 7.8899999999999997) + 1.2330000000000003e-16*pow(t, 14.199999999999999)))/vG6PDH_KG6PDHnadphg6pinh)*(cg6p + vG6PDH_KG6PDHg6p)*(-0.0055399999999999998*t/(0.01*pow(t, 2) - 0.27100000000000002*t + 2.7999999999999998) + vG6PDH_KG6PDHnadp*(1 + (0.062 + 0.33200000000000002*pow(2.718, -0.46400000000000002*t)*(1.3619999999999999e-13*pow(t, 11) + 0.0166*pow(t, 1.5800000000000001) + 0.000166*pow(t, 4.7300000000000004) + 1.3120000000000002e-10*pow(t, 7.8899999999999997) + 1.2330000000000003e-16*pow(t, 14.199999999999999)))/vG6PDH_KG6PDHnadphnadpinh) + 0.159 + 0.182/(0.52600000000000002*t + 4.8200000000000003)));
    w[10] = 1.0*cgap*vGAP_mu;
    w[11] = 1.0*vGAPDH_rmaxGAPDH*(cgap*(1.3140000000000001*pow(2.73, -0.043499999999999997*t - 0.34200000000000003) - pow(2.73, -0.0218*t - 0.17100000000000001)*(t + 7.8710000000000004)/(t + 8.4809999999999999) + 1.3140000000000001) - cpgp*(0.093399999999999997 + 0.0011100000000000001*pow(2.371, -0.123*t)*(0.104*pow(t, 3) + 0.84399999999999997*t))/vGAPDH_KGAPDHeq)/((cgap + vGAPDH_KGAPDHgap*(cpgp/vGAPDH_KGAPDHpgp + 1))*(1.3140000000000001*pow(2.73, -0.043499999999999997*t - 0.34200000000000003) - pow(2.73, -0.0218*t - 0.17100000000000001)*(t + 7.8710000000000004)/(t + 8.4809999999999999) + vGAPDH_KGAPDHnad*(1 + (0.093399999999999997 + 0.0011100000000000001*pow(2.371, -0.123*t)*(0.104*pow(t, 3) + 0.84399999999999997*t))/vGAPDH_KGAPDHnadh) + 1.3140000000000001));
    w[12] = 1.0*cg1p*vGLP_mu;
    w[13] = 1.0*vMURSyNTH_rmaxMurSynth;
    w[14] = 1.0*vMethSynth_rmaxMetSynth;
    w[15] = 1.0*pow(cpyr, vPDH_nPDH)*vPDH_rmaxPDH/(pow(cpyr, vPDH_nPDH) + vPDH_KPDHpyr);
    w[16] = 1.0*cpep*vPEP_mu;
    w[17] = 1.0*cf6p*vPFK_rmaxPFK*(-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/((cf6p + vPFK_KPFKf6ps*(cpep/vPFK_KPFKpep + 1 + (7.25*t/(0.17000000000000001*pow(t, 2) + 1.47*t + 7.25) + 0.123 + 1.073/(8.0500000000000007*t + 1.29))/vPFK_KPFKampb + (0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))/vPFK_KPFKadpb)/(1 + (7.25*t/(0.17000000000000001*pow(t, 2) + 1.47*t + 7.25) + 0.123 + 1.073/(8.0500000000000007*t + 1.29))/vPFK_KPFKampa + (0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))/vPFK_KPFKadpa))*(vPFK_LPFK*pow(cf6p*(1 + (7.25*t/(0.17000000000000001*pow(t, 2) + 1.47*t + 7.25) + 0.123 + 1.073/(8.0500000000000007*t + 1.29))/vPFK_KPFKampa + (0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))/vPFK_KPFKadpa)/(vPFK_KPFKf6ps*(cpep/vPFK_KPFKpep + 1 + (7.25*t/(0.17000000000000001*pow(t, 2) + 1.47*t + 7.25) + 0.123 + 1.073/(8.0500000000000007*t + 1.29))/vPFK_KPFKampb + (0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))/vPFK_KPFKadpb)) + 1, -vPFK_nPFK) + 1)*(-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + vPFK_KPFKatps*(1 + (0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))/vPFK_KPFKadpc) + 4.2699999999999996));
    w[18] = 1.0*cpg*vPG_mu;
    w[19] = 1.0*cpg3*vPG3_mu;
    w[20] = 1.0*cpg*vPGDH_rmaxPGDH*(-0.0055399999999999998*t/(0.01*pow(t, 2) - 0.27100000000000002*t + 2.7999999999999998) + 0.159 + 0.182/(0.52600000000000002*t + 4.8200000000000003))/((cpg + vPGDH_KPGDHpg)*(-0.0055399999999999998*t/(0.01*pow(t, 2) - 0.27100000000000002*t + 2.7999999999999998) + vPGDH_KPGDHnadp*(1 + (-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/vPGDH_KPGDHatpinh)*(1 + (0.062 + 0.33200000000000002*pow(2.718, -0.46400000000000002*t)*(1.3619999999999999e-13*pow(t, 11) + 0.0166*pow(t, 1.5800000000000001) + 0.000166*pow(t, 4.7300000000000004) + 1.3120000000000002e-10*pow(t, 7.8899999999999997) + 1.2330000000000003e-16*pow(t, 14.199999999999999)))/vPGDH_KPGDHnadphinh) + 0.159 + 0.182/(0.52600000000000002*t + 4.8200000000000003)));
    w[21] = 1.0*vPGI_rmaxPGI*(-cf6p/vPGI_KPGIeq + cg6p)/(cg6p + vPGI_KPGIg6p*(cf6p/(vPGI_KPGIf6p*(cpg/vPGI_KPGIf6ppginh + 1)) + cpg/vPGI_KPGIg6ppginh + 1));
    w[22] = 1.0*vPGK_rmaxPGK*(-cpg3*(-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/vPGK_KPGKeq + cpgp*(0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t)))/((cpgp + vPGK_KPGKpgp*(cpg3/vPGK_KPGKpg3 + 1))*(vPGK_KPGKadp*(1 + (-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/vPGK_KPGKatp) + 0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t)));
    w[23] = 1.0*vPGM_rmaxPGM*(-cg1p/vPGM_KPGMeq + cg6p)/(cg6p + vPGM_KPGMg6p*(cg1p/vPGM_KPGMg1p + 1));
    w[24] = 1.0*cpgp*vPGP_mu;
    w[25] = 1.0*cpep*vPK_rmaxPK*(0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t))*pow(cpep/vPK_KPKpep + 1, vPK_nPK - 1)/(vPK_KPKpep*(vPK_LPK*pow((1 + (-4.1630000000000003*t/(0.036400000000000002*pow(t, 2) + 1.4299999999999999*t + 0.65700000000000003) + 4.2699999999999996)/vPK_KPKatp)/(cfdp/vPK_KPKfdp + 1 + (7.25*t/(0.17000000000000001*pow(t, 2) + 1.47*t + 7.25) + 0.123 + 1.073/(8.0500000000000007*t + 1.29))/vPK_KPKamp), vPK_nPK) + pow(cpep/vPK_KPKpep + 1, vPK_nPK))*(vPK_KPKadp + 0.58199999999999996 + 1.73*pow(2.7309999999999999, -0.14999999999999999*t)*(0.000214*pow(t, 3) + 0.12*t)));
    w[26] = 1.0*crib5p*vPPK_rmaxRPPK/(crib5p + vPPK_KRPPKrib5p);
    w[27] = 1.0*cglcex*cpep*vPTS_rmaxPTS/(cpyr*(pow(cg6p, vPTS_nPTSg6p)/vPTS_KPTSg6p + 1)*(cglcex*cpep/cpyr + cglcex*vPTS_KPTSa3 + cpep*vPTS_KPTSa2/cpyr + vPTS_KPTSa1));
    w[28] = 1.0*vR5PI_rmaxR5PI*(-crib5p/vR5PI_KR5PIeq + cribu5p);
    w[29] = 1.0*crib5p*vRIB5P_mu;
    w[30] = 1.0*cribu5p*vRibu5p_mu;
    w[31] = 1.0*vRu5P_rmaxRu5P*(cribu5p - cxyl5p/vRu5P_KRu5Peq);
    w[32] = 1.0*csed7p*vSED7P_mu;
    w[33] = 1.0*cpep*vSynth1_rmaxSynth1/(cpep + vSynth1_KSynth1pep);
    w[34] = 1.0*cpyr*vSynth2_rmaxSynth2/(cpyr + vSynth2_KSynth2pyr);
    w[35] = 1.0*vTA_rmaxTA*(-ce4p*cf6p/vTA_KTAeq + cgap*csed7p);
    w[36] = 1.0*vTIS_rmaxTIS*(cdhap - cgap/vTIS_kTISeq)/(cdhap + vTIS_kTISdhap*(cgap/vTIS_kTISgap + 1));
    w[37] = 1.0*vTKA_rmaxTKa*(-cgap*csed7p/vTKA_KTKaeq + crib5p*cxyl5p);
    w[38] = 1.0*vTKB_rmaxTKb*(ce4p*cxyl5p - cf6p*cgap/vTKB_KTKbeq);
    w[39] = 1.0*vTRPSYNTH_rmaxTrpSynth;
    w[40] = 1.0*cxyl5p*vXYL5P_mu;
    w[41] = 1.0*cf6p*vf6P_mu;
    w[42] = 1.0*cfdp*vfdP_mu;
    w[43] = 1.0*cpep*vpepCxylase_rmaxpepCxylase*(pow(cfdp/vpepCxylase_KpepCxylasefdp, vpepCxylase_npepCxylasefdp) + 1)/(cpep + vpepCxylase_KpepCxylasepep);
    w[44] = 1.0*cpg2*vpg2_mu;
    w[45] = 1.0*cpyr*vpyr_mu;
    w[46] = 1.0*vrpGluMu_rmaxPGluMu*(-cpg2/vrpGluMu_KPGluMueq + cpg3)/(cpg3 + vrpGluMu_KPGluMupg3*(cpg2/vrpGluMu_KPGluMupg2 + 1));
    w[47] = 1.0*cpg3*vsersynth_rmaxSerSynth/(cpg3 + vsersynth_KSerSynthpg3);
}