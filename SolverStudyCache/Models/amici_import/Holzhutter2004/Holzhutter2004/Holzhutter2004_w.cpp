#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_Holzhutter2004(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*vGLT_Vmaxv0*(-Glcin/vGLT_Keqv0 + Glcout)/(vGLT_KMoutv0*(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1));
    w[1] = 1.0*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
    w[2] = 1.0*vGPI_Vmaxv2*(-Fru6P/vGPI_Keqv2 + Glc6P)/(Glc6P + vGPI_KGlc6Pv2*(Fru6P/vGPI_KFru6Pv2 + 1));
    w[3] = 1.0*vPFK_Vmaxv3*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/((Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1));
    w[4] = 1.0*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(vALD_KFru16P2v4*(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1));
    w[5] = 1.0*vTPI_Vmaxv5*(DHAP - GraP/vTPI_Keqv5)/(DHAP + vTPI_KDHAPv5*(GraP/vTPI_KGraPv5 + 1));
    w[6] = 1.0*vGAPDH_Vmaxv6*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*vGAPDH_KNADv6*vGAPDH_KPv6*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
    w[7] = 1.0*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)/(vPGK_K13P2Gv7*vPGK_KMgADPv7*((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1));
    w[8] = 1.0*vBPGM_kDPGMv8*(Gri13P2 - (Gri23P2f + MgGri23P2)/vBPGM_Keqv8)/(1 + (Gri23P2f + MgGri23P2)/vBPGM_K23P2Gv8);
    w[9] = 1.0*vBPGP_Vmaxv9*(Gri23P2f - Gri3P/vBPGP_Keqv9 + MgGri23P2)/(Gri23P2f + MgGri23P2 + vBPGP_K23P2Gv9);
    w[10] = 1.0*vPGM_Vmaxv10*(-Gri2P/vPGM_Keqv10 + Gri3P)/(Gri3P + vPGM_K3PGv10*(Gri2P/vPGM_K2PGv10 + 1));
    w[11] = 1.0*vENO_Vmaxv11*(Gri2P - PEP/vENO_Keqv11)/(Gri2P + vENO_K2PGv11*(PEP/vENO_KPEPv11 + 1));
    w[12] = 1.0*vPK_Vmaxv12*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/((MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1));
    w[13] = 1.0*vLDHNADH_Vmaxv13*(-Lac*NAD/vLDHNADH_Keqv13 + NADH*Pyr);
    w[14] = 1.0*vLDHNADPH_kLDHv14*(-Lac*NADPf/vLDHNADPH_Keqv14 + NADPHf*Pyr);
    w[15] = 1.0*MgATP*vATPase_kATPasev15;
    w[16] = 1.0*vAK_Vmaxv16*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(vAK_KAMPv16*vAK_KATPv16*(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16));
    w[17] = 1.0*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(vG6PDH_KG6Pv17*vG6PDH_KNADPv17*(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17));
    w[18] = 1.0*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(vPGLDH_K6PG1v18*vPGLDH_KNADPv18*(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18));
    w[19] = 1.0*vGSSGRD_Vmaxv19*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1);
    w[20] = 1.0*GSH*vGSHox_Kv20;
    w[21] = 1.0*vRibPepi_Vmaxv21*(Rul5P - Xul5P/vRibPepi_Keqv21)/(Rul5P + vRibPepi_KRu5Pv21*(Xul5P/vRibPepi_KX5Pv21 + 1));
    w[22] = 1.0*vRibPiso_Vmaxv22*(-Rib5P/vRibPiso_Keqv22 + Rul5P)/(Rul5P + vRibPiso_KRu5Pv22*(Rib5P/vRibPiso_KR5Pv22 + 1));
    w[23] = 1.0*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23));
    w[24] = 1.0*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24));
    w[25] = 1.0*vPPRPPS_Vmaxv25*(-MgAMP*PRPP/vPPRPPS_Keqv25 + MgATP*Rib5P)/((MgATP + vPPRPPS_KATPv25)*(Rib5P + vPPRPPS_KR5Pv25));
    w[26] = 1.0*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26));
    w[27] = 1.0*vPhiexch_Vmaxv27*(-Phi/vPhiexch_Keqv27 + Phiex);
    w[28] = 1.0*vLacexch_Vmaxv28*(-Lac/vLacexch_Keqv28 + Lacex);
    w[29] = 1.0*vPyrexch_Vmaxv29*(-Pyr/vPyrexch_Keqv29 + Pyrex);
    w[30] = 1.0*vMgATP_EqMult*(-ATPf*Mgf/vMgATP_KdATP + MgATP);
    w[31] = 1.0*vMgADP_EqMult*(-ADPf*Mgf/vMgADP_KdADP + MgADP);
    w[32] = 1.0*vMgAMP_EqMult*(-AMPf*Mgf/vMgAMP_KdAMP + MgAMP);
    w[33] = 1.0*vMgGri23P2_EqMult*(-Gri23P2f*Mgf/vMgGri23P2_Kd23P2G + MgGri23P2);
    w[34] = 1.0*vP1NADP_EqMult*(-NADPf*P1f/vP1NADP_Kd1 + P1NADP);
    w[35] = 1.0*vP1NADPH_EqMult*(-NADPHf*P1f/vP1NADPH_Kd3 + P1NADPH);
    w[36] = 1.0*vP2NADP_EqMult*(-NADPf*P2f/vP2NADP_Kd2 + P2NADP);
    w[37] = 1.0*vP2NADPH_EqMult*(-NADPHf*P2f/vP2NADPH_Kd4 + P2NADPH);
}