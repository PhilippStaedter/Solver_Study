#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Holzhutter2004(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*Glcin*Glcout*vGLT_Vmaxv0*(-Glcin/vGLT_Keqv0 + Glcout)/(vGLT_KMinv0*pow(vGLT_KMoutv0, 2)*pow(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1, 2));
            break;
        case 1:
            dwdp[0] = 1.0*vGLT_Vmaxv0*(-Glcin/vGLT_Keqv0 + Glcout)*(Glcin*Glcout*vGLT_alfav0/(pow(vGLT_KMinv0, 2)*vGLT_KMoutv0) + Glcin/pow(vGLT_KMinv0, 2))/(vGLT_KMoutv0*pow(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1, 2));
            break;
        case 2:
            dwdp[0] = 1.0*Glcin*vGLT_Vmaxv0/(vGLT_KMoutv0*pow(vGLT_Keqv0, 2)*(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1));
            break;
        case 3:
            dwdp[0] = 1.0*vGLT_Vmaxv0*(-Glcin/vGLT_Keqv0 + Glcout)*(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*pow(vGLT_KMoutv0, 2)) + Glcout/pow(vGLT_KMoutv0, 2))/(vGLT_KMoutv0*pow(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1, 2)) - 1.0*vGLT_Vmaxv0*(-Glcin/vGLT_Keqv0 + Glcout)/(pow(vGLT_KMoutv0, 2)*(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1));
            break;
        case 4:
            dwdp[0] = 1.0*(-Glcin/vGLT_Keqv0 + Glcout)/(vGLT_KMoutv0*(Glcin*Glcout*vGLT_alfav0/(vGLT_KMinv0*vGLT_KMoutv0) + Glcin/vGLT_KMinv0 + Glcout/vGLT_KMoutv0 + 1));
            break;
        case 5:
            dwdp[1] = 1.0*Glcin*Mgf*vHEX_Inhibv1*vHEX_Vmax1v1*(Gri23P2f + MgGri23P2)*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_KMg23P2Gv1, 2)*vHEX_KMgATPv1*vHEX_KMgv1*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2));
            break;
        case 6:
            dwdp[1] = 1.0*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(Gri23P2f + MgGri23P2)*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_K23P2Gv1, 2)*vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2));
            break;
        case 7:
            dwdp[1] = 1.0*Glc6P*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(Mgf/vHEX_KMgv1 + 1)*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_KGlc6Pv1, 2)*vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2));
            break;
        case 8:
            dwdp[1] = 1.0*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)*(Mgf*(Glc6P/vHEX_KGlc6Pv1 + 1.55)/pow(vHEX_KMgv1, 2) + Mgf/pow(vHEX_KMgv1, 2) + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*pow(vHEX_KMgv1, 2)))/(vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2));
            break;
        case 9:
            dwdp[1] = 1.0*Glc6P*Glcin*MgADP*vHEX_Inhibv1*vHEX_Vmax1v1/(vHEX_KMgATPv1*pow(vHEX_Keqv1, 2)*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 10:
            dwdp[1] = -1.0*Glcin*MgATP*Mgf*vHEX_Inhibv1*vHEX_Vmax2v1/(pow(vHEX_KMgATPMgv1, 2)*vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1)) + 1.0*Glcin*MgATP*Mgf*vHEX_Inhibv1*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_KMgATPMgv1, 2)*pow(vHEX_KMgATPv1, 2)*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2));
            break;
        case 11:
            dwdp[1] = 1.0*Glcin*MgATP*Mgf*vHEX_Inhibv1/(vHEX_KMgATPMgv1*vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 12:
            dwdp[1] = 1.0*Glcin*MgATP*vHEX_Inhibv1*vHEX_Vmax1v1*(Mgf/vHEX_KMgATPMgv1 + 1)*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_KMgATPv1, 3)*(Glcin + vHEX_KMGlcv1)*pow(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1, 2)) - 1.0*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(pow(vHEX_KMgATPv1, 2)*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 13:
            dwdp[1] = -1.0*Glcin*MgATP*Mgf*vHEX_Inhibv1*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_KMgATPv1*vHEX_Vmax1v1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1)) + 1.0*Glcin*vHEX_Inhibv1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 14:
            dwdp[1] = -1.0*Glcin*vHEX_Inhibv1*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(vHEX_KMgATPv1*pow(Glcin + vHEX_KMGlcv1, 2)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 15:
            dwdp[1] = 1.0*Glcin*vHEX_Vmax1v1*(-Glc6P*MgADP/vHEX_Keqv1 + MgATP*Mgf*vHEX_Vmax2v1/(vHEX_KMgATPMgv1*vHEX_Vmax1v1) + MgATP)/(vHEX_KMgATPv1*(Glcin + vHEX_KMGlcv1)*(MgATP*(Mgf/vHEX_KMgATPMgv1 + 1)/vHEX_KMgATPv1 + Mgf/vHEX_KMgv1 + Mgf*(Gri23P2f + MgGri23P2)/(vHEX_KMg23P2Gv1*vHEX_KMgv1) + (Glc6P/vHEX_KGlc6Pv1 + 1.55)*(Mgf/vHEX_KMgv1 + 1) + 1 + (Gri23P2f + MgGri23P2)/vHEX_K23P2Gv1));
            break;
        case 16:
            dwdp[2] = 1.0*Fru6P*vGPI_KGlc6Pv2*vGPI_Vmaxv2*(-Fru6P/vGPI_Keqv2 + Glc6P)/(pow(vGPI_KFru6Pv2, 2)*pow(Glc6P + vGPI_KGlc6Pv2*(Fru6P/vGPI_KFru6Pv2 + 1), 2));
            break;
        case 17:
            dwdp[2] = 1.0*vGPI_Vmaxv2*(-Fru6P/vGPI_KFru6Pv2 - 1)*(-Fru6P/vGPI_Keqv2 + Glc6P)/pow(Glc6P + vGPI_KGlc6Pv2*(Fru6P/vGPI_KFru6Pv2 + 1), 2);
            break;
        case 18:
            dwdp[2] = 1.0*Fru6P*vGPI_Vmaxv2/(pow(vGPI_Keqv2, 2)*(Glc6P + vGPI_KGlc6Pv2*(Fru6P/vGPI_KFru6Pv2 + 1)));
            break;
        case 19:
            dwdp[2] = 1.0*(-Fru6P/vGPI_Keqv2 + Glc6P)/(Glc6P + vGPI_KGlc6Pv2*(Fru6P/vGPI_KFru6Pv2 + 1));
            break;
        case 20:
            dwdp[3] = -4.0*vPFK_L0v3*vPFK_Vmaxv3*(AMPf + MgAMP)*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(vPFK_KAMPv3, 2)*pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 5)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)*pow(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1, 2));
            break;
        case 21:
            dwdp[3] = 4.0*Mgf*vPFK_L0v3*vPFK_Vmaxv3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 3)*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(vPFK_KMgv3, 2)*pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)*pow(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1, 2));
            break;
        case 22:
            dwdp[3] = 4.0*ATPf*vPFK_L0v3*vPFK_Vmaxv3*pow(ATPf/vPFK_KATPv3 + 1, 3)*pow(Mgf/vPFK_KMgv3 + 1, 4)*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(vPFK_KATPv3, 2)*pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)*pow(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1, 2));
            break;
        case 23:
            dwdp[3] = -1.0*vPFK_Vmaxv3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)*pow(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1, 2));
            break;
        case 24:
            dwdp[3] = -1.0*vPFK_Vmaxv3*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/((Fru6P + vPFK_KFru6Pv3)*pow(MgATP + vPFK_KMgATPv3, 2)*(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1));
            break;
        case 25:
            dwdp[3] = -4.0*Fru6P*vPFK_L0v3*vPFK_Vmaxv3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(vPFK_KFru6Pv3, 2)*pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 5)*pow(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1, 2)) - 1.0*vPFK_Vmaxv3*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/(pow(Fru6P + vPFK_KFru6Pv3, 2)*(MgATP + vPFK_KMgATPv3)*(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1));
            break;
        case 26:
            dwdp[3] = 1.0*Fru16P2*MgADP*vPFK_Vmaxv3/(pow(vPFK_Keqv3, 2)*(Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1));
            break;
        case 27:
            dwdp[3] = 1.0*(-Fru16P2*MgADP/vPFK_Keqv3 + Fru6P*MgATP)/((Fru6P + vPFK_KFru6Pv3)*(MgATP + vPFK_KMgATPv3)*(vPFK_L0v3*pow(ATPf/vPFK_KATPv3 + 1, 4)*pow(Mgf/vPFK_KMgv3 + 1, 4)/(pow(1 + (AMPf + MgAMP)/vPFK_KAMPv3, 4)*pow(Fru6P/vPFK_KFru6Pv3 + 1, 4)) + 1));
            break;
        case 28:
            dwdp[4] = 1.0*Fru16P2*GraP*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(pow(vALD_KFru16P2v4, 2)*pow(vALD_KiiGraPv4, 2)*pow(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1, 2));
            break;
        case 29:
            dwdp[4] = 1.0*DHAP*vALD_Vmaxv4*(GraP + vALD_KGraPv4)*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(pow(vALD_KDHAPv4, 2)*vALD_KFru16P2v4*vALD_KiGraPv4*pow(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1, 2));
            break;
        case 30:
            dwdp[4] = -1.0*DHAP*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(vALD_KDHAPv4*vALD_KFru16P2v4*vALD_KiGraPv4*pow(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1, 2));
            break;
        case 31:
            dwdp[4] = 1.0*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)*(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*pow(vALD_KiGraPv4, 2)) + GraP/pow(vALD_KiGraPv4, 2))/(vALD_KFru16P2v4*pow(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1, 2));
            break;
        case 32:
            dwdp[4] = 1.0*DHAP*GraP*vALD_Vmaxv4/(vALD_KFru16P2v4*pow(vALD_Keqv4, 2)*(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1));
            break;
        case 33:
            dwdp[4] = 1.0*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)*(Fru16P2*GraP/(pow(vALD_KFru16P2v4, 2)*vALD_KiiGraPv4) + Fru16P2/pow(vALD_KFru16P2v4, 2))/(vALD_KFru16P2v4*pow(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1, 2)) - 1.0*vALD_Vmaxv4*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(pow(vALD_KFru16P2v4, 2)*(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1));
            break;
        case 34:
            dwdp[4] = 1.0*(-DHAP*GraP/vALD_Keqv4 + Fru16P2)/(vALD_KFru16P2v4*(DHAP*(GraP + vALD_KGraPv4)/(vALD_KDHAPv4*vALD_KiGraPv4) + Fru16P2*GraP/(vALD_KFru16P2v4*vALD_KiiGraPv4) + Fru16P2/vALD_KFru16P2v4 + GraP/vALD_KiGraPv4 + 1));
            break;
        case 35:
            dwdp[5] = 1.0*GraP*vTPI_KDHAPv5*vTPI_Vmaxv5*(DHAP - GraP/vTPI_Keqv5)/(pow(vTPI_KGraPv5, 2)*pow(DHAP + vTPI_KDHAPv5*(GraP/vTPI_KGraPv5 + 1), 2));
            break;
        case 36:
            dwdp[5] = 1.0*vTPI_Vmaxv5*(DHAP - GraP/vTPI_Keqv5)*(-GraP/vTPI_KGraPv5 - 1)/pow(DHAP + vTPI_KDHAPv5*(GraP/vTPI_KGraPv5 + 1), 2);
            break;
        case 37:
            dwdp[5] = 1.0*GraP*vTPI_Vmaxv5/(pow(vTPI_Keqv5, 2)*(DHAP + vTPI_KDHAPv5*(GraP/vTPI_KGraPv5 + 1)));
            break;
        case 38:
            dwdp[5] = 1.0*(DHAP - GraP/vTPI_Keqv5)/(DHAP + vTPI_KDHAPv5*(GraP/vTPI_KGraPv5 + 1));
            break;
        case 39:
            dwdp[6] = 1.0*Gri13P2*vGAPDH_Vmaxv6*(NADH/vGAPDH_KNADHv6 + 1)*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(pow(vGAPDH_K13P2Gv6, 2)*vGAPDH_KGraPv6*vGAPDH_KNADv6*vGAPDH_KPv6*pow((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1, 2));
            break;
        case 40:
            dwdp[6] = 1.0*NADH*vGAPDH_Vmaxv6*(Gri13P2/vGAPDH_K13P2Gv6 + 1)*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*pow(vGAPDH_KNADHv6, 2)*vGAPDH_KNADv6*vGAPDH_KPv6*pow((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1, 2));
            break;
        case 41:
            dwdp[6] = 1.0*Gri13P2*NADH*vGAPDH_Vmaxv6/(vGAPDH_KGraPv6*vGAPDH_KNADv6*vGAPDH_KPv6*pow(vGAPDH_Keqv6, 2)*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
            break;
        case 42:
            dwdp[6] = 1.0*Phi*vGAPDH_Vmaxv6*(GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*vGAPDH_KNADv6*pow(vGAPDH_KPv6, 3)*pow((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1, 2)) - 1.0*vGAPDH_Vmaxv6*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*vGAPDH_KNADv6*pow(vGAPDH_KPv6, 2)*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
            break;
        case 43:
            dwdp[6] = 1.0*GraP*vGAPDH_Vmaxv6*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1)*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(pow(vGAPDH_KGraPv6, 3)*vGAPDH_KNADv6*vGAPDH_KPv6*pow((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1, 2)) - 1.0*vGAPDH_Vmaxv6*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(pow(vGAPDH_KGraPv6, 2)*vGAPDH_KNADv6*vGAPDH_KPv6*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
            break;
        case 44:
            dwdp[6] = 1.0*NAD*vGAPDH_Vmaxv6*(GraP/vGAPDH_KGraPv6 + 1)*(Phi/vGAPDH_KPv6 + 1)*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*pow(vGAPDH_KNADv6, 3)*vGAPDH_KPv6*pow((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1, 2)) - 1.0*vGAPDH_Vmaxv6*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*pow(vGAPDH_KNADv6, 2)*vGAPDH_KPv6*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
            break;
        case 45:
            dwdp[6] = 1.0*(GraP*NAD*Phi - Gri13P2*NADH/vGAPDH_Keqv6)/(vGAPDH_KGraPv6*vGAPDH_KNADv6*vGAPDH_KPv6*((GraP/vGAPDH_KGraPv6 + 1)*(NAD/vGAPDH_KNADv6 + 1)*(Phi/vGAPDH_KPv6 + 1) + (Gri13P2/vGAPDH_K13P2Gv6 + 1)*(NADH/vGAPDH_KNADHv6 + 1) - 1));
            break;
        case 46:
            dwdp[7] = 1.0*Gri3P*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)*(MgATP/vPGK_KMgATPv7 + 1)/(vPGK_K13P2Gv7*pow(vPGK_K3PGv7, 2)*vPGK_KMgADPv7*pow((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1, 2));
            break;
        case 47:
            dwdp[7] = 1.0*MgATP*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)*(Gri3P/vPGK_K3PGv7 + 1)/(vPGK_K13P2Gv7*vPGK_KMgADPv7*pow(vPGK_KMgATPv7, 2)*pow((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1, 2));
            break;
        case 48:
            dwdp[7] = 1.0*Gri3P*MgATP*vPGK_Vmaxv7/(vPGK_K13P2Gv7*vPGK_KMgADPv7*pow(vPGK_Keqv7, 2)*((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1));
            break;
        case 49:
            dwdp[7] = 1.0*Gri13P2*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)*(MgADP/vPGK_KMgADPv7 + 1)/(pow(vPGK_K13P2Gv7, 3)*vPGK_KMgADPv7*pow((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1, 2)) - 1.0*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)/(pow(vPGK_K13P2Gv7, 2)*vPGK_KMgADPv7*((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1));
            break;
        case 50:
            dwdp[7] = 1.0*MgADP*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)*(Gri13P2/vPGK_K13P2Gv7 + 1)/(vPGK_K13P2Gv7*pow(vPGK_KMgADPv7, 3)*pow((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1, 2)) - 1.0*vPGK_Vmaxv7*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)/(vPGK_K13P2Gv7*pow(vPGK_KMgADPv7, 2)*((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1));
            break;
        case 51:
            dwdp[7] = 1.0*(Gri13P2*MgADP - Gri3P*MgATP/vPGK_Keqv7)/(vPGK_K13P2Gv7*vPGK_KMgADPv7*((Gri13P2/vPGK_K13P2Gv7 + 1)*(MgADP/vPGK_KMgADPv7 + 1) + (Gri3P/vPGK_K3PGv7 + 1)*(MgATP/vPGK_KMgATPv7 + 1) - 1));
            break;
        case 52:
            dwdp[8] = 1.0*vBPGM_kDPGMv8*(Gri13P2 - (Gri23P2f + MgGri23P2)/vBPGM_Keqv8)*(Gri23P2f + MgGri23P2)/(pow(vBPGM_K23P2Gv8, 2)*pow(1 + (Gri23P2f + MgGri23P2)/vBPGM_K23P2Gv8, 2));
            break;
        case 53:
            dwdp[8] = -1.0*vBPGM_kDPGMv8*(-Gri23P2f - MgGri23P2)/(pow(vBPGM_Keqv8, 2)*(1 + (Gri23P2f + MgGri23P2)/vBPGM_K23P2Gv8));
            break;
        case 54:
            dwdp[8] = 1.0*(Gri13P2 - (Gri23P2f + MgGri23P2)/vBPGM_Keqv8)/(1 + (Gri23P2f + MgGri23P2)/vBPGM_K23P2Gv8);
            break;
        case 55:
            dwdp[9] = -1.0*vBPGP_Vmaxv9*(Gri23P2f - Gri3P/vBPGP_Keqv9 + MgGri23P2)/pow(Gri23P2f + MgGri23P2 + vBPGP_K23P2Gv9, 2);
            break;
        case 56:
            dwdp[9] = 1.0*Gri3P*vBPGP_Vmaxv9/(pow(vBPGP_Keqv9, 2)*(Gri23P2f + MgGri23P2 + vBPGP_K23P2Gv9));
            break;
        case 57:
            dwdp[9] = 1.0*(Gri23P2f - Gri3P/vBPGP_Keqv9 + MgGri23P2)/(Gri23P2f + MgGri23P2 + vBPGP_K23P2Gv9);
            break;
        case 58:
            dwdp[10] = 1.0*Gri2P*vPGM_K3PGv10*vPGM_Vmaxv10*(-Gri2P/vPGM_Keqv10 + Gri3P)/(pow(vPGM_K2PGv10, 2)*pow(Gri3P + vPGM_K3PGv10*(Gri2P/vPGM_K2PGv10 + 1), 2));
            break;
        case 59:
            dwdp[10] = 1.0*vPGM_Vmaxv10*(-Gri2P/vPGM_K2PGv10 - 1)*(-Gri2P/vPGM_Keqv10 + Gri3P)/pow(Gri3P + vPGM_K3PGv10*(Gri2P/vPGM_K2PGv10 + 1), 2);
            break;
        case 60:
            dwdp[10] = 1.0*Gri2P*vPGM_Vmaxv10/(pow(vPGM_Keqv10, 2)*(Gri3P + vPGM_K3PGv10*(Gri2P/vPGM_K2PGv10 + 1)));
            break;
        case 61:
            dwdp[10] = 1.0*(-Gri2P/vPGM_Keqv10 + Gri3P)/(Gri3P + vPGM_K3PGv10*(Gri2P/vPGM_K2PGv10 + 1));
            break;
        case 62:
            dwdp[11] = 1.0*PEP*vENO_K2PGv11*vENO_Vmaxv11*(Gri2P - PEP/vENO_Keqv11)/(pow(vENO_KPEPv11, 2)*pow(Gri2P + vENO_K2PGv11*(PEP/vENO_KPEPv11 + 1), 2));
            break;
        case 63:
            dwdp[11] = 1.0*vENO_Vmaxv11*(Gri2P - PEP/vENO_Keqv11)*(-PEP/vENO_KPEPv11 - 1)/pow(Gri2P + vENO_K2PGv11*(PEP/vENO_KPEPv11 + 1), 2);
            break;
        case 64:
            dwdp[11] = 1.0*PEP*vENO_Vmaxv11/(pow(vENO_Keqv11, 2)*(Gri2P + vENO_K2PGv11*(PEP/vENO_KPEPv11 + 1)));
            break;
        case 65:
            dwdp[11] = 1.0*(Gri2P - PEP/vENO_Keqv11)/(Gri2P + vENO_K2PGv11*(PEP/vENO_KPEPv11 + 1));
            break;
        case 66:
            dwdp[12] = -4.0*Fru16P2*vPK_L0v12*vPK_Vmaxv12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/(pow(vPK_KFru16P2v12, 2)*(MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*pow(Fru16P2/vPK_KFru16P2v12 + 1, 5)*pow(PEP/vPK_KPEPv12 + 1, 4)*pow(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1, 2));
            break;
        case 67:
            dwdp[12] = 4.0*vPK_L0v12*vPK_Vmaxv12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 3)*(ATPf + MgATP)*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/(pow(vPK_KATPv12, 2)*(MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)*pow(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1, 2));
            break;
        case 68:
            dwdp[12] = -1.0*vPK_Vmaxv12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/((MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)*pow(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1, 2));
            break;
        case 69:
            dwdp[12] = -1.0*vPK_Vmaxv12*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/(pow(MgADP + vPK_KMgADPv12, 2)*(PEP + vPK_KPEPv12)*(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1));
            break;
        case 70:
            dwdp[12] = -4.0*PEP*vPK_L0v12*vPK_Vmaxv12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/(pow(vPK_KPEPv12, 2)*(MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 5)*pow(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1, 2)) - 1.0*vPK_Vmaxv12*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/((MgADP + vPK_KMgADPv12)*pow(PEP + vPK_KPEPv12, 2)*(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1));
            break;
        case 71:
            dwdp[12] = 1.0*MgATP*Pyr*vPK_Vmaxv12/(pow(vPK_Keqv12, 2)*(MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1));
            break;
        case 72:
            dwdp[12] = 1.0*(MgADP*PEP - MgATP*Pyr/vPK_Keqv12)/((MgADP + vPK_KMgADPv12)*(PEP + vPK_KPEPv12)*(vPK_L0v12*pow(1 + (ATPf + MgATP)/vPK_KATPv12, 4)/(pow(Fru16P2/vPK_KFru16P2v12 + 1, 4)*pow(PEP/vPK_KPEPv12 + 1, 4)) + 1));
            break;
        case 73:
            dwdp[13] = 1.0*Lac*NAD*vLDHNADH_Vmaxv13/pow(vLDHNADH_Keqv13, 2);
            break;
        case 74:
            dwdp[13] = -1.0*Lac*NAD/vLDHNADH_Keqv13 + 1.0*NADH*Pyr;
            break;
        case 75:
            dwdp[14] = 1.0*Lac*NADPf*vLDHNADPH_kLDHv14/pow(vLDHNADPH_Keqv14, 2);
            break;
        case 76:
            dwdp[14] = -1.0*Lac*NADPf/vLDHNADPH_Keqv14 + 1.0*NADPHf*Pyr;
            break;
        case 77:
            dwdp[15] = 1.0*MgATP;
            break;
        case 78:
            dwdp[16] = 1.0*vAK_Vmaxv16*(2*ADPf*MgADP/pow(vAK_KADPv16, 3) + (ADPf + MgADP)/pow(vAK_KADPv16, 2))*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(vAK_KAMPv16*vAK_KATPv16*pow(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16, 2));
            break;
        case 79:
            dwdp[16] = 1.0*ADPf*MgADP*vAK_Vmaxv16/(vAK_KAMPv16*vAK_KATPv16*pow(vAK_Keqv16, 2)*(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16));
            break;
        case 80:
            dwdp[16] = 1.0*AMPf*vAK_Vmaxv16*(MgATP/vAK_KATPv16 + 1)*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(pow(vAK_KAMPv16, 3)*vAK_KATPv16*pow(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16, 2)) - 1.0*vAK_Vmaxv16*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(pow(vAK_KAMPv16, 2)*vAK_KATPv16*(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16));
            break;
        case 81:
            dwdp[16] = 1.0*MgATP*vAK_Vmaxv16*(AMPf/vAK_KAMPv16 + 1)*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(vAK_KAMPv16*pow(vAK_KATPv16, 3)*pow(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16, 2)) - 1.0*vAK_Vmaxv16*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(vAK_KAMPv16*pow(vAK_KATPv16, 2)*(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16));
            break;
        case 82:
            dwdp[16] = 1.0*(-ADPf*MgADP/vAK_Keqv16 + AMPf*MgATP)/(vAK_KAMPv16*vAK_KATPv16*(ADPf*MgADP/pow(vAK_KADPv16, 2) + (AMPf/vAK_KAMPv16 + 1)*(MgATP/vAK_KATPv16 + 1) + (ADPf + MgADP)/vAK_KADPv16));
            break;
        case 83:
            dwdp[17] = 1.0*vG6PDH_Vmaxv17*(Gri23P2f + MgGri23P2)*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(vG6PDH_KG6Pv17*vG6PDH_KNADPv17*pow(vG6PDH_KPGA23v17, 2)*pow(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17, 2));
            break;
        case 84:
            dwdp[17] = 1.0*NADPHf*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(vG6PDH_KG6Pv17*pow(vG6PDH_KNADPHv17, 2)*vG6PDH_KNADPv17*pow(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17, 2));
            break;
        case 85:
            dwdp[17] = 1.0*vG6PDH_Vmaxv17*(ATPf + MgATP)*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(pow(vG6PDH_KATPv17, 2)*vG6PDH_KG6Pv17*vG6PDH_KNADPv17*pow(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17, 2));
            break;
        case 86:
            dwdp[17] = 1.0*GlcA6P*NADPHf*vG6PDH_Vmaxv17/(vG6PDH_KG6Pv17*vG6PDH_KNADPv17*pow(vG6PDH_Keqv17, 2)*(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17));
            break;
        case 87:
            dwdp[17] = 1.0*NADPf*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)*(Glc6P/vG6PDH_KG6Pv17 + 1)/(vG6PDH_KG6Pv17*pow(vG6PDH_KNADPv17, 3)*pow(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17, 2)) - 1.0*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(vG6PDH_KG6Pv17*pow(vG6PDH_KNADPv17, 2)*(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17));
            break;
        case 88:
            dwdp[17] = 1.0*Glc6P*NADPf*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(pow(vG6PDH_KG6Pv17, 3)*pow(vG6PDH_KNADPv17, 2)*pow(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17, 2)) - 1.0*vG6PDH_Vmaxv17*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(pow(vG6PDH_KG6Pv17, 2)*vG6PDH_KNADPv17*(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17));
            break;
        case 89:
            dwdp[17] = 1.0*(Glc6P*NADPf - GlcA6P*NADPHf/vG6PDH_Keqv17)/(vG6PDH_KG6Pv17*vG6PDH_KNADPv17*(NADPHf/vG6PDH_KNADPHv17 + NADPf*(Glc6P/vG6PDH_KG6Pv17 + 1)/vG6PDH_KNADPv17 + 1 + (Gri23P2f + MgGri23P2)/vG6PDH_KPGA23v17 + (ATPf + MgATP)/vG6PDH_KATPv17));
            break;
        case 90:
            dwdp[18] = 1.0*NADPHf*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)*(GlcA6P/vPGLDH_K6PG2v18 + 1)/(vPGLDH_K6PG1v18*pow(vPGLDH_KNADPHv18, 2)*vPGLDH_KNADPv18*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2));
            break;
        case 91:
            dwdp[18] = 1.0*GlcA6P*NADPHf*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(vPGLDH_K6PG1v18*pow(vPGLDH_K6PG2v18, 2)*vPGLDH_KNADPHv18*vPGLDH_KNADPv18*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2));
            break;
        case 92:
            dwdp[18] = 1.0*vPGLDH_Vmaxv18*(ATPf + MgATP)*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(vPGLDH_K6PG1v18*pow(vPGLDH_KATPv18, 2)*vPGLDH_KNADPv18*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2));
            break;
        case 93:
            dwdp[18] = 1.0*vPGLDH_Vmaxv18*(Gri23P2f + MgGri23P2)*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)*(NADPf/vPGLDH_KNADPv18 + 1)/(vPGLDH_K6PG1v18*vPGLDH_KNADPv18*pow(vPGLDH_KPGA23v18, 2)*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2));
            break;
        case 94:
            dwdp[18] = 1.0*NADPHf*Rul5P*vPGLDH_Vmaxv18/(vPGLDH_K6PG1v18*vPGLDH_KNADPv18*pow(vPGLDH_Keqv18, 2)*(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18));
            break;
        case 95:
            dwdp[18] = 1.0*NADPf*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18)/(vPGLDH_K6PG1v18*pow(vPGLDH_KNADPv18, 3)*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2)) - 1.0*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(vPGLDH_K6PG1v18*pow(vPGLDH_KNADPv18, 2)*(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18));
            break;
        case 96:
            dwdp[18] = 1.0*GlcA6P*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)*(NADPf/vPGLDH_KNADPv18 + 1)/(pow(vPGLDH_K6PG1v18, 3)*vPGLDH_KNADPv18*pow(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18, 2)) - 1.0*vPGLDH_Vmaxv18*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(pow(vPGLDH_K6PG1v18, 2)*vPGLDH_KNADPv18*(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18));
            break;
        case 97:
            dwdp[18] = 1.0*(GlcA6P*NADPf - NADPHf*Rul5P/vPGLDH_Keqv18)/(vPGLDH_K6PG1v18*vPGLDH_KNADPv18*(NADPHf*(GlcA6P/vPGLDH_K6PG2v18 + 1)/vPGLDH_KNADPHv18 + (NADPf/vPGLDH_KNADPv18 + 1)*(GlcA6P/vPGLDH_K6PG1v18 + 1 + (Gri23P2f + MgGri23P2)/vPGLDH_KPGA23v18) + (ATPf + MgATP)/vPGLDH_KATPv18));
            break;
        case 98:
            dwdp[19] = 1.0*pow(GSH, 2)*NADPf*vGSSGRD_Vmaxv19/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*pow(vGSSGRD_Keqv19, 2)*(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1));
            break;
        case 99:
            dwdp[19] = 1.0*pow(GSH, 2)*NADPf*vGSSGRD_Vmaxv19/(pow(vGSSGRD_KGSHv19, 2)*pow(vGSSGRD_KNADPv19, 2)*vGSSGRD_Keqv19*(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1)) + 1.0*NADPf*vGSSGRD_Vmaxv19*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(pow(vGSSGRD_KNADPv19, 2)*pow(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1, 2));
            break;
        case 100:
            dwdp[19] = 2.0*pow(GSH, 2)*NADPf*vGSSGRD_Vmaxv19/(pow(vGSSGRD_KGSHv19, 3)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19*(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1)) - 1.0*NADPf*vGSSGRD_Vmaxv19*(-pow(GSH, 2)/pow(vGSSGRD_KGSHv19, 3) - GSH*(GSH/vGSSGRD_KGSHv19 + 1)/pow(vGSSGRD_KGSHv19, 2))*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(vGSSGRD_KNADPv19*pow(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1, 2));
            break;
        case 101:
            dwdp[19] = -1.0*GSSG*NADPHf*vGSSGRD_Vmaxv19/(vGSSGRD_KGSSGv19*pow(vGSSGRD_KNADPHv19, 2)*(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1)) + 1.0*NADPHf*vGSSGRD_Vmaxv19*(GSSG/vGSSGRD_KGSSGv19 + 1)*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(pow(vGSSGRD_KNADPHv19, 2)*pow(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1, 2));
            break;
        case 102:
            dwdp[19] = 1.0*GSSG*NADPHf*vGSSGRD_Vmaxv19*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(pow(vGSSGRD_KGSSGv19, 2)*vGSSGRD_KNADPHv19*pow(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1, 2)) - 1.0*GSSG*NADPHf*vGSSGRD_Vmaxv19/(pow(vGSSGRD_KGSSGv19, 2)*vGSSGRD_KNADPHv19*(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1));
            break;
        case 103:
            dwdp[19] = 1.0*(-pow(GSH, 2)*NADPf/(pow(vGSSGRD_KGSHv19, 2)*vGSSGRD_KNADPv19*vGSSGRD_Keqv19) + GSSG*NADPHf/(vGSSGRD_KGSSGv19*vGSSGRD_KNADPHv19))/(NADPHf*(GSSG/vGSSGRD_KGSSGv19 + 1)/vGSSGRD_KNADPHv19 + NADPf*(GSH*(GSH/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KGSHv19 + 1)/vGSSGRD_KNADPv19 + 1);
            break;
        case 104:
            dwdp[20] = 1.0*GSH;
            break;
        case 105:
            dwdp[21] = 1.0*Xul5P*vRibPepi_KRu5Pv21*vRibPepi_Vmaxv21*(Rul5P - Xul5P/vRibPepi_Keqv21)/(pow(vRibPepi_KX5Pv21, 2)*pow(Rul5P + vRibPepi_KRu5Pv21*(Xul5P/vRibPepi_KX5Pv21 + 1), 2));
            break;
        case 106:
            dwdp[21] = 1.0*vRibPepi_Vmaxv21*(Rul5P - Xul5P/vRibPepi_Keqv21)*(-Xul5P/vRibPepi_KX5Pv21 - 1)/pow(Rul5P + vRibPepi_KRu5Pv21*(Xul5P/vRibPepi_KX5Pv21 + 1), 2);
            break;
        case 107:
            dwdp[21] = 1.0*Xul5P*vRibPepi_Vmaxv21/(pow(vRibPepi_Keqv21, 2)*(Rul5P + vRibPepi_KRu5Pv21*(Xul5P/vRibPepi_KX5Pv21 + 1)));
            break;
        case 108:
            dwdp[21] = 1.0*(Rul5P - Xul5P/vRibPepi_Keqv21)/(Rul5P + vRibPepi_KRu5Pv21*(Xul5P/vRibPepi_KX5Pv21 + 1));
            break;
        case 109:
            dwdp[22] = 1.0*Rib5P*vRibPiso_KRu5Pv22*vRibPiso_Vmaxv22*(-Rib5P/vRibPiso_Keqv22 + Rul5P)/(pow(vRibPiso_KR5Pv22, 2)*pow(Rul5P + vRibPiso_KRu5Pv22*(Rib5P/vRibPiso_KR5Pv22 + 1), 2));
            break;
        case 110:
            dwdp[22] = 1.0*vRibPiso_Vmaxv22*(-Rib5P/vRibPiso_KR5Pv22 - 1)*(-Rib5P/vRibPiso_Keqv22 + Rul5P)/pow(Rul5P + vRibPiso_KRu5Pv22*(Rib5P/vRibPiso_KR5Pv22 + 1), 2);
            break;
        case 111:
            dwdp[22] = 1.0*Rib5P*vRibPiso_Vmaxv22/(pow(vRibPiso_Keqv22, 2)*(Rul5P + vRibPiso_KRu5Pv22*(Rib5P/vRibPiso_KR5Pv22 + 1)));
            break;
        case 112:
            dwdp[22] = 1.0*(-Rib5P/vRibPiso_Keqv22 + Rul5P)/(Rul5P + vRibPiso_KRu5Pv22*(Rib5P/vRibPiso_KR5Pv22 + 1));
            break;
        case 113:
            dwdp[23] = -1.0*GraP*Xul5P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 114:
            dwdp[23] = -1.0*Sed7P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 115:
            dwdp[23] = -1.0*GraP*Sed7P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 116:
            dwdp[23] = -1.0*GraP*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 117:
            dwdp[23] = -1.0*Rib5P*Sed7P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 118:
            dwdp[23] = -1.0*Rib5P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 119:
            dwdp[23] = -1.0*Xul5P*vTrKet1_Vmaxv23*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/pow(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23), 2);
            break;
        case 120:
            dwdp[23] = 1.0*GraP*Sed7P*vTrKet1_Vmaxv23/(pow(vTrKet1_Keqv23, 2)*(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23)));
            break;
        case 121:
            dwdp[23] = 1.0*(-GraP*Sed7P/vTrKet1_Keqv23 + Rib5P*Xul5P)/(GraP*Xul5P*vTrKet1_K7v23 + GraP*(Sed7P*vTrKet1_K5v23 + vTrKet1_K3v23) + Rib5P*(Sed7P*vTrKet1_K6v23 + vTrKet1_K2v23) + Sed7P*vTrKet1_K4v23 + Xul5P*(Rib5P + vTrKet1_K1v23));
            break;
        case 122:
            dwdp[24] = -1.0*E4P*Sed7P*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 123:
            dwdp[24] = -1.0*Fru6P*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 124:
            dwdp[24] = -1.0*E4P*Fru6P*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 125:
            dwdp[24] = -1.0*E4P*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 126:
            dwdp[24] = -1.0*Fru6P*GraP*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 127:
            dwdp[24] = -1.0*GraP*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 128:
            dwdp[24] = -1.0*Sed7P*vTrAld_Vmaxv24*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/pow(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24), 2);
            break;
        case 129:
            dwdp[24] = 1.0*E4P*Fru6P*vTrAld_Vmaxv24/(pow(vTrAld_Keqv24, 2)*(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24)));
            break;
        case 130:
            dwdp[24] = 1.0*(-E4P*Fru6P/vTrAld_Keqv24 + GraP*Sed7P)/(E4P*Sed7P*vTrAld_K7v24 + E4P*(Fru6P*vTrAld_K5v24 + vTrAld_K3v24) + Fru6P*vTrAld_K4v24 + GraP*(Fru6P*vTrAld_K6v24 + vTrAld_K2v24) + Sed7P*(GraP + vTrAld_K1v24));
            break;
        case 131:
            dwdp[25] = -1.0*vPPRPPS_Vmaxv25*(-MgAMP*PRPP/vPPRPPS_Keqv25 + MgATP*Rib5P)/((MgATP + vPPRPPS_KATPv25)*pow(Rib5P + vPPRPPS_KR5Pv25, 2));
            break;
        case 132:
            dwdp[25] = -1.0*vPPRPPS_Vmaxv25*(-MgAMP*PRPP/vPPRPPS_Keqv25 + MgATP*Rib5P)/(pow(MgATP + vPPRPPS_KATPv25, 2)*(Rib5P + vPPRPPS_KR5Pv25));
            break;
        case 133:
            dwdp[25] = 1.0*MgAMP*PRPP*vPPRPPS_Vmaxv25/(pow(vPPRPPS_Keqv25, 2)*(MgATP + vPPRPPS_KATPv25)*(Rib5P + vPPRPPS_KR5Pv25));
            break;
        case 134:
            dwdp[25] = 1.0*(-MgAMP*PRPP/vPPRPPS_Keqv25 + MgATP*Rib5P)/((MgATP + vPPRPPS_KATPv25)*(Rib5P + vPPRPPS_KR5Pv25));
            break;
        case 135:
            dwdp[26] = -1.0*GraP*Xul5P*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 136:
            dwdp[26] = -1.0*Fru6P*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 137:
            dwdp[26] = -1.0*Fru6P*GraP*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 138:
            dwdp[26] = -1.0*GraP*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 139:
            dwdp[26] = -1.0*E4P*Fru6P*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 140:
            dwdp[26] = -1.0*E4P*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 141:
            dwdp[26] = -1.0*Xul5P*vTrKet2_Vmaxv26*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/pow(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26), 2);
            break;
        case 142:
            dwdp[26] = 1.0*Fru6P*GraP*vTrKet2_Vmaxv26/(pow(vTrKet2_Keqv26, 2)*(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26)));
            break;
        case 143:
            dwdp[26] = 1.0*(E4P*Xul5P - Fru6P*GraP/vTrKet2_Keqv26)/(E4P*(Fru6P*vTrKet2_K6v26 + vTrKet2_K2v26) + Fru6P*vTrKet2_K4v26 + GraP*Xul5P*vTrKet2_K7v26 + GraP*(Fru6P*vTrKet2_K5v26 + vTrKet2_K3v26) + Xul5P*(E4P + vTrKet2_K1v26));
            break;
        case 144:
            dwdp[27] = 1.0*Phi*vPhiexch_Vmaxv27/pow(vPhiexch_Keqv27, 2);
            break;
        case 145:
            dwdp[27] = -1.0*Phi/vPhiexch_Keqv27 + 1.0*Phiex;
            break;
        case 146:
            dwdp[28] = 1.0*Lac*vLacexch_Vmaxv28/pow(vLacexch_Keqv28, 2);
            break;
        case 147:
            dwdp[28] = -1.0*Lac/vLacexch_Keqv28 + 1.0*Lacex;
            break;
        case 148:
            dwdp[29] = 1.0*Pyr*vPyrexch_Vmaxv29/pow(vPyrexch_Keqv29, 2);
            break;
        case 149:
            dwdp[29] = -1.0*Pyr/vPyrexch_Keqv29 + 1.0*Pyrex;
            break;
        case 150:
            dwdp[30] = 1.0*ATPf*Mgf*vMgATP_EqMult/pow(vMgATP_KdATP, 2);
            break;
        case 151:
            dwdp[30] = -1.0*ATPf*Mgf/vMgATP_KdATP + 1.0*MgATP;
            break;
        case 152:
            dwdp[31] = 1.0*ADPf*Mgf*vMgADP_EqMult/pow(vMgADP_KdADP, 2);
            break;
        case 153:
            dwdp[31] = -1.0*ADPf*Mgf/vMgADP_KdADP + 1.0*MgADP;
            break;
        case 154:
            dwdp[32] = 1.0*AMPf*Mgf*vMgAMP_EqMult/pow(vMgAMP_KdAMP, 2);
            break;
        case 155:
            dwdp[32] = -1.0*AMPf*Mgf/vMgAMP_KdAMP + 1.0*MgAMP;
            break;
        case 156:
            dwdp[33] = 1.0*Gri23P2f*Mgf*vMgGri23P2_EqMult/pow(vMgGri23P2_Kd23P2G, 2);
            break;
        case 157:
            dwdp[33] = -1.0*Gri23P2f*Mgf/vMgGri23P2_Kd23P2G + 1.0*MgGri23P2;
            break;
        case 158:
            dwdp[34] = 1.0*NADPf*P1f*vP1NADP_EqMult/pow(vP1NADP_Kd1, 2);
            break;
        case 159:
            dwdp[34] = -1.0*NADPf*P1f/vP1NADP_Kd1 + 1.0*P1NADP;
            break;
        case 160:
            dwdp[35] = 1.0*NADPHf*P1f*vP1NADPH_EqMult/pow(vP1NADPH_Kd3, 2);
            break;
        case 161:
            dwdp[35] = -1.0*NADPHf*P1f/vP1NADPH_Kd3 + 1.0*P1NADPH;
            break;
        case 162:
            dwdp[36] = 1.0*NADPf*P2f*vP2NADP_EqMult/pow(vP2NADP_Kd2, 2);
            break;
        case 163:
            dwdp[36] = -1.0*NADPf*P2f/vP2NADP_Kd2 + 1.0*P2NADP;
            break;
        case 164:
            dwdp[37] = 1.0*NADPHf*P2f*vP2NADPH_EqMult/pow(vP2NADPH_Kd4, 2);
            break;
        case 165:
            dwdp[37] = -1.0*NADPHf*P2f/vP2NADPH_Kd4 + 1.0*P2NADPH;
            break;
    }
}