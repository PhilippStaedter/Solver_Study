#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_teusink2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 1:
            dwdp[2] = -Fru16P2*Vlower*(-ATP + 5)/(pow(KADP, 2)*KFru16P2*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1)) + Fru16P2*Vlower*pow(-ATP + 5, 2)/(pow(KADP, 3)*KFru16P2*pow(1 + (-ATP + 5)/KADP, 2)*(Fru16P2/KFru16P2 + 1));
            break;
        case 2:
            dwdp[0] = pow(ATP, 2)*Glc*VHK/(pow(KATP, 3)*KGlc*pow(ATP/KATP + 1, 2)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1)) - ATP*Glc*VHK/(pow(KATP, 2)*KGlc*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
            break;
        case 3:
            dwdp[3] = -ATP*VATPase/pow(ATP + KATP3, 2);
            break;
        case 4:
            dwdp[2] = pow(Fru16P2, 2)*Vlower*(-ATP + 5)/(KADP*pow(KFru16P2, 3)*(1 + (-ATP + 5)/KADP)*pow(Fru16P2/KFru16P2 + 1, 2)) - Fru16P2*Vlower*(-ATP + 5)/(KADP*pow(KFru16P2, 2)*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1));
            break;
        case 5:
            dwdp[0] = ATP*pow(Glc, 2)*VHK/(KATP*pow(KGlc, 3)*(ATP/KATP + 1)*pow(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1, 2)) - ATP*Glc*VHK/(KATP*pow(KGlc, 2)*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
            break;
        case 6:
            dwdp[1] = ATP*HMP*VPFK*gR*(-ATP*HMP*gR/(pow(KRATP, 2)*KRHMP) - ATP/pow(KRATP, 2))/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2))) + ATP*HMP*VPFK*gR*(-L0*pow(ATP*ci/KiATP + 1, 2)*(-2*ATP*HMP*c1*c2*gT/(pow(KRATP, 2)*KRHMP) - 2*ATP*c2/pow(KRATP, 2))*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/pow(ATP/KiATP + 1, 2) - (-2*ATP*HMP*gR/(pow(KRATP, 2)*KRHMP) - 2*ATP/pow(KRATP, 2))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2)) - ATP*HMP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(pow(KRATP, 2)*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
            break;
        case 7:
            dwdp[1] = ATP*HMP*VPFK*gR*(-ATP*HMP*gR/(KRATP*pow(KRHMP, 2)) - HMP/pow(KRHMP, 2))/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2))) + ATP*HMP*VPFK*gR*(-L0*pow(ATP*ci/KiATP + 1, 2)*(-2*ATP*HMP*c1*c2*gT/(KRATP*pow(KRHMP, 2)) - 2*HMP*c1/pow(KRHMP, 2))*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/pow(ATP/KiATP + 1, 2) - (-2*ATP*HMP*gR/(KRATP*pow(KRHMP, 2)) - 2*HMP/pow(KRHMP, 2))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2)) - ATP*HMP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*pow(KRHMP, 2)*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
            break;
        case 8:
            dwdp[1] = ATP*HMP*VPFK*gR*(2*ATP*L0*ci*(ATP*ci/KiATP + 1)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(pow(KiATP, 2)*pow(ATP/KiATP + 1, 2)) - 2*ATP*L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(pow(KiATP, 2)*pow(ATP/KiATP + 1, 3)))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 9:
            dwdp[0] = ATP*Glc*pow(HMP, 2)*VHK*wildtype/(KATP*KGlc*pow(KiTre6P, 2)*(ATP/KATP + 1)*pow(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1, 2));
            break;
        case 10:
            dwdp[1] = -ATP*HMP*VPFK*gR*pow(ATP*ci/KiATP + 1, 2)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(KRATP*KRHMP*pow(ATP/KiATP + 1, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 11:
            dwdp[3] = ATP/(ATP + KATP3);
            break;
        case 12:
            dwdp[0] = ATP*Glc/(KATP*KGlc*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
            break;
        case 13:
            dwdp[1] = ATP*HMP*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
            break;
        case 14:
            dwdp[2] = Fru16P2*(-ATP + 5)/(KADP*KFru16P2*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1));
            break;
        case 15:
            dwdp[1] = -ATP*HMP*L0*VPFK*gR*pow(ATP*ci/KiATP + 1, 2)*(2*ATP*HMP*c2*gT/(KRATP*KRHMP) + 2*HMP/KRHMP)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/(KRATP*KRHMP*pow(ATP/KiATP + 1, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 16:
            dwdp[1] = -ATP*HMP*L0*VPFK*gR*pow(ATP*ci/KiATP + 1, 2)*(2*ATP*HMP*c1*gT/(KRATP*KRHMP) + 2*ATP/KRATP)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/(KRATP*KRHMP*pow(ATP/KiATP + 1, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 17:
            dwdp[1] = -2*pow(ATP, 2)*HMP*L0*VPFK*gR*(ATP*ci/KiATP + 1)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(KRATP*KRHMP*KiATP*pow(ATP/KiATP + 1, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 18:
            dwdp[1] = pow(ATP, 2)*pow(HMP, 2)*VPFK*gR/(pow(KRATP, 2)*pow(KRHMP, 2)*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2))) - 2*pow(ATP, 2)*pow(HMP, 2)*VPFK*gR*pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)/(pow(KRATP, 2)*pow(KRHMP, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2)) + ATP*HMP*VPFK*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
            break;
        case 19:
            dwdp[1] = -2*pow(ATP, 2)*pow(HMP, 2)*L0*VPFK*c1*c2*gR*pow(ATP*ci/KiATP + 1, 2)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/(pow(KRATP, 2)*pow(KRHMP, 2)*pow(ATP/KiATP + 1, 2)*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2));
            break;
        case 20:
            dwdp[0] = -ATP*Glc*pow(HMP, 2)*VHK/(KATP*KGlc*KiTre6P*(ATP/KATP + 1)*pow(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1, 2));
            break;
    }
}