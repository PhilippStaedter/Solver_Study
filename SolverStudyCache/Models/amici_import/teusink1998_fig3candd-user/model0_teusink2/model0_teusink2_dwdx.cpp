#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_teusink2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -ATP*Glc*VHK/(pow(KATP, 2)*KGlc*pow(ATP/KATP + 1, 2)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1)) + Glc*VHK/(KATP*KGlc*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
    dwdx[1] = ATP*HMP*VPFK*gR*(HMP*gR/(KRATP*KRHMP) + 1.0/KRATP)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2))) + ATP*HMP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)*(-L0*pow(ATP*ci/KiATP + 1, 2)*(2*HMP*c1*c2*gT/(KRATP*KRHMP) + 2*c2/KRATP)*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/pow(ATP/KiATP + 1, 2) - (2*HMP*gR/(KRATP*KRHMP) + 2/KRATP)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1) - 2*L0*ci*(ATP*ci/KiATP + 1)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(KiATP*pow(ATP/KiATP + 1, 2)) + 2*L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/(KiATP*pow(ATP/KiATP + 1, 3)))/(KRATP*KRHMP*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2)) + HMP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
    dwdx[2] = -Fru16P2*Vlower/(KADP*KFru16P2*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1)) + Fru16P2*Vlower*(-ATP + 5)/(pow(KADP, 2)*KFru16P2*pow(1 + (-ATP + 5)/KADP, 2)*(Fru16P2/KFru16P2 + 1));
    dwdx[3] = -ATP*VATPase/pow(ATP + KATP3, 2) + VATPase/(ATP + KATP3);
    dwdx[4] = -Fru16P2*Vlower*(-ATP + 5)/(KADP*pow(KFru16P2, 2)*(1 + (-ATP + 5)/KADP)*pow(Fru16P2/KFru16P2 + 1, 2)) + Vlower*(-ATP + 5)/(KADP*KFru16P2*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1));
    dwdx[5] = -ATP*Glc*VHK/(KATP*pow(KGlc, 2)*(ATP/KATP + 1)*pow(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1, 2)) + ATP*VHK/(KATP*KGlc*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
    dwdx[6] = -2*ATP*Glc*HMP*VHK*wildtype/(KATP*KGlc*KiTre6P*(ATP/KATP + 1)*pow(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1, 2));
    dwdx[7] = ATP*HMP*VPFK*gR*(ATP*gR/(KRATP*KRHMP) + 1.0/KRHMP)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2))) + ATP*HMP*VPFK*gR*(-L0*pow(ATP*ci/KiATP + 1, 2)*(2*ATP*c1*c2*gT/(KRATP*KRHMP) + 2*c1/KRHMP)*(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1)/pow(ATP/KiATP + 1, 2) - (2*ATP*gR/(KRATP*KRHMP) + 2/KRHMP)*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1))*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*pow(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2), 2)) + ATP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
}