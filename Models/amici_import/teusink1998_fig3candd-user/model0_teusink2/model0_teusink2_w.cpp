#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_teusink2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = ATP*Glc*VHK/(KATP*KGlc*(ATP/KATP + 1)*(Glc/KGlc + pow(HMP, 2)*wildtype/KiTre6P + 1));
    w[1] = ATP*HMP*VPFK*gR*(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1)/(KRATP*KRHMP*(L0*pow(ATP*ci/KiATP + 1, 2)*pow(ATP*HMP*c1*c2*gT/(KRATP*KRHMP) + ATP*c2/KRATP + HMP*c1/KRHMP + 1, 2)/pow(ATP/KiATP + 1, 2) + pow(ATP*HMP*gR/(KRATP*KRHMP) + ATP/KRATP + HMP/KRHMP + 1, 2)));
    w[2] = Fru16P2*Vlower*(-ATP + 5)/(KADP*KFru16P2*(1 + (-ATP + 5)/KADP)*(Fru16P2/KFru16P2 + 1));
    w[3] = ATP*VATPase/(ATP + KATP3);
}