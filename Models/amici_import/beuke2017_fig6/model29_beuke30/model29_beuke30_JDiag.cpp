#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model29_beuke30(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -76923.076923076922*dwdx0 - 76923.076923076922*dwdx1 - 76923.076923076922*dwdx2 - 76923.076923076922*dwdx3;
    JDiag[1] = -76923.076923076922*dwdx4;
    JDiag[2] = -1250000.0*dwdx6 - 1250000.0*dwdx7;
    JDiag[3] = -1250000.0*dwdx8;
    JDiag[4] = -1250000.0*dwdx9;
    JDiag[5] = -76923.076923076922*dwdx10 - 76923.076923076922*dwdx11;
    JDiag[6] = -76923.076923076922*dwdx12;
    JDiag[7] = -386100.38610038609*dwdx13 - 386100.38610038609*dwdx14;
    JDiag[8] = -1250000.0*dwdx17;
    JDiag[9] = -1250000.0*dwdx18;
    JDiag[10] = 386100.38610038609*dwdx21 - 386100.38610038609*dwdx22;
    JDiag[11] = -76923.076923076922*dwdx23;
    JDiag[12] = 386100.38610038609*dwdx26 - 386100.38610038609*dwdx27;
    JDiag[13] = -76923.076923076922*dwdx29;
    JDiag[14] = -2127659.5744680851*dwdx31;
    JDiag[15] = -2857142.8571428573*dwdx33;
    JDiag[16] = -2127659.5744680851*dwdx34;
    JDiag[17] = -2857142.8571428573*dwdx35;
    JDiag[18] = -76923.076923076922*dwdx36 - 76923.076923076922*dwdx37;
    JDiag[19] = -76923.076923076922*dwdx39;
    JDiag[20] = -1250000.0*dwdx40 - 1250000.0*dwdx41;
    JDiag[21] = -1250000.0*dwdx43 - 1250000.0*dwdx44 - 1250000.0*dwdx45 + 1250000.0*dwdx46;
    JDiag[22] = -76923.076923076922*dwdx47 - 76923.076923076922*dwdx48;
    JDiag[24] = 1250000.0*dwdx50 - 1250000.0*dwdx51 - 1250000.0*dwdx53;
    JDiag[25] = -1250000.0*dwdx54 + 1250000.0*dwdx55 - 1250000.0*dwdx56;
    JDiag[26] = -76923.076923076922*dwdx57 - 76923.076923076922*dwdx58;
    JDiag[27] = -76923.076923076922*dwdx61;
    JDiag[28] = -386100.38610038609*dwdx62 - 386100.38610038609*dwdx63;
    JDiag[29] = 76923.076923076922*dwdx64 - 76923.076923076922*dwdx65 - 76923.076923076922*dwdx66 + 76923.076923076922*dwdx67 - 76923.076923076922*dwdx68;
    JDiag[30] = -76923.076923076922*dwdx69 - 76923.076923076922*dwdx70 - 76923.076923076922*dwdx71;
    JDiag[31] = -76923.076923076922*dwdx72;
    JDiag[32] = -1250000.0*dwdx73 - 1250000.0*dwdx74 - 1250000.0*dwdx75 + 1250000.0*dwdx77;
}