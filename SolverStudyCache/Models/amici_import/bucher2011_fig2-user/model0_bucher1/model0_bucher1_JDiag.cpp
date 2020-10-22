#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_bucher1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 70.422535211267601*dwdx0;
    JDiag[1] = -70.422535211267601*dwdx1 - 70.422535211267601*dwdx2 - 70.422535211267601*dwdx3 - 70.422535211267601*dwdx6 - 70.422535211267601*dwdx7;
    JDiag[2] = -0.5*dwdx8 - 0.5*dwdx9;
    JDiag[3] = 70.422535211267601*dwdx10;
    JDiag[4] = -70.422535211267601*dwdx11 - 70.422535211267601*dwdx12 - 70.422535211267601*dwdx13;
    JDiag[5] = -0.5*dwdx14 - 0.5*dwdx15;
    JDiag[6] = 70.422535211267601*dwdx16;
    JDiag[7] = -70.422535211267601*dwdx17 - 70.422535211267601*dwdx18 - 70.422535211267601*dwdx19;
    JDiag[8] = -0.5*dwdx20 - 0.5*dwdx21;
    JDiag[9] = 70.422535211267601*dwdx22;
    JDiag[10] = -70.422535211267601*dwdx23 - 70.422535211267601*dwdx26 - 70.422535211267601*dwdx27 - 70.422535211267601*dwdx28 - 70.422535211267601*dwdx29;
    JDiag[11] = -0.5*dwdx30;
    JDiag[12] = 70.422535211267601*dwdx31;
    JDiag[13] = -70.422535211267601*dwdx32 - 70.422535211267601*dwdx33;
    JDiag[14] = -0.5*dwdx34;
    JDiag[15] = 70.422535211267601*dwdx35;
    JDiag[16] = -70.422535211267601*dwdx36 - 70.422535211267601*dwdx37;
    JDiag[17] = -0.5*dwdx38;
}