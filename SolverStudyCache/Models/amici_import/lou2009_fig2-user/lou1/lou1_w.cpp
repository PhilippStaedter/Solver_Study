#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_lou1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = rt;
    w[1] = Sb*(It*betaTB + Iv*betaVB)/(It + Iv + St + Sv);
    w[2] = Sb*dM;
    w[3] = Ib*dI;
    w[4] = St*(Ib*betaBT + Iv*betaVT)/(Ib + Iv + Sb + Sv);
    w[5] = St*dM;
    w[6] = It*dI;
    w[7] = rv;
    w[8] = Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/(Ib + It + Iv + Sb + St + Sv);
    w[9] = Sv*dM;
    w[10] = Iv*dI;
    w[11] = rb;
}