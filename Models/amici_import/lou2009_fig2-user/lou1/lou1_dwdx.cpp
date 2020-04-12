#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_lou1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = dI;
    dwdx[1] = St*betaBT/(Ib + Iv + Sb + Sv) - St*(Ib*betaBT + Iv*betaVT)/pow(Ib + Iv + Sb + Sv, 2);
    dwdx[2] = Sv*betaBV/(Ib + It + Iv + Sb + St + Sv) - Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2);
    dwdx[3] = Sb*betaTB/(It + Iv + St + Sv) - Sb*(It*betaTB + Iv*betaVB)/pow(It + Iv + St + Sv, 2);
    dwdx[4] = dI;
    dwdx[5] = Sv*betaTV/(Ib + It + Iv + Sb + St + Sv) - Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2);
    dwdx[6] = Sb*betaVB/(It + Iv + St + Sv) - Sb*(It*betaTB + Iv*betaVB)/pow(It + Iv + St + Sv, 2);
    dwdx[7] = St*betaVT/(Ib + Iv + Sb + Sv) - St*(Ib*betaBT + Iv*betaVT)/pow(Ib + Iv + Sb + Sv, 2);
    dwdx[8] = Sv*betaVV/(Ib + It + Iv + Sb + St + Sv) - Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2);
    dwdx[9] = dI;
    dwdx[10] = (It*betaTB + Iv*betaVB)/(It + Iv + St + Sv);
    dwdx[11] = dM;
    dwdx[12] = -St*(Ib*betaBT + Iv*betaVT)/pow(Ib + Iv + Sb + Sv, 2);
    dwdx[13] = -Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2);
    dwdx[14] = -Sb*(It*betaTB + Iv*betaVB)/pow(It + Iv + St + Sv, 2);
    dwdx[15] = (Ib*betaBT + Iv*betaVT)/(Ib + Iv + Sb + Sv);
    dwdx[16] = dM;
    dwdx[17] = -Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2);
    dwdx[18] = -Sb*(It*betaTB + Iv*betaVB)/pow(It + Iv + St + Sv, 2);
    dwdx[19] = -St*(Ib*betaBT + Iv*betaVT)/pow(Ib + Iv + Sb + Sv, 2);
    dwdx[20] = -Sv*(Ib*betaBV + It*betaTV + Iv*betaVV)/pow(Ib + It + Iv + Sb + St + Sv, 2) + (Ib*betaBV + It*betaTV + Iv*betaVV)/(Ib + It + Iv + Sb + St + Sv);
    dwdx[21] = dM;
}