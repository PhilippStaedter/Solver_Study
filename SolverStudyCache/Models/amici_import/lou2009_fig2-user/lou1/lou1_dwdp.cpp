#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_lou1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[4] = Ib*St/(Ib + Iv + Sb + Sv);
            break;
        case 1:
            dwdp[8] = Ib*Sv/(Ib + It + Iv + Sb + St + Sv);
            break;
        case 2:
            dwdp[1] = It*Sb/(It + Iv + St + Sv);
            break;
        case 3:
            dwdp[8] = It*Sv/(Ib + It + Iv + Sb + St + Sv);
            break;
        case 4:
            dwdp[1] = Iv*Sb/(It + Iv + St + Sv);
            break;
        case 5:
            dwdp[4] = Iv*St/(Ib + Iv + Sb + Sv);
            break;
        case 6:
            dwdp[8] = Iv*Sv/(Ib + It + Iv + Sb + St + Sv);
            break;
        case 7:
            dwdp[3] = Ib;
            dwdp[6] = It;
            dwdp[10] = Iv;
            break;
        case 8:
            dwdp[2] = Sb;
            dwdp[5] = St;
            dwdp[9] = Sv;
            break;
        case 9:
            dwdp[11] = 1;
            break;
        case 10:
            dwdp[0] = 1;
            break;
        case 11:
            dwdp[7] = 1;
            break;
    }
}