#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void dydx_Ouzounoglou2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[91] = 1;
    dydx[182] = 1;
    dydx[273] = 1;
    dydx[364] = 1;
    dydx[455] = 1;
    dydx[546] = 1;
    dydx[637] = 1;
    dydx[728] = 1;
    dydx[819] = 1;
    dydx[910] = 1;
    dydx[1001] = 1;
    dydx[1092] = 1;
    dydx[1183] = 1;
    dydx[1274] = 1;
    dydx[1365] = 1;
    dydx[1456] = 1;
    dydx[1547] = 1;
    dydx[1638] = 1;
    dydx[1729] = 1;
    dydx[1820] = 1;
    dydx[1911] = 1;
    dydx[2002] = 1;
    dydx[2093] = 1;
    dydx[2184] = 1;
    dydx[2275] = 1;
    dydx[2366] = 1;
    dydx[2457] = 1;
    dydx[2548] = 1;
    dydx[2639] = 1;
    dydx[2730] = 1;
    dydx[2821] = 1;
    dydx[2912] = 1;
    dydx[3003] = 1;
    dydx[3094] = 1;
    dydx[3185] = 1;
    dydx[3276] = 1;
    dydx[3367] = 1;
    dydx[3458] = 1;
    dydx[3549] = 1;
    dydx[3640] = 1;
    dydx[3731] = 1;
    dydx[3822] = 1;
    dydx[3913] = 1;
    dydx[4004] = 1;
    dydx[4095] = 1;
    dydx[4186] = 1;
    dydx[4277] = 1;
    dydx[4368] = 1;
    dydx[4459] = 1;
    dydx[4550] = 1;
    dydx[4641] = 1;
    dydx[4732] = 1;
    dydx[4823] = 1;
    dydx[4914] = 1;
    dydx[5005] = 1;
    dydx[5096] = 1;
    dydx[5187] = 1;
    dydx[5278] = 1;
    dydx[5369] = 1;
    dydx[5460] = 1;
    dydx[5551] = 1;
    dydx[5642] = 1;
    dydx[5733] = 1;
    dydx[5824] = 1;
    dydx[5915] = 1;
    dydx[6006] = 1;
    dydx[6097] = 1;
    dydx[6188] = 1;
    dydx[6279] = 1;
    dydx[6370] = 1;
    dydx[6461] = 1;
    dydx[6552] = 1;
    dydx[6643] = 1;
    dydx[6734] = 1;
    dydx[6825] = 1;
    dydx[6916] = 1;
    dydx[7007] = 1;
    dydx[7098] = 1;
    dydx[7189] = 1;
    dydx[7280] = 1;
    dydx[7371] = 1;
    dydx[7462] = 1;
    dydx[7553] = 1;
    dydx[7644] = 1;
    dydx[7735] = 1;
    dydx[7826] = 1;
    dydx[7917] = 1;
    dydx[8008] = 1;
    dydx[8099] = 1;
}