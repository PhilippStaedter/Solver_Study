#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_Nakakuki2010(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = 10.0;
    x0[2] = 182.35239999999999;
    x0[4] = 25.38702;
    x0[6] = 13.09262;
    x0[9] = 570.41790000000003;
    x0[17] = 82.66574;
    x0[20] = 637.32119999999998;
    x0[23] = 353.0;
    x0[25] = 247.40350000000001;
    x0[27] = 1000.0;
    x0[29] = 1624.9000000000001;
    x0[32] = 1510.0;
}