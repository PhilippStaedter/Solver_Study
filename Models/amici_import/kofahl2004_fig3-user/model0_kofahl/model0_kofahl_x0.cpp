#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model0_kofahl(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 200.0;
    x0[3] = 300.0;
    x0[4] = 500.0;
    x0[7] = 686.39970164051294;
    x0[11] = 1666.6666667;
    x0[14] = 158.33176608789;
    x0[15] = 200.0;
    x0[17] = 1666.6666667;
    x0[18] = 1000.0;
    x0[20] = 158.33176608789;
    x0[21] = 36.399701640514103;
    x0[22] = 100.0;
    x0[23] = 105.943298120207;
    x0[24] = 77.875362567582897;
    x0[25] = 235.72493579190299;
}