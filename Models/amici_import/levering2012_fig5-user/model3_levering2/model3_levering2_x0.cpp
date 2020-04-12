#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model3_levering2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.071888170346324706;
    x0[1] = species10init;
    x0[2] = 9.9999598735169304;
    x0[4] = 1.0000204776599499;
    x0[8] = 0.16;
    x0[9] = 4.5514359175799699;
    x0[11] = 0.246;
    x0[12] = 8.1919627519824392;
    x0[13] = 2.3399999999999999;
    x0[14] = 50.0;
    x0[17] = 1.0;
    x0[18] = 2.6549280134581998;
    x0[20] = 3.1062223058351002;
    x0[21] = 8.9296639189104301;
    x0[23] = 4.1379171389201597;
    x0[24] = parameter_1 - species10init;
}