#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "p.h"
#include "k.h"

void x0_model5_levchenko2(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = initial_value_C1_Fraction*initial_value_Total_Scaffold;
    x0[1] = initial_value_C2_Fraction*initial_value_Total_Scaffold;
    x0[3] = initial_value_C4_Fraction*initial_value_Total_Scaffold;
    x0[5] = initial_value_C6_Fraction*initial_value_Total_Scaffold;
    x0[9] = initial_value_MAPK;
    x0[11] = 0.29999999999999999;
    x0[17] = initial_value_MEK;
    x0[18] = 0.20000000000000001;
    x0[25] = 0.29999999999999999;
    x0[26] = initial_value_RAFK;
    x0[27] = 0.29999999999999999;
}