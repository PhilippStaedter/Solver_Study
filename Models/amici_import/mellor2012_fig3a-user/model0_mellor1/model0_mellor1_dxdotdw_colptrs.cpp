#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_mellor1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 9;
    colptrs[2] = 18;
    colptrs[3] = 27;
    colptrs[4] = 29;
    colptrs[5] = 31;
}