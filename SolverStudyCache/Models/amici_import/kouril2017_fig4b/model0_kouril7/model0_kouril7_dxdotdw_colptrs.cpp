#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_kouril7(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 4;
    colptrs[2] = 8;
    colptrs[3] = 12;
    colptrs[4] = 15;
}