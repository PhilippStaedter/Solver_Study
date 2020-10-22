#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_montanez1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 3;
    colptrs[3] = 4;
    colptrs[4] = 5;
    colptrs[5] = 6;
}