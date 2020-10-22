#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_goldbeter6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
    colptrs[5] = 5;
    colptrs[6] = 6;
}