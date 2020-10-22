#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_lee2(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 3;
    colptrs[2] = 6;
    colptrs[3] = 9;
    colptrs[4] = 12;
    colptrs[5] = 15;
    colptrs[6] = 18;
    colptrs[7] = 21;
    colptrs[8] = 24;
}