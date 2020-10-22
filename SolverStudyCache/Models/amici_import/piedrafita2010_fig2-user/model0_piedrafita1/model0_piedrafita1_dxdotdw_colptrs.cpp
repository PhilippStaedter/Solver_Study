#include "sundials/sundials_types.h"

void dxdotdw_colptrs_model0_piedrafita1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 2;
    colptrs[2] = 5;
    colptrs[3] = 6;
    colptrs[4] = 8;
    colptrs[5] = 11;
    colptrs[6] = 12;
    colptrs[7] = 15;
    colptrs[8] = 17;
    colptrs[9] = 20;
    colptrs[10] = 21;
    colptrs[11] = 23;
}