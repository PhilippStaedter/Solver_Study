#include "sundials/sundials_types.h"

void JSparseB_colptrs_kolodkin6(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 4;
    colptrs[3] = 8;
    colptrs[4] = 12;
    colptrs[5] = 18;
    colptrs[6] = 22;
    colptrs[7] = 25;
    colptrs[8] = 28;
}