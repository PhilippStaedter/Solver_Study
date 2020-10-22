#include "sundials/sundials_types.h"

void JSparseB_colptrs_kolodkin3(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 6;
    colptrs[3] = 10;
    colptrs[4] = 14;
    colptrs[5] = 19;
    colptrs[6] = 23;
    colptrs[7] = 27;
    colptrs[8] = 30;
    colptrs[9] = 33;
    colptrs[10] = 36;
}