#include "sundials/sundials_types.h"

void JSparseB_colptrs_model3_mellor1(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 1;
    colptrs[2] = 2;
    colptrs[3] = 3;
    colptrs[4] = 4;
    colptrs[5] = 5;
    colptrs[6] = 6;
    colptrs[7] = 8;
    colptrs[8] = 10;
    colptrs[9] = 12;
    colptrs[10] = 13;
}