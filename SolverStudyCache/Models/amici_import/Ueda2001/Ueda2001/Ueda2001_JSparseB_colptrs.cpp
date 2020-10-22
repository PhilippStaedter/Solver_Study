#include "sundials/sundials_types.h"

void JSparseB_colptrs_Ueda2001(sunindextype *colptrs){
    colptrs[0] = 0;
    colptrs[1] = 0;
    colptrs[2] = 4;
    colptrs[3] = 6;
    colptrs[4] = 10;
    colptrs[5] = 13;
    colptrs[6] = 18;
    colptrs[7] = 21;
    colptrs[8] = 25;
    colptrs[9] = 27;
    colptrs[10] = 31;
    colptrs[11] = 34;
    colptrs[12] = 34;
    colptrs[13] = 34;
}