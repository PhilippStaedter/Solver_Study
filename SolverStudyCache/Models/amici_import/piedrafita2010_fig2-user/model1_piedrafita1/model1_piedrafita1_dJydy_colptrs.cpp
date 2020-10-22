#include "sundials/sundials_types.h"

void dJydy_colptrs_model1_piedrafita1(sunindextype *colptrs, int index){
    switch(index) {
        case 0:
                colptrs[0] = 0;
                colptrs[1] = 1;
                colptrs[2] = 1;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 1:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 1;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 2:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 1;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 3:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 1;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 4:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 1;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 5:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 1;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 6:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 1;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 7:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 0;
                colptrs[8] = 1;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 8:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 0;
                colptrs[8] = 0;
                colptrs[9] = 1;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 9:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 0;
                colptrs[8] = 0;
                colptrs[9] = 0;
                colptrs[10] = 1;
                colptrs[11] = 1;
            break;
        case 10:
                colptrs[0] = 0;
                colptrs[1] = 0;
                colptrs[2] = 0;
                colptrs[3] = 0;
                colptrs[4] = 0;
                colptrs[5] = 0;
                colptrs[6] = 0;
                colptrs[7] = 0;
                colptrs[8] = 0;
                colptrs[9] = 0;
                colptrs[10] = 0;
                colptrs[11] = 1;
            break;
    }
}