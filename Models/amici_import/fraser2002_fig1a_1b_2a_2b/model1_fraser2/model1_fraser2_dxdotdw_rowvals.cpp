#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_fraser2(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 5;
    rowvals[3] = 5;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 5;
    rowvals[7] = 0;
    rowvals[8] = 3;
    rowvals[9] = 3;
    rowvals[10] = 1;
    rowvals[11] = 1;
    rowvals[12] = 4;
    rowvals[13] = 4;
    rowvals[14] = 2;
}