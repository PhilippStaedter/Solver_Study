#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model1_perelson2(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 4;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 4;
    rowvals[7] = 2;
    rowvals[8] = 3;
    rowvals[9] = 0;
    rowvals[10] = 3;
    rowvals[11] = 1;
    rowvals[12] = 3;
    rowvals[13] = 0;
    rowvals[14] = 3;
    rowvals[15] = 0;
    rowvals[16] = 3;
}