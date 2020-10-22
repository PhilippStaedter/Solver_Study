#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kowald1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 1;
    rowvals[3] = 0;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 2;
    rowvals[7] = 1;
    rowvals[8] = 1;
    rowvals[9] = 3;
    rowvals[10] = 6;
    rowvals[11] = 3;
    rowvals[12] = 5;
    rowvals[13] = 6;
    rowvals[14] = 3;
    rowvals[15] = 0;
    rowvals[16] = 3;
    rowvals[17] = 5;
    rowvals[18] = 0;
    rowvals[19] = 2;
    rowvals[20] = 4;
    rowvals[21] = 2;
    rowvals[22] = 4;
    rowvals[23] = 2;
    rowvals[24] = 4;
    rowvals[25] = 2;
    rowvals[26] = 6;
    rowvals[27] = 4;
    rowvals[28] = 6;
    rowvals[29] = 5;
}