#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_gorlich1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 4;
    rowvals[3] = 8;
    rowvals[4] = 1;
    rowvals[5] = 2;
    rowvals[6] = 1;
    rowvals[7] = 4;
    rowvals[8] = 2;
    rowvals[9] = 8;
    rowvals[10] = 7;
    rowvals[11] = 5;
    rowvals[12] = 6;
    rowvals[13] = 3;
    rowvals[14] = 3;
    rowvals[15] = 4;
    rowvals[16] = 5;
    rowvals[17] = 7;
    rowvals[18] = 0;
    rowvals[19] = 6;
    rowvals[20] = 7;
    rowvals[21] = 0;
    rowvals[22] = 8;
}