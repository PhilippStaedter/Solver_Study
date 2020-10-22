#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kholodenko2(sunindextype *rowvals){
    rowvals[0] = 4;
    rowvals[1] = 5;
    rowvals[2] = 4;
    rowvals[3] = 5;
    rowvals[4] = 3;
    rowvals[5] = 6;
    rowvals[6] = 6;
    rowvals[7] = 7;
    rowvals[8] = 6;
    rowvals[9] = 7;
    rowvals[10] = 3;
    rowvals[11] = 6;
    rowvals[12] = 0;
    rowvals[13] = 1;
    rowvals[14] = 1;
    rowvals[15] = 2;
    rowvals[16] = 1;
    rowvals[17] = 2;
    rowvals[18] = 0;
    rowvals[19] = 1;
}