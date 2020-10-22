#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_kouril8(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 4;
    rowvals[4] = 0;
    rowvals[5] = 1;
    rowvals[6] = 9;
    rowvals[7] = 11;
    rowvals[8] = 2;
    rowvals[9] = 4;
    rowvals[10] = 2;
    rowvals[11] = 5;
    rowvals[12] = 7;
    rowvals[13] = 8;
    rowvals[14] = 6;
    rowvals[15] = 7;
    rowvals[16] = 8;
    rowvals[17] = 5;
}