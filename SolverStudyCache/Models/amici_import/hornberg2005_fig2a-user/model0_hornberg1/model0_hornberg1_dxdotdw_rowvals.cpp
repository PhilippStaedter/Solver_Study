#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_hornberg1(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 2;
    rowvals[2] = 1;
    rowvals[3] = 2;
    rowvals[4] = 3;
    rowvals[5] = 4;
    rowvals[6] = 3;
    rowvals[7] = 4;
    rowvals[8] = 5;
    rowvals[9] = 6;
    rowvals[10] = 5;
    rowvals[11] = 6;
    rowvals[12] = 7;
    rowvals[13] = 8;
    rowvals[14] = 7;
    rowvals[15] = 8;
}