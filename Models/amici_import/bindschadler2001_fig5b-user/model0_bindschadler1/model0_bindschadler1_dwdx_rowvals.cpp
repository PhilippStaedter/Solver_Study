#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_bindschadler1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 4;
    rowvals[2] = 6;
    rowvals[3] = 8;
    rowvals[4] = 10;
    rowvals[5] = 1;
    rowvals[6] = 5;
    rowvals[7] = 7;
    rowvals[8] = 9;
    rowvals[9] = 10;
    rowvals[10] = 0;
    rowvals[11] = 6;
    rowvals[12] = 8;
    rowvals[13] = 1;
    rowvals[14] = 7;
    rowvals[15] = 9;
}