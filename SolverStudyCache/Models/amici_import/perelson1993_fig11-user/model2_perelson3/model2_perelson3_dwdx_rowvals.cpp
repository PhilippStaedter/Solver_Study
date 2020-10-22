#include "sundials/sundials_types.h"

void dwdx_rowvals_model2_perelson3(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 5;
    rowvals[2] = 7;
    rowvals[3] = 2;
    rowvals[4] = 6;
    rowvals[5] = 0;
    rowvals[6] = 4;
    rowvals[7] = 7;
    rowvals[8] = 1;
    rowvals[9] = 3;
    rowvals[10] = 8;
}