#include "sundials/sundials_types.h"

void dwdx_rowvals_model1_alexander1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 5;
    rowvals[3] = 6;
    rowvals[4] = 7;
    rowvals[5] = 8;
    rowvals[6] = 4;
    rowvals[7] = 6;
    rowvals[8] = 10;
    rowvals[9] = 2;
    rowvals[10] = 3;
    rowvals[11] = 11;
    rowvals[12] = 1;
    rowvals[13] = 9;
}