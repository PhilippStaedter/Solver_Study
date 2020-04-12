#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_li1(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 4;
    rowvals[2] = 5;
    rowvals[3] = 6;
    rowvals[4] = 7;
    rowvals[5] = 10;
    rowvals[6] = 1;
    rowvals[7] = 5;
    rowvals[8] = 8;
    rowvals[9] = 9;
    rowvals[10] = 10;
    rowvals[11] = 2;
    rowvals[12] = 6;
    rowvals[13] = 7;
}