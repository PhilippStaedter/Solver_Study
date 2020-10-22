#include "sundials/sundials_types.h"

void dwdx_rowvals_kolodkin3(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 6;
    rowvals[2] = 1;
    rowvals[3] = 2;
    rowvals[4] = 3;
    rowvals[5] = 5;
    rowvals[6] = 6;
    rowvals[7] = 3;
    rowvals[8] = 5;
    rowvals[9] = 0;
    rowvals[10] = 2;
    rowvals[11] = 4;
    rowvals[12] = 6;
    rowvals[13] = 3;
    rowvals[14] = 4;
    rowvals[15] = 2;
    rowvals[16] = 0;
    rowvals[17] = 0;
}