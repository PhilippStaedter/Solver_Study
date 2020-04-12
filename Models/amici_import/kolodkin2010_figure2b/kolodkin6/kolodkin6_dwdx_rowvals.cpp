#include "sundials/sundials_types.h"

void dwdx_rowvals_kolodkin6(sunindextype *rowvals){
    rowvals[0] = 1;
    rowvals[1] = 3;
    rowvals[2] = 1;
    rowvals[3] = 2;
    rowvals[4] = 2;
    rowvals[5] = 4;
    rowvals[6] = 3;
    rowvals[7] = 5;
    rowvals[8] = 0;
    rowvals[9] = 2;
    rowvals[10] = 5;
    rowvals[11] = 3;
    rowvals[12] = 4;
    rowvals[13] = 0;
    rowvals[14] = 0;
}