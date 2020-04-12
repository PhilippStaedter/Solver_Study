#include "sundials/sundials_types.h"

void dwdx_rowvals_model0_hornberg1(sunindextype *rowvals){
    rowvals[0] = 7;
    rowvals[1] = 0;
    rowvals[2] = 2;
    rowvals[3] = 1;
    rowvals[4] = 2;
    rowvals[5] = 3;
    rowvals[6] = 4;
    rowvals[7] = 4;
    rowvals[8] = 5;
    rowvals[9] = 6;
    rowvals[10] = 6;
    rowvals[11] = 7;
}