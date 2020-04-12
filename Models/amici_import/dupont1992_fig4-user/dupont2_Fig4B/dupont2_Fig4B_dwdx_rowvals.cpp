#include "sundials/sundials_types.h"

void dwdx_rowvals_dupont2_Fig4B(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 0;
    rowvals[2] = 4;
    rowvals[3] = 6;
    rowvals[4] = 0;
    rowvals[5] = 3;
    rowvals[6] = 4;
    rowvals[7] = 5;
}