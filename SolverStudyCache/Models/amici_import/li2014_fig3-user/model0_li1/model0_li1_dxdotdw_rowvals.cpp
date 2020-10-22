#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_li1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 2;
    rowvals[2] = 2;
    rowvals[3] = 0;
    rowvals[4] = 0;
    rowvals[5] = 0;
    rowvals[6] = 0;
    rowvals[7] = 1;
    rowvals[8] = 1;
    rowvals[9] = 1;
    rowvals[10] = 1;
}