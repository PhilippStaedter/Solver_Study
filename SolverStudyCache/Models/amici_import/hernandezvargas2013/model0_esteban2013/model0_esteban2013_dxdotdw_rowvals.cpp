#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model0_esteban2013(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 1;
    rowvals[2] = 4;
    rowvals[3] = 4;
    rowvals[4] = 4;
    rowvals[5] = 2;
    rowvals[6] = 2;
    rowvals[7] = 3;
    rowvals[8] = 2;
    rowvals[9] = 3;
    rowvals[10] = 0;
    rowvals[11] = 0;
    rowvals[12] = 0;
    rowvals[13] = 1;
    rowvals[14] = 0;
}