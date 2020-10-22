#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model2_(sunindextype *rowvals){
    rowvals[0] = 3;
    rowvals[1] = 4;
    rowvals[2] = 3;
    rowvals[3] = 4;
    rowvals[4] = 6;
    rowvals[5] = 3;
    rowvals[6] = 6;
    rowvals[7] = 6;
    rowvals[8] = 4;
    rowvals[9] = 6;
}