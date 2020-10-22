#include "sundials/sundials_types.h"

void dxdotdw_rowvals_model2_balagadde1(sunindextype *rowvals){
    rowvals[0] = 2;
    rowvals[1] = 2;
    rowvals[2] = 3;
    rowvals[3] = 3;
    rowvals[4] = 0;
    rowvals[5] = 0;
    rowvals[6] = 1;
    rowvals[7] = 1;
}