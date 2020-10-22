#include "sundials/sundials_types.h"

void dxdotdw_rowvals_kolodkin1(sunindextype *rowvals){
    rowvals[0] = 0;
    rowvals[1] = 1;
    rowvals[2] = 2;
    rowvals[3] = 0;
}