#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_ferreira1_Fig2F(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*v2b_k2b*v2b_p2;
    dwdx[1] = 1.0*v4_k4*v4_ox*v4_p4;
    dwdx[2] = 1.0*Lysine*v1a_k1a*v1a_p1;
    dwdx[3] = 0.59298553243265262*pow(Glucose, -0.64000000000000001)*v3_k3*v3_ox*v3_p3;
    dwdx[4] = 1.0*Lysine*v5_k5*v5_ox*v5_p5;
    dwdx[5] = 1.0*v5b_k5b;
    dwdx[6] = 1.0*Glucose*v1a_k1a*v1a_p1;
    dwdx[7] = 1.0*Glyoxal*v5_k5*v5_ox*v5_p5;
    dwdx[8] = 1.0*v1b_k1b;
    dwdx[9] = 1.0*v2a_k2a*v2a_p2;
    dwdx[10] = 0.59298553243265262*pow(Schiff, -0.64000000000000001)*v6_k3*v6_ox*v6_p6;
    dwdx[11] = 0.02964927662163263*pow(Schiff, -0.64000000000000001)*v7a_k3*v7a_ox*v7a_p7;
    dwdx[12] = 0.002964927662163263*pow(Schiff, -0.64000000000000001)*v7b_k3*v7b_ox*v7b_p7;
    dwdx[13] = 0.59298553243265262*pow(Schiff, -0.64000000000000001)*v7c_k3*v7c_ox*v7c_p7;
}