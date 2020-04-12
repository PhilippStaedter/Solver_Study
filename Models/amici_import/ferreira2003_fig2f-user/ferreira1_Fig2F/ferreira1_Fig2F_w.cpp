#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_ferreira1_Fig2F(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Glucose*Lysine*v1a_k1a*v1a_p1;
    w[1] = 1.0*Schiff*v1b_k1b;
    w[2] = 1.0*Schiff*v2a_k2a*v2a_p2;
    w[3] = 1.0*Amadori*v2b_k2b*v2b_p2;
    w[4] = 1.6471820345351462*pow(Glucose, 0.35999999999999999)*v3_k3*v3_ox*v3_p3;
    w[5] = 1.0*Amadori*v4_k4*v4_ox*v4_p4;
    w[6] = 1.0*Glyoxal*Lysine*v5_k5*v5_ox*v5_p5;
    w[7] = 1.0*Glyoxal*v5b_k5b;
    w[8] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v6_k3*v6_ox*v6_p6;
    w[9] = 0.082359101726757311*pow(Schiff, 0.35999999999999999)*v7a_k3*v7a_ox*v7a_p7;
    w[10] = 0.0082359101726757304*pow(Schiff, 0.35999999999999999)*v7b_k3*v7b_ox*v7b_p7;
    w[11] = 1.6471820345351462*pow(Schiff, 0.35999999999999999)*v7c_k3*v7c_ox*v7c_p7;
}