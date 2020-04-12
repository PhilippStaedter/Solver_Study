#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_hald(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 3:
            dwdp[12] = -0.69999999999999996*Pyr*V1PDC*sc*(G6P/KxPDC + 1)/(pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2)*(G6P*bet/(KxPDC*alp) + 1));
            break;
        case 4:
            dwdp[1] = -2*NADH*pow(V1adh1, 2)*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)/(Keq*pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2));
            break;
        case 5:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-Kiadh1NAD*V2adh1 - NAD*V2adh1)*(-ACA*NADH/Keq + EtOH*NAD)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
            break;
        case 6:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-EtOH*V2adh1 - EtOH*NADH*V2adh1/Kiadh1NADH)*(-ACA*NADH/Keq + EtOH*NAD)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
            break;
        case 7:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*V1adh1/Keq - ACA*NAD*V1adh1/(Keq*Kiadh1NAD))/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
            break;
        case 8:
            dwdp[3] = sc*(DHAP*pow(GAP, 2)*pow(Valdf, 2)/(pow(Kaldeq, 2)*Valdr*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) - FBP*GAP*pow(Valdf, 2)/(Kaldeq*Valdr*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)));
            break;
        case 9:
            dwdp[3] = sc*(DHAP*GAP*Valdf/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) - FBP*Valdf/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2));
            break;
        case 10:
            dwdp[3] = sc*(pow(DHAP, 2)*GAP*pow(Valdf, 2)/(pow(Kaldeq, 2)*Valdr*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) - DHAP*FBP*pow(Valdf, 2)/(Kaldeq*Valdr*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)));
            break;
        case 11:
            dwdp[3] = sc*(-DHAP*FBP*pow(GAP, 2)*Valdf/(pow(KaldIGAP, 2)*Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) + pow(FBP, 2)*GAP*Valdf/(pow(KaldIGAP, 2)*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)));
            break;
        case 12:
            dwdp[3] = sc*(-DHAP*GAP*Valdf*(DHAP*GAP*Valdf/(pow(Kaldeq, 2)*Valdr) + DHAP*KaldGAP*Valdf/(pow(Kaldeq, 2)*Valdr) + GAP*KaldDHAP*Valdf/(pow(Kaldeq, 2)*Valdr))/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) + DHAP*GAP*Valdf/(pow(Kaldeq, 2)*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)) + FBP*Valdf*(DHAP*GAP*Valdf/(pow(Kaldeq, 2)*Valdr) + DHAP*KaldGAP*Valdf/(pow(Kaldeq, 2)*Valdr) + GAP*KaldDHAP*Valdf/(pow(Kaldeq, 2)*Valdr))/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2));
            break;
        case 13:
            dwdp[5] = -0.83199999999999996*ATP*Vatpconm*sc/pow(ATP + KatpconATP, 2);
            break;
        case 14:
            dwdp[1] = 2*ACA*NADH*V1adh1*V2adh1*sc/(pow(Keq, 2)*(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1)) + 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(ACA*EtOH*NADH*V1adh1/(pow(Keq, 2)*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/pow(Keq, 2) + ACA*Kadh1NADH*NAD*V1adh1/(pow(Keq, 2)*Kiadh1NAD) + ACA*NADH*V1adh1/pow(Keq, 2) + Kadh1ACA*NADH*V1adh1/pow(Keq, 2))/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
            break;
        case 15:
            dwdp[7] = sc*(-pow(DPG, 2)*NADH*Vgapdh/(pow(KgapdhDPG, 2)*KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) + DPG*GAP*NAD*Vgapdh/(pow(KgapdhDPG, 2)*KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)));
            break;
        case 16:
            dwdp[7] = sc*(-DPG*GAP*NADH*Vgapdh/(pow(KgapdhGAP, 3)*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) + DPG*NADH*Vgapdh/(pow(KgapdhGAP, 2)*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + pow(GAP, 2)*NAD*Vgapdh/(pow(KgapdhGAP, 3)*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*pow(DPG/KgapdhDPG + GAP/KgapdhGAP + 1, 2)) - GAP*NAD*Vgapdh/(pow(KgapdhGAP, 2)*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
            break;
        case 17:
            dwdp[7] = sc*(DPG*NADH*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 2)*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) - DPG*NAD*NADH*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 3)*Kgapdheq*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) - GAP*NAD*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 2)*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*pow(NAD, 2)*Vgapdh/(KgapdhGAP*pow(KgapdhNAD, 3)*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
            break;
        case 18:
            dwdp[7] = sc*(-DPG*pow(NADH, 2)*Vgapdh/(KgapdhGAP*KgapdhNAD*pow(KgapdhNADH, 2)*Kgapdheq*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*NAD*NADH*Vgapdh/(KgapdhGAP*KgapdhNAD*pow(KgapdhNADH, 2)*pow(1 + NADH/KgapdhNADH + NAD/KgapdhNAD, 2)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
            break;
        case 19:
            dwdp[7] = DPG*NADH*Vgapdh*sc/(KgapdhGAP*KgapdhNAD*pow(Kgapdheq, 2)*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1));
            break;
        case 20:
            dwdp[8] = sc*(-Glc*Vglctr*(G6P*Glc/(pow(KglctrGlc, 2)*KglctrIIG6P) + Glc*Pglctr*(Glc0/KglctrGlc + 1)/(pow(KglctrGlc, 2)*(Glc0*Pglctr/KglctrGlc + 1)) + Glc/pow(KglctrGlc, 2) - Glc0*Pglctr*(Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(pow(KglctrGlc, 2)*pow(Glc0*Pglctr/KglctrGlc + 1, 2)) + Glc0*(Glc*Pglctr/KglctrGlc + 1)/(pow(KglctrGlc, 2)*(Glc0*Pglctr/KglctrGlc + 1)))/(KglctrGlc*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + Glc*Vglctr/(pow(KglctrGlc, 2)*Vratio*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)) + Glc0*Vglctr*(-Glc*Pglctr*(Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(pow(KglctrGlc, 2)*pow(Glc*Pglctr/KglctrGlc + 1, 2)) + Glc0*Pglctr*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(pow(KglctrGlc, 2)*(Glc*Pglctr/KglctrGlc + 1)) + Glc0/pow(KglctrGlc, 2) - (Glc0*Pglctr/KglctrGlc + 1)*(-G6P*Glc/(pow(KglctrGlc, 2)*KglctrIIG6P) - Glc/pow(KglctrGlc, 2))/(Glc*Pglctr/KglctrGlc + 1))/(KglctrGlc*Vratio*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)) - Glc0*Vglctr/(pow(KglctrGlc, 2)*Vratio*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1))));
            break;
        case 21:
            dwdp[8] = sc*(-G6P*Glc*Vglctr/(KglctrGlc*pow(KglctrIG6P, 2)*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + G6P*Glc0*Vglctr*(Glc0*Pglctr/KglctrGlc + 1)/(KglctrGlc*pow(KglctrIG6P, 2)*Vratio*(Glc*Pglctr/KglctrGlc + 1)*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)));
            break;
        case 22:
            dwdp[8] = sc*(-G6P*pow(Glc, 2)*Vglctr/(pow(KglctrGlc, 2)*pow(KglctrIIG6P, 2)*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + G6P*Glc*Glc0*Vglctr*(Glc0*Pglctr/KglctrGlc + 1)/(pow(KglctrGlc, 2)*pow(KglctrIIG6P, 2)*Vratio*(Glc*Pglctr/KglctrGlc + 1)*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)));
            break;
        case 23:
            dwdp[10] = ATP*Glc*Vhkm*sc*(-Glc - KhkDGlc)/pow(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc, 2);
            break;
        case 24:
            dwdp[10] = -ATP*Glc*KhkATP*Vhkm*sc/pow(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc, 2);
            break;
        case 25:
            dwdp[10] = -pow(ATP, 2)*Glc*Vhkm*sc/pow(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc, 2);
            break;
        case 26:
            dwdp[1] = 2*ACA*EtOH*NAD*V1adh1*pow(V2adh1, 2)*sc*(-ACA*NADH/Keq + EtOH*NAD)/(pow(Kiadh1ACA, 2)*pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2));
            break;
        case 27:
            dwdp[1] = 2*ACA*EtOH*NADH*pow(V1adh1, 2)*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)/(Keq*pow(Kiadh1EtOH, 2)*pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2));
            break;
        case 28:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(ACA*Kadh1NADH*NAD*V1adh1/(Keq*pow(Kiadh1NAD, 2)) - Kadh1EtOH*V2adh1)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2);
            break;
        case 29:
            dwdp[1] = 2*EtOH*Kadh1NAD*NADH*V1adh1*pow(V2adh1, 2)*sc*(-ACA*NADH/Keq + EtOH*NAD)/(pow(Kiadh1NADH, 2)*pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2));
            break;
        case 30:
            dwdp[25] = DHAP*Vlpglycm*sc*(-KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH - 1)/pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2);
            break;
        case 31:
            dwdp[25] = DHAP*Vlpglycm*sc*(DHAP*KlpglycNADH*NAD/(pow(KlpglycINAD, 2)*NADH) + KlpglycDHAP*KlpglycINADH*NAD/(pow(KlpglycINAD, 2)*NADH))/pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2);
            break;
        case 32:
            dwdp[25] = -DHAP*KlpglycDHAP*Vlpglycm*sc*(1 + NAD/KlpglycINAD)/(NADH*pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2));
            break;
        case 33:
            dwdp[25] = -pow(DHAP, 2)*Vlpglycm*sc*(1 + NAD/KlpglycINAD)/(NADH*pow(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1), 2));
            break;
        case 34:
            dwdp[17] = -1.1499999999999999*Pyr*VPYRd*sc/pow(KmPYRd + Pyr, 2);
            break;
        case 36:
            dwdp[13] = ATP*pow(F6P, 2)*Vpfkm*sc*(-1 - pow(ATP, 2)*kappapfk/pow(AMP, 2))/((ATP + KpfkmATP)*pow(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)), 2));
            break;
        case 37:
            dwdp[13] = -ATP*pow(F6P, 2)*Vpfkm*sc/(pow(ATP + KpfkmATP, 2)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2))));
            break;
        case 38:
            dwdp[14] = sc*(-pow(F6P, 2)*KpgiG6P*Vpgi/(pow(KpgiF6P, 2)*Kpgieq*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)) + F6P*G6P*KpgiG6P*Vpgi/(pow(KpgiF6P, 2)*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)));
            break;
        case 39:
            dwdp[14] = sc*(-F6P*Vpgi*(-F6P/KpgiF6P - 1)/(Kpgieq*pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2)) + G6P*Vpgi*(-F6P/KpgiF6P - 1)/pow(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P, 2));
            break;
        case 40:
            dwdp[14] = F6P*Vpgi*sc/(pow(Kpgieq, 2)*(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P));
            break;
        case 41:
            dwdp[15] = -ADP*PEP*Vpkm*sc/(pow(ADP + KpkADP, 2)*(KpkPEP + PEP));
            break;
        case 43:
            dwdp[15] = -ADP*PEP*Vpkm*sc/((ADP + KpkADP)*pow(KpkPEP + PEP, 2));
            break;
        case 44:
            dwdp[19] = sc*(DHAP*Vtimm*(-GAP/KtimGAP - 1)/pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2) - GAP*Vtimm*(-GAP/KtimGAP - 1)/(Ktimeq*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)));
            break;
        case 45:
            dwdp[19] = sc*(DHAP*GAP*KtimDHAP*Vtimm/(pow(KtimGAP, 2)*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)) - pow(GAP, 2)*KtimDHAP*Vtimm/(pow(KtimGAP, 2)*Ktimeq*pow(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP, 2)));
            break;
        case 46:
            dwdp[19] = GAP*Vtimm*sc/(pow(Ktimeq, 2)*(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP));
            break;
        case 47:
            dwdp[12] = 0.69999999999999996*Pyr*V1PDC*sc*(G6P*KaPDC/(pow(KxPDC, 2)*(G6P*bet/(KxPDC*alp) + 1)) - G6P*KaPDC*bet*(G6P/KxPDC + 1)/(pow(KxPDC, 2)*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)) - G6P*Pyr*bet*(G6P/(KxPDC*alp) + 1)/(pow(KxPDC, 2)*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)) + G6P*Pyr/(pow(KxPDC, 2)*alp*(G6P*bet/(KxPDC*alp) + 1)))/pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2);
            break;
        case 49:
            dwdp[8] = sc*(-Glc*Vglctr*(-Glc*(Glc0/KglctrGlc + 1)/(KglctrGlc*(Glc0*Pglctr/KglctrGlc + 1)) + Glc0*(Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(KglctrGlc*pow(Glc0*Pglctr/KglctrGlc + 1, 2)))/(KglctrGlc*Vratio*pow(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1, 2)) + Glc0*Vglctr*(Glc*(Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(KglctrGlc*pow(Glc*Pglctr/KglctrGlc + 1, 2)) - Glc0*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(KglctrGlc*(Glc*Pglctr/KglctrGlc + 1)))/(KglctrGlc*Vratio*pow(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1), 2)));
            break;
        case 50:
            dwdp[12] = 0.69999999999999996*Pyr*sc/(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1));
            break;
        case 51:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*EtOH*NADH/(Keq*Kiadh1EtOH) - ACA*Kadh1NADH/Keq - ACA*Kadh1NADH*NAD/(Keq*Kiadh1NAD) - ACA*NADH/Keq - Kadh1ACA*NADH/Keq)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2) + 2*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1);
            break;
        case 52:
            dwdp[1] = 2*V1adh1*V2adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)*(-ACA*EtOH*NAD/Kiadh1ACA - EtOH*Kadh1NAD - EtOH*Kadh1NAD*NADH/Kiadh1NADH - EtOH*NAD - Kadh1EtOH*Kiadh1NAD - Kadh1EtOH*NAD)/pow(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1, 2) + 2*V1adh1*sc*(-ACA*NADH/Keq + EtOH*NAD)/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1);
            break;
        case 53:
            dwdp[17] = 1.1499999999999999*Pyr*sc/(KmPYRd + Pyr);
            break;
        case 54:
            dwdp[3] = sc*(-DHAP*GAP*Valdf*(-DHAP*GAP/(Kaldeq*Valdr) - DHAP*KaldGAP/(Kaldeq*Valdr) - GAP*KaldDHAP/(Kaldeq*Valdr))/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) - DHAP*GAP/(Kaldeq*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)) + FBP*Valdf*(-DHAP*GAP/(Kaldeq*Valdr) - DHAP*KaldGAP/(Kaldeq*Valdr) - GAP*KaldDHAP/(Kaldeq*Valdr))/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2) + FBP/(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP));
            break;
        case 55:
            dwdp[3] = sc*(-DHAP*GAP*Valdf*(DHAP*GAP*Valdf/(Kaldeq*pow(Valdr, 2)) + DHAP*KaldGAP*Valdf/(Kaldeq*pow(Valdr, 2)) + GAP*KaldDHAP*Valdf/(Kaldeq*pow(Valdr, 2)))/(Kaldeq*pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2)) + FBP*Valdf*(DHAP*GAP*Valdf/(Kaldeq*pow(Valdr, 2)) + DHAP*KaldGAP*Valdf/(Kaldeq*pow(Valdr, 2)) + GAP*KaldDHAP*Valdf/(Kaldeq*pow(Valdr, 2)))/pow(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP, 2));
            break;
        case 56:
            dwdp[5] = 0.83199999999999996*ATP*sc/(ATP + KatpconATP);
            break;
        case 57:
            dwdp[7] = sc*(-DPG*NADH/(KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*NAD/(KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)));
            break;
        case 58:
            dwdp[8] = sc*(-Glc/(KglctrGlc*Vratio*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)) + Glc0/(KglctrGlc*Vratio*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1))));
            break;
        case 59:
            dwdp[10] = ATP*Glc*sc/(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc);
            break;
        case 60:
            dwdp[25] = DHAP*sc/(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1));
            break;
        case 62:
            dwdp[13] = ATP*pow(F6P, 2)*sc/((ATP + KpfkmATP)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2))));
            break;
        case 63:
            dwdp[14] = sc*(-F6P/(Kpgieq*(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P)) + G6P/(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P));
            break;
        case 64:
            dwdp[15] = ADP*PEP*sc/((ADP + KpkADP)*(KpkPEP + PEP));
            break;
        case 65:
            dwdp[8] = sc*(Glc*Vglctr/(KglctrGlc*pow(Vratio, 2)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)) - Glc0*Vglctr/(KglctrGlc*pow(Vratio, 2)*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1))));
            dwdp[20] = sc*(-ACA*dACA/pow(Vratio, 2) + ACA0*dACA/pow(Vratio, 2));
            dwdp[21] = sc*(-EtOH*dEtOH/pow(Vratio, 2) + EtOH0*dEtOH/pow(Vratio, 2));
            dwdp[22] = sc*(-Glyc*dGlyc/pow(Vratio, 2) + Glyc0*dGlyc/pow(Vratio, 2));
            dwdp[23] = sc*(-HCN*dHCN/pow(Vratio, 2) + HCN0*dHCN/pow(Vratio, 2));
            dwdp[24] = sc*(-OAc*dOAc/pow(Vratio, 2) + OAc0*dOAc/pow(Vratio, 2));
            break;
        case 66:
            dwdp[19] = sc*(DHAP/(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP) - GAP/(Ktimeq*(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP)));
            break;
        case 67:
            dwdp[12] = 0.69999999999999996*Pyr*V1PDC*sc*(-G6P*KaPDC*bet*(G6P/KxPDC + 1)/(KxPDC*pow(alp, 2)*pow(G6P*bet/(KxPDC*alp) + 1, 2)) - G6P*Pyr*bet*(G6P/(KxPDC*alp) + 1)/(KxPDC*pow(alp, 2)*pow(G6P*bet/(KxPDC*alp) + 1, 2)) + G6P*Pyr/(KxPDC*pow(alp, 2)*(G6P*bet/(KxPDC*alp) + 1)))/pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2);
            break;
        case 68:
            dwdp[12] = 0.69999999999999996*Pyr*V1PDC*sc*(G6P*KaPDC*(G6P/KxPDC + 1)/(KxPDC*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)) + G6P*Pyr*(G6P/(KxPDC*alp) + 1)/(KxPDC*alp*pow(G6P*bet/(KxPDC*alp) + 1, 2)))/pow(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1), 2);
            break;
        case 69:
            dwdp[20] = sc*(ACA/Vratio - ACA0/Vratio);
            break;
        case 70:
            dwdp[21] = sc*(EtOH/Vratio - EtOH0/Vratio);
            break;
        case 71:
            dwdp[22] = sc*(Glyc/Vratio - Glyc0/Vratio);
            break;
        case 72:
            dwdp[23] = sc*(HCN/Vratio - HCN0/Vratio);
            break;
        case 73:
            dwdp[24] = sc*(OAc/Vratio - OAc0/Vratio);
            break;
        case 74:
            dwdp[4] = ACA*NAD*sc;
            break;
        case 75:
            dwdp[0] = 0.5*ACA0;
            dwdp[9] = HCN0;
            break;
        case 76:
            dwdp[11] = -lacto;
            break;
        case 77:
            dwdp[11] = ACA0*HCN0;
            break;
        case 78:
            dwdp[2] = AMP*ATP*sc;
            break;
        case 79:
            dwdp[2] = -pow(ADP, 2)*sc;
            break;
        case 80:
            dwdp[13] = -pow(ATP, 3)*pow(F6P, 2)*Kpfk*Vpfkm*sc/(pow(AMP, 2)*(ATP + KpfkmATP)*pow(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2)), 2));
            break;
        case 81:
            dwdp[6] = -DHAPCN;
            break;
        case 82:
            dwdp[6] = DHAP*HCN;
            break;
        case 83:
            dwdp[18] = ATP*G6P*sc;
            break;
        case 84:
            dwdp[26] = ADP*DPG*sc;
            break;
        case 85:
            dwdp[26] = -ATP*PEP*sc;
            break;
        case 86:
            dwdp[16] = -PyrCN;
            break;
        case 87:
            dwdp[16] = HCN*Pyr;
            break;
        case 88:
            dwdp[1] = 2*V1adh1*V2adh1*(-ACA*NADH/Keq + EtOH*NAD)/(ACA*EtOH*NAD*V2adh1/Kiadh1ACA + ACA*EtOH*NADH*V1adh1/(Keq*Kiadh1EtOH) + ACA*Kadh1NADH*V1adh1/Keq + ACA*Kadh1NADH*NAD*V1adh1/(Keq*Kiadh1NAD) + ACA*NADH*V1adh1/Keq + EtOH*Kadh1NAD*V2adh1 + EtOH*Kadh1NAD*NADH*V2adh1/Kiadh1NADH + EtOH*NAD*V2adh1 + Kadh1ACA*NADH*V1adh1/Keq + Kadh1EtOH*Kiadh1NAD*V2adh1 + Kadh1EtOH*NAD*V2adh1);
            dwdp[2] = -pow(ADP, 2)*kakr + AMP*ATP*kakf;
            dwdp[3] = -DHAP*GAP*Valdf/(Kaldeq*(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP)) + FBP*Valdf/(DHAP*GAP*Valdf/(Kaldeq*Valdr) + DHAP*KaldGAP*Valdf/(Kaldeq*Valdr) + FBP*GAP/KaldIGAP + FBP + GAP*KaldDHAP*Valdf/(Kaldeq*Valdr) + KaldFBP);
            dwdp[4] = ACA*NAD*kOAc;
            dwdp[5] = 0.83199999999999996*ATP*Vatpconm/(ATP + KatpconATP);
            dwdp[7] = -DPG*NADH*Vgapdh/(KgapdhGAP*KgapdhNAD*Kgapdheq*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1)) + GAP*NAD*Vgapdh/(KgapdhGAP*KgapdhNAD*(1 + NADH/KgapdhNADH + NAD/KgapdhNAD)*(DPG/KgapdhDPG + GAP/KgapdhGAP + 1));
            dwdp[8] = -Glc*Vglctr/(KglctrGlc*Vratio*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + (Glc0/KglctrGlc + 1)*(Glc*Pglctr/KglctrGlc + 1)/(Glc0*Pglctr/KglctrGlc + 1) + 1)) + Glc0*Vglctr/(KglctrGlc*Vratio*(Glc0/KglctrGlc + 1 + (Glc0*Pglctr/KglctrGlc + 1)*(G6P*Glc/(KglctrGlc*KglctrIIG6P) + G6P/KglctrIG6P + Glc/KglctrGlc + 1)/(Glc*Pglctr/KglctrGlc + 1)));
            dwdp[10] = ATP*Glc*Vhkm/(ATP*Glc + ATP*KhkGlc + Glc*KhkATP + KhkATP*KhkDGlc);
            dwdp[12] = 0.69999999999999996*Pyr*V1PDC/(KaPDC*(G6P/KxPDC + 1)/(G6P*bet/(KxPDC*alp) + 1) + Pyr*(G6P/(KxPDC*alp) + 1)/(G6P*bet/(KxPDC*alp) + 1));
            dwdp[13] = ATP*pow(F6P, 2)*Vpfkm/((ATP + KpfkmATP)*(pow(F6P, 2) + Kpfk*(1 + pow(ATP, 2)*kappapfk/pow(AMP, 2))));
            dwdp[14] = -F6P*Vpgi/(Kpgieq*(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P)) + G6P*Vpgi/(F6P*KpgiG6P/KpgiF6P + G6P + KpgiG6P);
            dwdp[15] = ADP*PEP*Vpkm/((ADP + KpkADP)*(KpkPEP + PEP));
            dwdp[17] = 1.1499999999999999*Pyr*VPYRd/(KmPYRd + Pyr);
            dwdp[18] = ATP*G6P*kg6pst;
            dwdp[19] = DHAP*Vtimm/(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP) - GAP*Vtimm/(Ktimeq*(DHAP + GAP*KtimDHAP/KtimGAP + KtimDHAP));
            dwdp[20] = ACA*dACA/Vratio - ACA0*dACA/Vratio;
            dwdp[21] = EtOH*dEtOH/Vratio - EtOH0*dEtOH/Vratio;
            dwdp[22] = Glyc*dGlyc/Vratio - Glyc0*dGlyc/Vratio;
            dwdp[23] = HCN*dHCN/Vratio - HCN0*dHCN/Vratio;
            dwdp[24] = OAc*dOAc/Vratio - OAc0*dOAc/Vratio;
            dwdp[25] = DHAP*Vlpglycm/(DHAP*(KlpglycNADH*(1 + NAD/KlpglycINAD)/NADH + 1) + KlpglycDHAP*(KlpglycINADH*(1 + NAD/KlpglycINAD)/NADH + 1));
            dwdp[26] = ADP*DPG*klppepf - ATP*PEP*klppepr;
            break;
    }
}