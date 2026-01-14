//  Element.cpp
//
//  Created by Lele Shu on 7/17/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "Element.hpp"
void Triangle::printHeader(FILE *fp){
    for(int i = 0; i < 3; i++){
        fprintf(fp, "node%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "nabr%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "edge%d\t", i);
    }
    fprintf(fp, "%s\t", "area");
    fprintf(fp, "%s\t", "x");
    fprintf(fp, "%s\t", "y");
    fprintf(fp, "%s\t", "zmin");
    fprintf(fp, "%s\t", "zmax");
}
void Triangle::printInfo(FILE *fp){
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%d\t", node[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%d\t", nabr[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", edge[i]);
    }
    fprintf(fp, "%g\t", area);
    fprintf(fp, "%g\t", x);
    fprintf(fp, "%g\t", y);
    fprintf(fp, "%g\t", z_bottom);
    fprintf(fp, "%g\t", z_surf);
}
void AttriuteIndex::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "iSoil");
    fprintf(fp, "%s\t", "iGeol");
    fprintf(fp, "%s\t", "iLC");
    fprintf(fp, "%s\t", "IC");
    fprintf(fp, "%s\t", "iForc");
    fprintf(fp, "%s\t", "iMF");
    fprintf(fp, "%s\t", "iBC");
    fprintf(fp, "%s\t", "iLake");
}
void AttriuteIndex::printInfo(FILE *fp){
    fprintf(fp, "%d\t", iSoil);
    fprintf(fp, "%d\t", iGeol);
    fprintf(fp, "%d\t", iLC);
    fprintf(fp, "%d\t", IC);
    fprintf(fp, "%d\t", iForc);
    fprintf(fp, "%d\t", iMF);
    fprintf(fp, "%d\t", iBC);
    fprintf(fp, "%d\t", iSS);
    fprintf(fp, "%d\t", iLake);
}

void _Element::applyGeometry(_Node *Node){
    double  x1, y1, x2, y2, x3, y3;
    double  zmin1, zmax1, zmin2, zmax2, zmin3, zmax3;
    double  aqd, z0;
    double  d1, d2, d3;
    double  px1, py1, pz1,
            px2, py2, pz2,
            px3, py3, pz3;
    x1 = Node[node[0] - 1].x;
    x2 = Node[node[1] - 1].x;
    x3 = Node[node[2] - 1].x;
    y1 = Node[node[0] - 1].y;
    y2 = Node[node[1] - 1].y;
    y3 = Node[node[2] - 1].y;
    
    zmin1 = Node[node[0] - 1].zmin;
    zmin2 = Node[node[1] - 1].zmin;
    zmin3 = Node[node[2] - 1].zmin;
    zmax1 = Node[node[0] - 1].zmax;
    zmax2 = Node[node[1] - 1].zmax;
    zmax3 = Node[node[2] - 1].zmax;
    
    area = 0.5 * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1));
    
    z_surf = (zmax1 + zmax2 + zmax3) / 3.0;
    z_bottom = (zmin1 + zmin2 + zmin3) / 3.0;
    if(z_bottom >= z_surf){
        printf("WARNING: zmin(%f) >= zmax(%f) at element %d.", z_bottom, z_surf, index);
    }
    if(zcentroid != NA_VALUE){
        z0 = z_surf;
        aqd = z_surf - z_bottom;
        z_surf = (zmax1 + zmax2 + zmax3) / 3.0;
        z_bottom = z_surf - aqd;
        if(fabs(z0 - z_surf) > 10.){
//            printf("DZ(%d) = %f\n", index, z0 - zmax);
        }
    }
    /* calculate centroid of triangle */
    x = (x1 + x2 + x3) / 3.0;
    y = (y1 + y2 + y3) / 3.0;
    edge[0] = Eudist(x2, y2, x3, y3);
    edge[1] = Eudist(x3, y3, x1, y1);
    edge[2] = Eudist(x1, y1, x2, y2);
    
    PointPerpdicularOnLine(&px1, &py1, x, y, x2, y2, x3, y3);
    PointPerpdicularOnLine(&px2, &py2, x, y, x3, y3, x1, y1);
    PointPerpdicularOnLine(&px3, &py3, x, y, x1, y1, x2, y2);
    d1 = Eudist(px1, py1, x, y);
    d2 = Eudist(px2, py2, x, y);
    d3 = Eudist(px3, py3, x, y);
    
    Dist2Edge[0] = d1;
    Dist2Edge[1] = d2;
    Dist2Edge[2] = d3;
//    Dist2Edge[0] = sqpow2(edge[0] * edge[1] * edge[2] / (4 * area), edge[0] / 2 );
//    Dist2Edge[1] = sqpow2(edge[0] * edge[1] * edge[2] / (4 * area), edge[1] / 2 );
//    Dist2Edge[2] = sqpow2(edge[0] * edge[1] * edge[2] / (4 * area), edge[2] / 2 );
    
    pz1 = ZOnLine(x2, y2, zmax2, x3, y3, zmax3, px1, py1);
    pz2 = ZOnLine(x3, y3, zmax3, x1, y1, zmax1, px2, py2);
    pz3 = ZOnLine(x1, y1, zmax1, x2, y2, zmax2, px3, py3);
    slope[0] =  (z_surf - pz1)/ d1;
    slope[1] =  (z_surf - pz2)/ d1;
    slope[2] =  (z_surf - pz3)/ d1;

    /*
     * ================================================================
     * Terrain geometry (normal vector / slope angle / aspect)
     * ------------------------------------------------
     * 计算步骤（与需求一一对应）：
     * 1) 获取三个节点坐标 (x, y, z)：这里使用地表高程 zmax 作为地形面。
     * 2) 计算两条边向量 v1, v2
     * 3) 叉积得到法向量 n = v1 × v2
     * 4) 单位化并确保 nz >= 0（若 nz < 0 则翻转法向量）
     * 5) slopeAngle = atan2(hypot(nx,ny), nz)
     * 6) aspect = atan2(nx, ny)，并调整到 [0, 2*pi)
     *    这里采用 North=0（+y 方向为北），East=pi/2 的定义：
     *      - 若法向量水平投影指向 +y，则 aspect=0
     *      - 若法向量水平投影指向 +x，则 aspect=pi/2
     * 7) 水平面（几乎无坡度）时 aspect = 0（用 slopeAngle 阈值判断以避免近水平面噪声）
     *
     * 注：
     * - 若三点共线（退化三角形），叉积长度为 0，此时回退为水平面：
     *   (nx, ny, nz) = (0, 0, 1)，slopeAngle=0，aspect=0。
     * ================================================================
     */
    {
        /* Step 1: three node coordinates (surface points) */
        const double p1x = x1, p1y = y1, p1z = zmax1;
        const double p2x = x2, p2y = y2, p2z = zmax2;
        const double p3x = x3, p3y = y3, p3z = zmax3;

        /* Step 2: edge vectors v1 = p2 - p1, v2 = p3 - p1 */
        const double v1x = p2x - p1x;
        const double v1y = p2y - p1y;
        const double v1z = p2z - p1z;
        const double v2x = p3x - p1x;
        const double v2y = p3y - p1y;
        const double v2z = p3z - p1z;

        /* Step 3: cross product n = v1 x v2 */
        double nx_raw = v1y * v2z - v1z * v2y;
        double ny_raw = v1z * v2x - v1x * v2z;
        double nz_raw = v1x * v2y - v1y * v2x;

        /* Step 4: normalize and ensure nz >= 0 */
        const double nlen = sqrt(nx_raw * nx_raw + ny_raw * ny_raw + nz_raw * nz_raw);
        if (nlen <= ZERO) {
            /* Degenerate triangle (or extremely small): fall back to flat terrain */
            nx = 0.0;
            ny = 0.0;
            nz = 1.0;
        } else {
            nx = nx_raw / nlen;
            ny = ny_raw / nlen;
            nz = nz_raw / nlen;
            if (nz < 0.0) {
                nx = -nx;
                ny = -ny;
                nz = -nz;
            }
        }

        /*
         * Step 5: slopeAngle = atan2(hypot(nx,ny), nz)
         *
         * 相比 acos(nz)，atan2 在 nz≈1（近水平面、小坡度）时数值更稳定：
         * acos(1-ε) 会因为浮点精度导致有效位损失，从而引入 slopeAngle 噪声。
         * 这里仍对 nz 做 [0,1] 截断以避免归一化误差带来的越界。
         */
        const double nz_clamped = min(1.0, max(0.0, nz));
        slopeAngle = atan2(hypot(nx, ny), nz_clamped);

        /*
         * Step 7: horizontal plane (or near horizontal) => aspect = 0
         *
         * 近水平面时，nx/ny 的微小数值噪声会被 atan2 放大，导致 aspect 在 [0,2π)
         * 内抖动。与其使用过小的 nx/ny 阈值（如 1e-12），这里改为基于坡度角判断：
         * 当 slopeAngle < 1e-6 时视为水平面并固定 aspect=0。
         */
        const double aspect_flat_slope_eps = 1e-6;
        if (slopeAngle < aspect_flat_slope_eps) {
            aspect = 0.0;
        } else {
            /* Step 6: aspect = atan2(nx, ny) and map to [0, 2*pi) */
            aspect = atan2(nx, ny);
            if (aspect < 0.0) {
                aspect += 2.0 * PI;
            }
            if (aspect >= 2.0 * PI) {
                aspect -= 2.0 * PI;
            }
        }
    }
}
void _Element::InitElement(){
    AquiferDepth = z_surf - z_bottom;
    WetlandLevel = AquiferDepth - infD;
    RootReachLevel = AquiferDepth - RzD;
    FixPressure = PressureElevation(z_surf); // P atmospheric pressure [kPa]
//    FixGamma = PsychrometricConstant(FixPressure);
    MacporeLevel = AquiferDepth - macD;
    
    if (AquiferDepth < macD)
        macD = AquiferDepth;
    
    CheckNonZero(RootReachLevel, index-1, "RootReachLevel");
    CheckNonZero(FixPressure, index-1, "FixPressure");
    CheckNonZero(area, index-1, "area");
    CheckNA(z_surf, "zmax");
    CheckNA(z_bottom, "zmin");
    for(int i = 0; i < 3; i++){
        CheckNonZero(edge[i], index-1, "edge");
    }
}
void _Element::applyNabor(_Node *Node, _Element *Ele){
    for (int j = 0; j < 3; j++) {
        if(nabr[j]>0){
            for(int k = 0; k < 3; k++){
                if(Ele[nabr[j] - 1].nabr[k] == this->index){
                    this->nabrToMe[j] = k + 1;
//                    printf("%d, %d <-> %d, %d\n", this->index, j+1, nabr[j], k+1);
                    /* Neighbor's K direction and My J direction shared a edge*/
                }
            }
        }else{
            /* VOID */
        }
    }

    for(int j = 0; j < 3; j++){
//        Dist2Edge[j] = sqpow2(edge[0] * edge[1] * edge[2] / (4 * area), edge[j] / 2 );
        //          Dist2Edge[j] =  sqrt(pow(edge[0] * edge[1] * edge[2] / (4 * area), 2) - pow(edge[j] / 2, 2));
        if(nabr[j] > 0){ /* Neighbor exist */
            Dist2Nabor[j] = Eudist(x, y,
                                   Ele[nabr[j] - 1].x, Ele[nabr[j] - 1].y);
            avgRough[j] = 0.5 * (Rough + Ele[nabr[j] - 1].Rough);
        }else if(nabr[j] < 0){ /* Next to LAKE Element*/
            Dist2Nabor[j] = Dist2Edge[j];
            avgRough[j] = Rough;
        }else{  /* No Neighbor, Not Lake element*/
            Dist2Nabor[j] = 0.;
            avgRough[j] = Rough;
        }
        //            Dist2Nabor[j] = sqrt(pow((x - Ele[eNabor].x), 2) + pow((y - Ele[eNabor].y), 2));
    }
    
}
void _Element::Flux_Infiltration(double Ysurf, double Yunsat, double Ygw, double netprcp){
    double av = Ysurf + netprcp, grad=0;
    if(Ygw + Yunsat > AquiferDepth || u_deficit < Yunsat){
        /* GW reach the surface */
        u_qex =  fabs(Ygw + Yunsat - AquiferDepth) / AquiferDepth * Kmax; /* Exfiltration must be positive (upward). */
        u_qi = 0. ;
    }else{
        /* GW level is lower than Surface */
        u_qex = 0.;
        if(av > 0. && u_deficit > infD){
            grad = 1. + av / infD;
//            grad = (av -  u_phius )/ infD;
            if( av > Kmax){
                /* Heavy rainfall, Macropore works */
                u_effkInfi = infKsatV * (1 - hAreaF) + hAreaF * macKsatV * u_satn;
            }else if( av > infKsatV ){
                /* Medium rainfall, Macropore works */
                u_effkInfi = u_satKr * infKsatV * (1 - hAreaF) + hAreaF * macKsatV * u_satn;
            }else{
                /* Light rainfall */
                u_effkInfi = u_satKr * infKsatV * (1 - hAreaF);
            }
            u_qi = grad * u_effkInfi;
            u_qi = min(av, max(0., u_qi) );
        }else{
            u_qi = 0;
        }
    }
//    if(u_qi > 1.){ /* DEBUG ONLY*/
//        u_qi=u_qi;
//    }
//    CheckNANi(u_qi, 0, "u_qi");
}
double _Element::Flux_Recharge(double Yunsat, double Ygw){
    double ke=0., grad, ku;
    if(Ygw > AquiferDepth - infD && Yunsat < u_deficit){
        u_qr  =  0.;
        return u_qr;
    }
    if (u_theta > ThetaR) {
        if(Yunsat <= EPSILON){
            grad = 0.;
        }else{
            grad = (u_theta - ThetaR) / (ThetaFC - ThetaR);
//            grad = (u_theta - ThetaFC) / (ThetaS - ThetaFC);
//            grad = (0.5 * u_deficit + u_phius) / (0.5 * u_deficit);
            grad = max(grad, 0.);
        }
    }else{
        /* DRY condition. */
        grad = 0.;
    }
    if( infKsatV <= 0. || KsatV <= 0.){
        u_qr = 0.;
    }else{
        ku = infKsatV * u_satKr; //max(Yunsat / u_deficit, 0.);
        ke = meanHarmonic(ku, KsatV, u_deficit, Ygw);
//        ke = meanArithmetic(ku, KsatV, u_deficit, Ygw);
        u_qr = grad * ke;
    }
#ifdef DEBUG
    CheckNANi(u_qr, 0, "_Element::Flux_Recharge():u_qr");
#endif
    return u_qr;
}
void _Element::updateLakeElement(){
    u_effKH = KsatH;
    u_deficit = 0;
    Kmax = infKsatV;
    u_deficit = 0.;
    u_satn = 1.;
    u_theta = ThetaS;
    u_satKr = 1.0;
    u_phius = 0.;
    u_effkInfi = infKsatV;
}
void _Element::updateElement(double Ysurf, double Yunsat, double Ygw){
    u_effKH = effKH(Ygw,  AquiferDepth,  macD,  macKsatH,  geo_vAreaF,  KsatH);
    u_deficit = AquiferDepth - Ygw;
    Kmax = infKsatV * (1. - hAreaF) + macKsatV * hAreaF ;
    if(u_deficit <= 0. ){
        u_deficit = 0.;
        u_satn = 1.;
        u_theta = ThetaS;
    }else{
        u_theta = Yunsat / u_deficit * ThetaS;
        u_satn = (u_theta - ThetaR) / (ThetaS - ThetaR) ;
    }
    if(u_satn > 0.99 ){
        u_satn = 1.0;
        u_satKr = 1.0;
        u_phius = 0.;
        u_theta = ThetaS;
    }else if(u_satn <= ZERO){
        u_satn = 0.;
        u_satKr = 0.;
        u_phius = MINPSI;
        u_theta = ThetaR;
    }else{
        u_satKr = satKfun(u_satn, Beta);
        u_phius = sat2psi(u_satn, Alpha, Beta);
        u_phius = max(MINPSI, u_phius);
    }
    u_effkInfi = infKsatV * (1 - hAreaF) + u_satn * macKsatV * hAreaF ;
#ifdef DEBUG
//    CheckNANi(u_satKr, 0, "_Element::updateElement():u_satKr");
//    if (u_effkInfi < ZERO){
//        printf("WARNING: Negative effective conductivity for infiltration.\n");
//    }
//    if (u_effKH < ZERO){
//        printf("WARNING: Negative effective conductivity for infiltration.\n");
//    }
#endif
}

void _Element::copyGeol(Geol_Layer *g){
    KsatH        = g[iGeol - 1].KsatH;
    KsatV        = g[iGeol - 1].KsatV;
    geo_ThetaS    = g[iGeol - 1].geo_ThetaS;
    geo_ThetaR    = g[iGeol - 1].geo_ThetaR;
    geo_vAreaF    = g[iGeol - 1].geo_vAreaF;
    macKsatH    = g[iGeol - 1].macKsatH;
    macD        = g[iGeol - 1].macD;
    Sy    = g[iGeol - 1].Sy;
}
void _Element::copySoil(Soil_Layer *g){
    infKsatV = g[iSoil - 1].infKsatV;
    ThetaS = g[iSoil - 1].ThetaS;
    ThetaFC = ThetaS * FieldCapacityRatio;
    ThetaR = g[iSoil - 1].ThetaR;
    Alpha = g[iSoil - 1].Alpha;
    Beta = g[iSoil - 1].Beta;
    hAreaF = g[iSoil - 1].hAreaF;
    macKsatV = g[iSoil - 1].macKsatV;
    infD = g[iSoil - 1].infD;
    
    CheckNonZero(ThetaS, index-1, "ThetaS");
    CheckNonNegative(ThetaR, index-1, "ThetaR");
    CheckNonZero(infD, index-1, "infD");
}
void _Element::copyLandc(Landcover *g){
//    LAImax   = g[iLC - 1].LAImax;
    VegFrac  = g[iLC - 1].VegFrac;
    Albedo   = g[iLC - 1].Albedo;
//    Rs_ref   = g[iLC - 1].Rs_ref;
//    Rmin     = g[iLC - 1].Rmin;
    Rough    = g[iLC - 1].Rough;
    RzD      = g[iLC - 1].RzD;
    SoilDgrd = g[iLC - 1].SoilDgrd;
    ImpAF    = g[iLC - 1].ImpAF;
}
//void _Element::updateWF(double qi, double dt){
//    double dh = 0.;
//    double grad;
//    double effk;
//    if( qi > 0. ){
//        /* Infiltration occur
//         wf increases, bottom loss
//         */
//        dh += qi;
//        if( u_wf > 0.){
//            /* available water in wf */
//            effk = u_satKr * infKsatV;
//            if(u_deficit - u_wf <= 0){
//                grad = 0;
//            }else{
//                grad = (u_deficit * 0.5 - u_phius) / ( 0.5 * (u_deficit - u_wf) );
//            }
//            dh += -effk * grad;
//        }
//        u_wf +=  dh * dt;
//        if(u_wf > 0.01){
//            dh = dh;
//        }
//    }else{
//        /* No infiltration: sat layer moves forward */
//        u_wf = infD;
//    }
//    CheckNA(u_wf, "Weting Front");
//}
void _Element::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "index");
    AttriuteIndex::printHeader(fp);
    Triangle::printHeader(fp);
    Soil_Layer::printHeader(fp);
    Geol_Layer::printHeader(fp);
    Landcover::printHeader(fp);
    for(int i = 0; i < 3; i++){
        fprintf(fp, "Dist2Edge%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "Dist2Nabor%d\t", i);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "avgRough%d\t", i);
    }
    fprintf(fp, "%s\t", "windH");
    fprintf(fp, "%s\t", "FixPressure");
    fprintf(fp, "%s\t", "AquiferDepth");
    fprintf(fp, "%s\t", "WetlandLevel");
    fprintf(fp, "%s\t", "RootReachLevel");
    fprintf(fp, "%s\t", "MacporeLevel");
    fprintf(fp, "\n");
}
void _Element::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
    AttriuteIndex::printInfo(fp);
    Triangle::printInfo(fp);
    Soil_Layer::printInfo(fp);
    Geol_Layer::printInfo(fp);
    Landcover::printInfo(fp);
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", Dist2Edge[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", Dist2Nabor[i]);
    }
    for(int i = 0; i < 3; i++){
        fprintf(fp, "%g\t", avgRough[i]);
    }
    fprintf(fp, "%g\t", windH);
    fprintf(fp, "%g\t", FixPressure);
    fprintf(fp, "%g\t", AquiferDepth);
    fprintf(fp, "%g\t", WetlandLevel);
    fprintf(fp, "%g\t", RootReachLevel);
    fprintf(fp, "%g\t", MacporeLevel);
    fprintf(fp, "\n");
}
