#pragma once

#include "TMath.h"
#include <iostream>

Double_t DGOneMuOneTotal(Double_t *v, Double_t *par) {

//First Gaussian function
//Par[0]=ScalingFactor, Par[1]=TotalSignal, Par[2]=MeanValue, Par[3]=StandardDeviation
    Double_t arg1 = 0;
    if (par[3] != 0)
        arg1 = (v[0] - par[2])/par[3];

    Double_t fitval1 = 0;
    if (par[0] < 1.)
        fitval1 = par[0]*par[1]*TMath::Exp(-0.5*arg1*arg1);

    fitval1 = (fitval1*2.00)/(TMath::Sqrt(TMath::TwoPi())*par[3]);

//Second Gaussian function (Same Mean Value)
//Par[2]=MeanValue, Par[4]=StandardDeviation
    Double_t arg2 = 0;
    if (par[4] != 0)
        arg2 = (v[0] - par[2])/par[4];

    Double_t fitval2 = 0;
    if (par[0] < 1.)
        fitval2 = (1 - par[0])*par[1]*TMath::Exp(-0.5*arg2*arg2);

fitval2 = (fitval2*2.00)/(TMath::Sqrt(TMath::TwoPi())*par[4]);

//Par[5]=LinearBackgroundIntercept, Par[6]=BackgroundSlope
//Adding Gaussians to linear background
Double_t fitval = fitval1 +fitval2 + par[5] + v[0]*par[6];

    return fitval;
}
