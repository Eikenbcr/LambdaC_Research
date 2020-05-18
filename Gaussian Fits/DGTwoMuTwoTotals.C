#pragma once

#include "TMath.h"
#include <iostream>

Double_t DGTwoMuTwoTotal(Double_t *v, Double_t *par) {

//First Gaussian function
//Par[0]=SignalEvents, Par[1]=MeanValue, Par[2]=StandardDeviation
    Double_t arg1 = 0;
    if (par[2] != 0)
        arg1 = (v[0] - par[1])/par[2];

    Double_t fitval1 = par[0]*TMath::Exp(-0.5*arg1*arg1);

    fitval1 = (fitval1*2.00)/(TMath::Sqrt(TMath::TwoPi())*par[2]);

//Second Gaussian function (Same Mean Value)
//Par[3]=SignalEvents, Par[4]=MeanValue, Par[5]=StandardDeviation
Double_t arg2 = 0;
if (par[5] != 0)
    arg2 = (v[0] - par[4])/par[5];

Double_t fitval2 = par[3]*TMath::Exp(-0.5*arg2*arg2);

fitval2 = (fitval2*2.00)/(TMath::Sqrt(TMath::TwoPi())*par[5]);

//Par[5]=LinearBackgroundIntercept, Par[6]=BackgroundSlope
//Adding Gaussians to linear background
Double_t fitval = fitval1 +fitval2 + par[6] + v[0]*par[7];

    return fitval;
}
