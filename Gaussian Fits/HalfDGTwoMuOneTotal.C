#pragma once

#include "TMath.h"
#include <iostream>

Double_t HalfDGTwoMuOneTotal(Double_t *v, Double_t *par) {

//First Gaussian function
//Par[0]=ScalingFactor, Par[1]=TotalSignal, Par[2]=MeanValue, Par[3]=StandardDeviation
    Double_t arg1 = 0;
    if (par[3] != 0)
        arg1 = (v[0] - par[2])/par[3];

    Double_t fitval1 = 0;
    if (par[0] < 1. && par[0] > 0.)
        fitval1 = par[0]*par[1]*TMath::Exp(-0.5*arg1*arg1);

    fitval1 = (fitval1*(0.5))/(TMath::Sqrt(TMath::TwoPi())*par[3]);

//Second Gaussian function (Same Mean Value)
//Par[4]=MeanValue, Par[5]=StandardDeviation
    Double_t arg2 = 0;
    if (par[5] != 0)
        arg2 = (v[0] - par[4])/par[5];

    Double_t fitval2 = 0;
    if (par[0] < 1. && par[0] > 0.)
        fitval2 = (1 - par[0])*par[1]*TMath::Exp(-0.5*arg2*arg2);

fitval2 = (fitval2*(0.5))/(TMath::Sqrt(TMath::TwoPi())*par[5]);

//Par[5]=LinearBackgroundIntercept, Par[6]=BackgroundSlope
//Adding Gaussians to linear background
Double_t fitval = fitval1 +fitval2 + par[6] + v[0]*par[7];
  Double_t MuAVG = par[0]*par[2] + (1-par[0])*par[4];
    return fitval;
}
