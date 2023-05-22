#pragma once

#include "evolvent.h"

struct TThread{
    TEvolvent* myClassEvovlent;
    double result_image[2];
};

TThread threads[4];

extern "C" void CreateEvovlent(int idx);
extern "C" void SetBounds(int idx, double A0, double A1, double B0, double B1);
extern "C" double GetImage0(int idx, double x);
extern "C" double GetImage1(int idx, double x);
extern "C" void FreeEvolvent(int idx);