#pragma once

#include "evolvent.h"

TEvolvent* myClassEvovlent;
double result_image[2];

extern "C" void CreateEvovlent();
extern "C" void SetBounds(double A0, double A1, double B0, double B1);
extern "C" double GetImage0(double x);
extern "C" double GetImage1(double x);
extern "C" void FreeEvolvent();