#pragma once

#include "grishagin_function.h"

vagrish::GrishaginFunction* myClassGrishaginFunction;

extern "C" void CreateGrishaginFunction();
extern "C" void MySetFunctionNumber(int i);
extern "C" double GetBoundsA0();
extern "C" double GetBoundsA1();
extern "C" double GetBoundsB0();
extern "C" double GetBoundsB1();
extern "C" double MyCalculate(double y0, double y1);
extern "C" void FreeGrishaginFunction();