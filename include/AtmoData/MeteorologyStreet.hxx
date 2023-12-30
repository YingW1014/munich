// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Alice Maison - 02/2022
//
// This file is part of AtmoData library, a tool for data processing in
// atmospheric sciences.
//
// AtmoData is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// AtmoData is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoData home page:
//      http://cerea.enpc.fr/polyphemus/atmodata.html


#ifndef ATMODATA_FILE_METEOROLOGYSTREET_HXX
//#ifndef METEOROLOGYSTREET_HXX

//#include <list>
#include<cmath>

namespace AtmoData
{

  template<class T>
  T ComputeMacdonaldD(T Hmean, T Lp, T A);

  template<class T>
  T ComputeMacdonaldZ0(T Hmean, T Lf, T Lp, T d, T B, T Cd);

  template<class T>
  T ComputeMacdonaldUstar(T Uref, T Zref, T d, T z0);

  template<class T>
  T ComputeMacdonaldUH(T H, T ustar, T d, T z0);

  template<class T>
  T ComputeSiraneC(T z0s, T delta);

  template<class T>
  T ComputeSiraneUm(T C, T ustar, T z0s);

  template<class T>
  T ComputeUmtoUHFactor(T C, T z0s);

  template<class T>
  T ComputeSiraneUstreet(T H, T W, T C, T z0s, T delta, T Um);

  template<class T>
  T ComputeSiraneUstarProfile(T z, T H, T C, T delta, T Um);

  template<class T>
  T ComputeExpUstreet(T H, T W, T z0s, T UH);

  template<class T>
  T ComputeExpUstarProfile(T z, T H, T W, T ustar_H);

  template<class T>
  T BESSI0(T X);

  template<class T>
  T BESSI1(T X);

  template<class T>
  T BESSK0(T X);

  template<class T>
  T BESSK1(T X);

  template<class T>
  T ComputeWangsH(T H, T W, T L, T hmax, T htrunk, T LAI, int nbtree, T Cdt);

  template<class T>
  T ComputeWangUstreet(T H, T W, T L, T phi, T z0s, T sH, int nz, T UH, T hmax, T htrunk, T LAI, int nbtree, T Cdt);

  template<class T>
  T ComputeWangUstarProfile(T z, T H, T W, T L, T z0s, T sH, T ustar_H, T hmax, T htrunk, T LAI, int nbtree, T Cdt);

  template<class T>
  T ComputeSchulteLm(T H, T W);


} //namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGYSTREET_HXX
#endif
