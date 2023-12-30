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


#ifndef ATMODATA_FILE_METEOROLOGYSTREET_CXX

#include "MeteorologyStreet.hxx"

namespace AtmoData
{

  // Compute transport functions for horizontal and vertical transfers
  // and friction velocity above the road and ground surfaces for deposition
  /*!
     2 parametrizations are available for wind velocity at roof level : Macdonald and Sirane
     3 parametrizations are available for horizontal wind velocity : Wang, Sirane and Exponential
     3 parametrizations are available for vertical transfer velocity : Wang, Sirane and Schulte
  */

  //! Compute friction velocity and wind speed at roof level according to Macdonald et al. (1998)
  /*!
    Compute city displacement height (Macdonald et al., 1998)
    \param Hmean average building height over the city (m)
    \param Lp plan area density (-)
    \param A empiric constant (-) (Macdonald et al., 1998)
    \return d displacement height of the urban canopy (m)
  */
  template<class T>
  T ComputeMacdonaldD(T Hmean, T Lp, T A)
  {
    T d = Hmean * ( 1. + pow(A, -Lp) * (Lp - 1.));
    return d;
  }  

  //! Compute friction velocity and wind speed at roof level according to Macdonald et al. (1998)
  /*!
    Compute city roughness height (Macdonald et al., 1998)
    \param Hmean average building height over the city (m)
    \param Lf frontal area density (-)
    \param Lp plan area density (-)
    \param d displacement height of the urban canopy (m)
    \param B empiric constant (-) (Macdonald et al., 1998)
    \pram Cd building drag coefficient (-)
    \return z0 roughness height of the urban canopy (m)
  */
  template<class T>
  T ComputeMacdonaldZ0(T Hmean, T Lf, T Lp, T d, T B, T Cd)
  {
    const T karman = 0.41; // Von Karman constant (-)
    T tmp = 1. - (d / Hmean);
    T z0 = Hmean * tmp * exp(-pow(0.5 * B * Cd / pow(karman, 2.) * tmp * Lf, -0.5));
    return z0;
  }  

  //! Compute friction velocity and wind speed at roof level according to Macdonald et al. (1998)
  /*!
    Compute friction velocity above the city from the wind speed at a reference height (Macdonald et al., 1998)
    \param Uref Wind speed at the reference height Zref (m/s)
    \param d displacement height of the urban canopy (m)
    \param z0 roughness height of the urban canopy (m)
    \return ustar friction velocity above the urban canopy (m/s)
  */
  template<class T>
  T ComputeMacdonaldUstar(T Uref, T Zref, T d, T z0)
  {
    const T karman = 0.41; // Von Karman constant (-)
    T ustar = Uref * karman / log((Zref - d) / z0);
    return ustar;
  }

  //! Compute friction velocity and wind speed at roof level according to Macdonald et al. (1998)
  /*!
    Compute wind speed at roof level
    \param H building height (m)
    \param ustar friction velocity above the urban canopy (m/s)
    \param d displacement height of the urban canopy (m)
    \param z0 roughness height of the urban canopy (m)
    \return UH wind speed at z = H (m/s)
  */
  template<class T>
  T ComputeMacdonaldUH(T H, T ustar, T d, T z0)
  {
    const T karman = 0.41; // Von Karman constant (-)
    T UH = ustar / karman * log((H - d) / z0);
    return UH;
  }

  //! Compute average wind speed in the street according to Souhlac et al. (2008, 2011)
  /*! 
    Compute C coefficient (Soulhac et al., 2008)
    \param z0s surface (street ground or walls) roughness height (m)
    \param delta = min(H, W/2) boundary-layer thickness (m)
    \return solutionC constant function of street characteristics (-)
  */
  template<class T>
  T ComputeSiraneC(T z0s, T delta)
  {
    const T pi = acos(-1);
    T gamma = 0.577; // Euler constant
    int nc = 100;
    T maxC = 1.;    
    T step = maxC / nc;
    Array<T, 1> temp(nc);
    Array<T, 1> listC(nc);
    T tempC = 0.0;
    T solutionC;

    if ((z0s / delta) > 0.01)
      throw string("Please reduce the surface roughness length, z0_surface (current value: ") + to_str(z0s) + ").";
    
    for (int i = 0; i < nc; ++i)
      {
        tempC += step;
        listC(i) = tempC;
        T y1C = y1(tempC);
        T j1C = j1(tempC);
        temp(i) = abs(2. / tempC * exp(pi / 2. * y1C / j1C - gamma) - z0s / delta);
      }
    T minValue = temp(0);
    int index = 0;
    for (int i = 1; i < nc; ++i)
      {
        if (minValue > temp(i))
          {
            minValue = temp(i);
            index = i;
          }
      }
    solutionC = listC(index);
    T threshold = 0.001;
    if (minValue > threshold)
      throw string("Fail to find a solution. Please adjust maxC (current value: ")
        + to_str(maxC) + ").";
    return solutionC;
  }

  //! Compute average wind speed in the street according to Souhlac et al. (2008, 2011)
  /*! 
    Compute wind speed at roof level in the middle of the street (Soulhac et al., 2008)
    \param C constant function of street characteristics (-)
    \param ustar friction velocity above the urban canopy (m/s)
    \param z0s surface (street ground or walls) roughness height (m)
    \return Um wind speed at roof level in the middle of the street (z = H and y = delta) (m/s)
  */
  template<class T>
  T ComputeSiraneUm(T C, T ustar, T z0s)
  {
    const T pi = acos(-1);
    const T karman = 0.41; // Von Karman constant (-)
    T term1 = pi / (sqrt(2.) * pow(karman, 2.) * C);
    T term2 = y0(C) - (j0(C) * y1(C)) / j1(C);
    T Um = ustar * sqrt(term1 * term2);
    return Um;
  }

  //! Compute average wind speed in the street according to Souhlac et al. (2008, 2011)
  /*! 
  Compute the ratio UH/Um by averaging the horizontal variation of the wind speed (function f(y) in Soulhac et al., 2008)
    \param C constant function of street characteristics (-)
    \param z0s surface (street ground or walls) roughness height (m)
    \return f_mean = UH/Um wind speed at roof level ratio (-)
  */
  template<class T>
  T ComputeUmtoUHFactor(T C, T z0s)
  {
    T y, f_y;
    int ny = 100;
    T f_mean = 0.;
    int i = 0;
    for (int k = 1; k <= ny; ++k)
      {
        y = 1. * k / ny;
        if (y >= z0s)
          {
            f_y = (j1(C) * y0(C*y) - j0(C*y) * y1(C)) / (j1(C) * y0(C) - j0(C) * y1(C));
            f_mean += f_y;
            i++;
          }
      }
    f_mean /= i * 1.;
    return f_mean;
  }

  //! Compute average wind speed in the street according to Souhlac et al. (2008, 2011)
  /*! 
    Compute average wind speed in the street (Soulhac et al., 2008)
    \param H building height (m)
    \param W building width (m)
    \param C constant function of street characteristics (-)
    \param z0s surface (street ground or walls) roughness height (m)
    \param delta = min(H, W/2) boundary-layer thickness (m)
    \param Um wind speed at roof level in the middle of the street (z = H and y = delta) (m/s)
    \return Ustreet average wind speed in the street (m/s)
  */
  template<class T>
  T ComputeSiraneUstreet(T H, T W, T C, T z0s, T delta, T Um)
  {
    T alpha = log(delta / z0s);
    T beta = exp(C / sqrt(2.) * (1. - H / delta));
    T tmp1 = pow(delta, 2.) / (H * W);
    T tmp2 = 2. * sqrt(2.) / C * (1. - beta);
    T tmp3 = 1. - pow(C, 2.) / 3. + pow(C, 4.) / 45.;
    T tmp4 = beta * (2. * alpha - 3.) / alpha;
    T tmp5 = (W / delta - 2.) * (alpha - 1.) / alpha;
    T Ustreet = Um * tmp1 * (tmp2 * tmp3 + tmp4 + tmp5);
    return Ustreet;
  }

  //!Compute ustar at a given altitude in the street from Souhlac et al. (2008, 2011) to compute deposition
  /*! 
    Compute ustar at a given altitude to compute deposition (Soulhac et al., 2008)
    \param z altitude (m)
    \param H building height (m)
    \param C constant function of street characteristics (-)
    \param z0s surface (street ground or walls) roughness height (m)
    \param delta = min(H, W/2) boundary-layer thickness (m)
    \param Um wind speed at roof level in the middle of the street (z = H and y = delta) (m/s)
    \return ustar friction velocity at the altitude z in the street (m/s)
  */
  template<class T>
  T ComputeSiraneUstarProfile(T z, T H, T C, T delta, T Um)
  {
    const T pi = acos(-1);
    const T karman = 0.41; // Von Karman constant (-)
    T ustar = Um * 2. / pi * exp(C / sqrt(2.) * ((z - H) / delta)) * karman * j1(C) / (j1(C) * y0(C) - j0(C) * y1(C));  
    return ustar;
  }

  //! Compute exponential wind speed profile (Masson et al., 2000; Lemonsu et al., 2004; Cherin et al., 2014)
  /*!
    Compute average wind speed in the street from an exponential attenuation wind profile
    \param H building height (m)
    \param W building width (m)
    \param z0s surface (street ground or walls) roughness height (m)
    \param UH wind speed at roof level (z = H) (m/s)
    \return Ustreet average wind speed in the street (m/s)
  */
  template<class T>
  T ComputeExpUstreet(T H, T W, T z0s, T UH)
  {
    T a = 0.5 * H / W;
    T Ustreet = UH * 1. / a * (1. - exp(a * ((z0s / H) - 1.)));
    return Ustreet;
  }

  //! Compute exponential ustar profile (Inoue, 1963; Masson et al., 2000; Lemonsu et al., 2004; Cherin et al., 2014)
  /*!
    Compute ustar at a given altitude in the street from an exp. profile to compute deposition
    \param z altitude (m)
    \param H building height (m)
    \param ustar_H friction velocity above the street (m/s)
    \return ustar friction velocity at the altitude z in the street (m/s)
  */
  template<class T>
  T ComputeExpUstarProfile(T z, T H, T W, T ustar_H)
  {
    T ustar;
    T a = 0.5 * H / W;
    if (z > H)
      throw string("The altitude to compute friction velocity for deposition is superior to the building height.");
    else
      ustar = ustar_H * exp(a * ((z / H) - 1.));
    return ustar;
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Preliminary function: compute first type modified Bessel functions of order 0
    \param X variable
    \return I0(X)
  */
  template<class T>
  T BESSI0(T X)
  {
    T Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
    P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067492;
    P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
    Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
    Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
    Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
    if (abs(X) < 3.75) 
    {
      Y=(X/3.75)*(X/3.75);
      return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
    }
    else 
    {
      AX=abs(X);
      Y=3.75/AX;
      BX=exp(AX)/sqrt(AX);
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
      return (AX*BX);
    }
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Preliminary function: compute first type modified Bessel functions of order 1
    \param X variable
    \return I1(X)
  */
  template<class T>
  T BESSI1(T X)
  {
    T Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
    P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
    P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
    Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
    Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
    Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
    if (fabs(X) < 3.75)
    {
      Y=(X/3.75)*(X/3.75);
      return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
    }
    else
    {
      AX=fabs(X);
      Y=3.75/AX;
      BX=exp(AX)/sqrt(AX);
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
      return (AX*BX);    
    }
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Preliminary function: compute second type modified Bessel functions of order 0
    \param X variable
    \return K0(X)
  */
  template<class T>
  T BESSK0(T X)
  {
    T Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,tmp;
    P1=-0.57721566; P2= 0.42278420; P3=0.23069756; P4=0.3488590e-1;
    P5= 0.262698e-2; P6=0.10750e-3; P7=0.74e-5;
    Q1= 1.25331414; Q2=-0.7832358e-1; Q3=0.2189568e-1; Q4=-0.1062446e-1;
    Q5= 0.587872e-2; Q6=-0.251540e-2; Q7=0.53208e-3;
    if (X == 0.0) return 1e30;  //arbitrary big value
    if (X <= 2.0)
    {
      Y=X*X/4.0;
      AX=-log(X/2.0)*BESSI0(X);
      T tmp = AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      return tmp;
    }
    else
    {
      Y=2.0/X;
      AX=exp(-X)/sqrt(X);
      T tmp = AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))));
      return tmp;
    }
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Preliminary function: compute second type modified Bessel functions of order 1
    \param X variable
    \return K1(X)
  */
  template<class T>
  T BESSK1(T X)
  {
    T Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7;
    P1=1.0; P2=0.15443144; P3=-0.67278579; P4=-0.18156897;
    P5=-0.1919402e-1; P6=-0.110404e-2; P7=-0.4686e-4;
    Q1=1.25331414; Q2=0.23498619; Q3=-0.3655620e-1; Q4=0.1504268e-1;
    Q5=-0.780353e-2; Q6=0.325614e-2; Q7=-0.68245e-3;
    if (X == 0.0) return 1e30;  //arbitrary big value
    if (X <= 2.0)
    {
      Y=X*X/4.0;
      AX=log(X/2.0)*BESSI1(X);
      return (AX+(1.0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
    }
    else
    {
      Y=2.0/X;
      AX=exp(-X)/sqrt(X);
      return (AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7)))))));
    }
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Compute sH factor representing the effect of the canopy on the mixing length in the street (Wang 2012, 2014)
    \param H building height (m)
    \param W street width (m)
    \param L street length (m)
    //! Arguments for trees:
    \param hmax tree height (m)
    \param htrunk tree trunk height (m)
    \param LAI tree Leaf Area Index (-)
    \param nbtree number of trees (-)
    \param Cdt tree drag coefficient (-)
    \return sH factor (-)
  */
  template<class T>
  T ComputeWangsH(T H, T W, T L, T hmax, T htrunk, T LAI, int nbtree, T Cdt)
  {
    const T karman = 0.41; // Von Karman constant (-)
    T sH;
    T lcb = 0.5 * W;
    if (LAI == 0.0 or hmax == 0.0)
      sH = lcb / (karman * H + lcb);
    else
      {
        const T pi = acos(-1);
        T r = (hmax - htrunk)/2.; // tree radius (m)
        T Sleaf = LAI * pi * pow(r,2.); // leaf surface for one tree (m2)
        T LAIstreet = Sleaf * nbtree / (W*L);
        const T Et = 0.054;
        const T a0 = 3.26;
        const T a1 = 0.0256;
        const T a2 = 6.70;
        T lct = Et * H / (Cdt * 0.5 * LAIstreet);
        T fbt = (a0 + a1 * exp(a2 * H / W))/(pow(hmax / H, 2.));
        sH = (lcb * lct * fbt) / ((karman * H * (lcb + lct * fbt)) + (lcb * lct * fbt));
      }
    return sH;
  }

  //! Compute wind speed in the street according to Wang (2012, 2014)
  /*!
    Compute average wind speed in the street from Wang (2012, 2014) wind profile
    \param H building height (m)
    \param W street width (m)
    \param L street length (m)
    \param phi angle between the wind direction and the street orientation (rad)
    \param z0s surface (street ground or walls) roughness height (m)
    \param sH factor representing the effect of the canopy on the mixing length in the street (-)
    \param nz number of points used to integrate the vertical wind profile (-)
    \param UH wind speed at roof level (z = H) (m/s)
    //! Arguments for trees:
    \param hmax tree height (m)
    \param htrunk tree trunk height (m)
    \param LAI tree Leaf Area Index (-)
    \param nbtree number of trees (-)
    \param Cdt tree drag coefficient (-)
    \return Ustreet average wind speed in the street (m/s)
  */
  template<class T>
  T ComputeWangUstreet(T H, T W, T L, T phi, T z0s, T sH, int nz, T UH, T hmax, T htrunk, T LAI, int nbtree, T Cdt)
  {
    const T pi = acos(-1);
    const T karman = 0.41; // Von Karman constant (-)
    T f_phi = 0.;
    if (((phi >= 0.) and (phi <= pi/4)) or ((phi >= 3*pi/4) and (phi <= 5*pi/4)) or ((phi >= 7*pi/4) and (phi <= 2*pi)))
      f_phi = pow(abs(cos(2 * phi)),3.);
    else if (((phi > pi/4) and (phi < 3*pi/4)) or ((phi > 5*pi/4) and (phi < 7*pi/4)))
      f_phi = pow(abs(cos(2 * pi / 4)),3);
    else
      throw string("wind angle is out of range ");
    T Cb = 0.31 * (1. - exp(-1.6 * (H / W))) * f_phi;
    T alpha;
    if (LAI == 0.0)
      alpha = Cb * (H / W)/(karman * sH);
    else
      {
        const T Cu = 6.7;
        T r = (hmax - htrunk)/2.; // tree radius (m)
        T Sleaf = LAI * pi * pow(r,2.); // leaf surface for one tree (m2)
        T LAIstreet = Sleaf * nbtree / (W*L);
        alpha = (Cb * (H / W) + Cdt * Cu * 0.5 * LAIstreet)/(karman * sH);
      }
    T g_z0 = 2. * sqrt(alpha * z0s / H); // g(z0s)
    T g_H = 2. * sqrt(alpha * H / H); // g(H)
    T C1 = 1. / (BESSI0(g_H) - BESSI0(g_z0) * BESSK0(g_H) / BESSK0(g_z0));
    T C2 = -C1 * BESSI0(g_z0) / BESSK0(g_z0);
    T Ustreet = 0.;
    Array<T, 1> z, u_z, function_g_z;
    z.resize(nz + 1);
    u_z.resize(nz + 1);
    function_g_z.resize(nz + 1);
    for (int k = 0; k <= nz; ++k)
    {
      z(k) = z0s + k * (H - z0s) / nz;
      function_g_z(k) = 2. * sqrt(alpha * z(k) / H);
      u_z(k) = UH * (C1 * BESSI0(function_g_z(k)) + C2 * BESSK0(function_g_z(k)));
      Ustreet += u_z(k);
    }
    Ustreet /= (nz + 1);
    return Ustreet;
  }

  //! Compute ustar profile according to Wang (2012, 2014)
  /*!
    Compute ustar at a given altitude in the street from Wang (2012, 2014) to compute deposition
    \param z altitude (m)
    \param H building height (m)
    \param W street width (m)
    \param L street length (m)
    \param z0s surface (street ground or walls) roughness height (m)
    \param sH factor representing the effect of the canopy on the mixing length in the street (-)
    \param ustar_H friction velocity above the street (m/s)
    // Optional arguments for trees:
    \param hmax tree height (m)
    \param htrunk tree trunk height (m)
    \param LAI tree Leaf Area Index (-)
    \param nbtree number of trees (-)
    \param Cdt tree drag coefficient (-)
    \return ustar friction velocity at the altitude z in the street (m/s)
  */
  template<class T>
  T ComputeWangUstarProfile(T z, T H, T W, T L, T z0s, T sH, T ustar_H, T hmax, T htrunk, T LAI, int nbtree, T Cdt)
  {
    const T karman = 0.41; // Von Karman constant (-)
    T f_phi = 1.; // wind angle has no effect on ustar
    T Cb = 0.31 * (1. - exp(-1.6 * (H / W))) * f_phi;
    T alpha;
    if (LAI == 0.0)   // Tree option is deactivated
      alpha = Cb * (H / W)/(karman * sH);
    else             // Tree option is activated
      {
        const T Cu = 6.7;
        const T pi = acos(-1);
        T r = (hmax - htrunk)/2.; // tree radius (m)
        T Sleaf = LAI * pi * pow(r,2.); // leaf surface for one tree (m2)
        T LAIstreet = Sleaf * nbtree / (W*L);
        alpha = (Cb * (H / W) + Cdt * Cu * 0.5 * LAIstreet)/(karman * sH);
      }
    T g_z0 = 2. * sqrt(alpha * z0s / H);
    T g_H = 2. * sqrt(alpha * H / H);
    T C1 = (2. * ustar_H / (karman * sH * g_H)) / (BESSI1(g_H) + BESSI0(g_z0) * BESSK1(g_H) / BESSK0(g_z0)); // Wang (2014) profile version
    T C2 = -C1 * BESSI0(g_z0) / BESSK0(g_z0);
    T gz = 2. * sqrt(alpha * z / H);
    T dU_dz = gz / (2. * z) * (C1 * BESSI1(gz) - C2 * BESSK1(gz));
    T Kz = ustar_H * karman * z * sH;
    T ustar = sqrt(Kz * dU_dz);
    return ustar;
  }

  //! Compute the mixing length used in vertical transfer velocity calculation according to Schulte et al. (2015)
  /*!
    Compute the mixing length in the street according to Schulte et al. (2015)
    \param H building height (m)
    \param W building width (m)
    \return lm mixing length in the street (m)
  */
  template<class T>
  T ComputeSchulteLm(T H, T W)
  {
    const T pi = acos(-1);
    T beta = 2. / (sqrt(2.) * pi);
    T lm = beta / (1 + H / W);
    return lm;
  }


} // namespace AtmoData.

#define ATMODATA_FILE_METEOROLOGYSTREET_CXX
#endif

