// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
//            Alice Maison - 02/2022 
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


#ifndef ATMODATA_FILE_DEPOSITION_CXX

#include "Deposition.hxx"

namespace AtmoData
{

  //! Computes the bulk surface resistance according to Wesely (1989) and
  //! its revision in Walmsley and Wesely (1996).
  /*!
    \param surface_temperature surface temperature (K).
    \param solar_radiation solar radiation (W.m^{-2}).
    \param species species name, useful only if species=="O3"
    or species=="SO2".
    \param reactivity reactivity factor for the species.
    \param diffusivity molecular diffusivity for the species in air.
    \param Henry Henry constant for the species.
    \param Ri minimum bulk canopy stomatal resistance for water vapour.
    \param Rlu resistance of the outer surfaces in the upper canopy.
    \param Rac in-canopy aerodynamic resistance.
    \param RgsS soil resistance for SO2.
    \param RgsO soil resistance for O3.
    \param RclS resistance due to exposed surface in the lower canopy,
    for SO2.
    \param RclO resistance due to exposed surface in the lower canopy, for O3.
    \param limit (optional) upper limit for resistances. Default: 1.e10.
    \param D_H2O (optional) molecular diffusivity for water vapor.
    Default: 0.25.
    \return the bulk surface resistance.
  */
  template <class T>
  inline T ComputeWesely(T surface_temperature, T solar_radiation,
                         string species, T reactivity, T diffusivity, T Henry,
                         T Ri, T Rlu, T Rac, T RgsS, T RgsO, T RclS, T RclO,
                         T limit, T D_H2O)
  {

    T Rsx, Rmx, Rlux, Rdc, Rclx, Rgx;

    surface_temperature -= 273.15;

    // Stomatal resistance.
    if ((surface_temperature < 0.01) || (surface_temperature > 39.99))
      Rsx = 1.e30;
    else
      Rsx = Ri * (1. + 400. / ((solar_radiation + 0.1)
                               * (solar_radiation + 0.1)))
        * (400. / (surface_temperature * (40. - surface_temperature)));
    Rsx *= D_H2O / diffusivity;

    // Mesophyll resistance.
    Rmx = 1. / (Henry / 3000. + 100. * reactivity);

    // Cuticle resistance.
    Rlux = Rlu / (1.e-5 * Henry + reactivity);

    // Buoyant-convection resistance.
    Rdc = 100. * (1. + 1000. / (solar_radiation + 10.));

    // Rclx.
    if (species == "O3")
      Rclx = RclO;
    else if (species == "SO2")
      Rclx = RclS;
    else
      Rclx = 1. / (1.e-5 * Henry / RclS + reactivity / RclO);

    // Rgx.
    if (species == "O3")
      Rgx = RgsO;
    else if (species == "SO2")
      Rgx = RgsS;
    else
      Rgx = 1. / (1.e-5 * Henry / RgsS + reactivity / RgsO);

    // Limit.
    if (Rsx > limit) Rsx = limit;
    if (Rmx > limit) Rmx = limit;
    if (Rlux > limit) Rlux = limit;
    if (Rdc > limit) Rdc = limit;
    if (Rclx > limit) Rclx = limit;
    if (Rgx > limit) Rgx = limit;

    // Rc.
    return  1. / (1. / (Rsx + Rmx) + 1. / Rlux
                  + 1. / (Rdc + Rclx) + 1. / (Rac + Rgx));

  }


  //! Compute quasi-laminar boundary layer resistance for gaseous species
  /*!
    \param ustar surface friction velocity (m/s)
    \param Sc Schmidt number (-)
    \return Rb quasi-laminar boundary layer resistance (s/m)
  */
  template<class T>
  T ComputeBoundaryLayerResistance(T ustar, T Sc)
  {
    const T karman = 0.41; // Von Karman constant (-)
    const T Pr = 0.74; // Prandtl number (-)
    T p = 2. / 3.; // Baldocchi et al. 1987
    T Rb = 1. / (karman * ustar) * pow(Sc / Pr, p);
    return Rb;
  }

  //! Compute ground and walls surface resistances for gaseous species
  /*!
    \param RgSO2 surface resistance for SO2 (s/m)
    \param RgO3 surface resistance for O3 (s/m)
    \alpha species-dependent scaling factor (-)
    \beta species-dependent scaling factor (-)
    \param Tair air temperature (°C)
    \return Rg surface resistance for gaseous species (s/m)
  */
  template<class T>
  T ComputeGroundResistance(T RgSO2, T RgO3, T alpha, T beta, T Tair)
  {
    T Rg;
    if (alpha == 0. && beta == 0.)
      Rg = 1e30;
    else
      Rg = 1. / ((alpha / RgSO2) + (beta / RgO3));
    if (Tair < -1.)
      Rg =  Rg * exp(-0.2 * (1. + Tair));
    return Rg;
  }

  //! Gas deposition on tree leaves
  //! Compute mesophyll resistance for gaseous species
  /*!
    \param Henry Henry constant for the species (mol/(L*atm))
    \param f0 reactivity factor for the species (-)
    \return Rmes mesophyll resistance for gaseous species (s/m)
  */
  template<class T>
  T ComputeMesophyllResistance(T Henry, T f0)
  {
    T Rmes = pow(Henry / 3000. + 100. * f0, -1.);
    return Rmes;
  }

  //! Compute cuticular resistance for gaseous species
  //! on a dry leaf surface (Zhang et al., 2003)
  /*!
    \param RcutSO2 cuticular resistance for SO2 (s/m)
    \param RcutO3 cuticular resistance for O3 (s/m)
    \alpha species-dependent scaling factor (-)
    \beta species-dependent scaling factor (-)
    \param Tair air temperature (°C)
    \param RH relative humidity (-)
    \param LAI tree Leaf Area Index (-)
    \return Rcut cuticular resistance for gaseous species (s/m)
  */
  template<class T>
  T ComputeZhangCuticularResistance(T RcutSO2_d0, T RcutO3_d0, T alpha, T beta, T Tair, T RH, T ustar, T LAI)
  {
    T RcutSO2 = RcutSO2_d0 / (exp(0.03 * RH) * pow(LAI, 1/4) * ustar);
    T RcutO3 = RcutO3_d0 / (exp(0.03 * RH) * pow(LAI, 1/4) * ustar);
    T Rcut;
    if (alpha == 0. && beta == 0.)
      Rcut = 1e30;
    else
      Rcut = 1. / ((alpha / RcutSO2) + (beta / RcutO3));
    if (Tair < -1.)
      Rcut =  Rcut * exp(-0.2 * (1. + Tair));
    return Rcut;
  }

  //! Compute cuticular resistance for gaseous species
  //! on a dry leaf surface (Wesely, 1989; Walmsley and Wesely, 1996)
  /*!
    \param H Henry constant for the species (mol/(L*atm))
    \param f0 reactivity factor for the species (-)
    \param Tair air temperature (°C)
    \param LAI tree Leaf Area Index (-)
    \return Rcut cuticular resistance for gaseous species (s/m)
  */
  template<class T>
  T ComputeWeselyCuticularResistance(T H, T f0, T Tair, T LAI)
  {
    T Rcut_ = 6000. - 4000. * tanh(1.6 * (LAI - 1.6)); // formula used in SURFEX to replace Wesely table
    T Rcut = Rcut_ * pow(H * 1e-5 + f0, -1.);
    if (Tair < -2.)
      Rcut = Rcut * 1000. * exp(-(Tair + 4.));
    return Rcut;
  }

  //! Compute stomatal resistance for gaseous species
  //! on a dry leaf surface (Wesely, 1989; Walmsley and Wesely, 1996)
  /*!
    \param Dm molecular diffusivity for the species (cm2/s)
    \param DmH2O molecular diffusivity for water vapor (cm2/s)
    \param Rsmin minimum stomatal resistance for water vapor (s/m)
    \param G incoming solar radiation (W/m2)
    \param Ts leaf surface temperature (°C)
    \return Rsto stomatal resistance for gaseous species (s/m)
  */
  template<class T>
  T ComputeWeselyStomatalResistance(T Dm, T DmH2O, T Rsmin, T G, T Ts)
  {
    T RstoH2O = Rsmin * (1 + pow(200. / (G + 0.1), 2.)) * (400. / (Ts * (40. - Ts)));
    T Rsto = RstoH2O * DmH2O / Dm;
    if (Ts < 0. or Ts > 40.)
      Rsto = 1e30;
    return Rsto;
  }

  //! Particles deposition
  //! Compute mean free path of air molecules
  /*!
    \param Tair air temperature (K)
    \param Pair air pressure (Pa)
    \param Rair air specific constant (J/(K.kg))
    \param dyn_vis dynamic viscosity of air (kg/(m.s))
    \return L mean free path of air molecules (m)
  */
  template<class T>
  T ComputeMeanFreePath(T Tair, T Pair, T dyn_vis)
  {
    const T pi = acos(-1);
    const T Rair = 287.058; // air specific constant (J/K/kg)
    T L = 2. * dyn_vis / Pair * pow(8. / (pi * Rair * Tair), -0.5);
    return L;
  }

  //! Compute Cunningham number
  /*!
    \param dp particle diameter (m)
    \param L mean free path of air molecules (m)
    \return Cu Cunningham number (-)
  */
  template<class T>
  T ComputeCunninghamNumber(T dp, T L)
  {
    // Cunningham factor depends on 3 coefficienits (Knudsen and Weber, 1911),
    // which are experimentally determined (Davies, 1945):
    const T a1 = 1.257;
    const T a2 = 0.400;
    const T a3 = 0.550;
    T Cu = 1. + 2. * L / dp * (a1 + a2 * exp(-a3 * dp / L));
    return Cu;
  }

  //! Compute particle relaxation time
  /*!
    \param dp particle diameter (m)
    \param Cu Cunningham number (-)
    \param dyn_vis dynamic viscosity of air (kg/(m.s))
    \param rho_p particle density (kg/m3)
    \return tau particle relaxation time (s)
  */
  template<class T>
  T ComputeRelaxationTime(T dp, T Cu, T dyn_vis, T rho_p)
  {
    T tau = rho_p * pow(dp, 2.) * Cu / (18. * dyn_vis);
    return tau;
  }

  //! Compute Reynolds number Re*
  /*!
    \param ustar surface friction velocity (m/s)
    \param kin_vis kinematic viscosity of air (m2/s)
    \param z0s surface (street ground or walls) roughness height (m)
    \return Re_star Reynolds number (-)
  */
  template<class T>
  T ComputeReynoldsNumber(T ustar, T z0s, T kin_vis)
  {
    T Re_star = ustar * z0s / kin_vis;
    return Re_star;
  }

  //! Compute aerosol molecular diffusivity in air
  /*!
    \param dp particle diameter (m)
    \param Cu Cunningham number (-)
    \param Tair air temperature (K)
    \param dyn_vis dynamic viscosity of air (kg/(m.s))
    \return Dm molecular diffusivity of particle in air (m2/s)
  */
  template<class T>
  T ComputeMolecularDiffusivity(T dp, T Cu, T Tair, T dyn_vis)
  {
    const T pi = acos(-1);
    const T Kb = 1.3806503e-23; // Boltzmann constant (in J/K or m^2.kg/s^2/K)
    T Dm = Kb * Tair * Cu / (3. * pi * dp * dyn_vis);
    return Dm;
  }

  //! Compute Schmidt number
  /*!
    \param kin_vis kinematic viscosity of air (m2/s)
    \param Dm molecular diffusivity of particle in air (m2/s)
    \return Schmidt number (-)
  */
  template<class T>
  T ComputeSchmidtNumber(T kin_vis, T Dm)
  {
    T Sc = kin_vis / Dm;
    return Sc;
  }

  //! Compute sedimentation velocity
  /*!
    \param dp particle diameter (m)
    \param Cu Cunningham number (-)
    \param dyn_vis dynamic viscosity of air (kg/(m.s))
    \param rho_p particle density (kg/m3)
    \return Vs sedimentation velocity (m/s)
  */
  template<class T>
  T ComputeSedimentationVelocity(T dp, T Cu, T dyn_vis, T rho_p)
  {
    const T g = 9.81; // Earth gravitational acceleration (m/s2)
    T Vs = rho_p * pow(dp, 2.) * g * Cu / (18. * dyn_vis);
    return Vs;
  }

  //! Compute Stokes number
  /*!
    \param ustar surface friction velocity (m/s)
    \param df characteristic collector radius (m)
    \param Vs sedimentation velocity (m/s)
    \param kin_vis kinematic viscosity of air (m2/s)
    \return St Stokes number (-)
  */
  template<class T>
  T ComputeStokesNumber(T ustar, T df, T Vs, T kin_vis)
  {
    const T g = 9.81; // Earth gravitational acceleration (m/s2)
    T St;
    if (df == 0.) // smooth surface
      St = Vs * pow(ustar, 2.) / (g * kin_vis);
    else // rough surface
      St = Vs * ustar / (g * df);
    return St;
  }

  //! Compute impaction efficiency coefficient (Zhang et al., 2001)
  /*!
    \param St Stokes number (-)
    \param alpha constant of Zhang et al. (2001) model which depends on land-use (-)
    \param beta constant of Zhang et al. (2001) model which depends on land-use (-)
    \return Eim impaction efficiency coefficient (-)
  */
  template<class T>
  T ComputeZhangImpactionEfficiency(T St, T alpha, T beta)
  {
    T Eim = pow(St / (alpha + St), beta);
    return Eim;
  }

  //! Compute particle rebound coefficient
  /*!
    \param St Stokes number (-)
    \return Rr rebound coefficient (-)
  */
  template<class T>
  T ComputeReboundCoefficient(T St)
  {
    T Rr = exp(-2. * sqrt(St));
    return Rr;
  }

  //! Zhang et al. (2001) model
  //! Compute interception efficiency coefficient (Zhang et al., 2001)
  /*!
    \param dp particle diameter (m)
    \param df characteristic collector radius (m)
    \return Ein interception efficiency coefficient (-)
  */
  template<class T>
  T ComputeZhangInterceptionEfficiency(T dp, T df)
  {
    T Ein;
    if (df == 0.)
      Ein = 0.;
    else
      Ein = 0.5 * pow(dp / df, 2.);
    return Ein;
  }

  //! Compute Brownian efficiency coefficient (Zhang et al., 2001)
  /*!
    \param Schmidt number (-)
    \param gamma constant of Zhang et al. (2001) model which depends on land-use (-)
    \return Eb Brownian diffusion efficiency coefficient (-)
  */
  template<class T>
  T ComputeZhangBrownianEfficiency(T Sc, T gamma)
  {
    T Eb = pow(Sc, -gamma);
    return Eb;
  }

  //! Compute surface resistance for particles according to Zhang et al. (2001)
  /*!
    \param ustar surface friction velocity (m/s)
    \param Eim impaction efficiency coefficient (-)
    \param Ein interception efficiency coefficient (-)
    \param Eb Brownian diffusion efficiency coefficient (-)
    \param Rr rebound coefficient (-)
    \return Rs surface resistance for particles (s/m)
  */
  template<class T>
  T ComputeZhangSurfaceResistance(T ustar, T Eim, T Ein, T Eb, T Rr)
  {
    const T epsilon0 = 3.; // Zhang et al. (2001) model constant (idependent of land-use)
    T Rs = 1. / (epsilon0 * ustar * (Eim + Ein + Eb) * Rr);
    return Rs;
  }

  //! Giardina and Buffa (2018) model
  //! Compute inertial impact resistance (Giardina and Buffa, 2018)
  /*!
    \param St Stokes number (-) // Equation with df = 0
    \param ustar surface friction velocity (m/s)
    \param df characteristic collector radius (m)
    \return Rii inertial impact resistance (s/m)
  */
  template<class T>
  T ComputeGiardinaInertialImpactResistance(T St, T ustar, T df)
  {
    T Rii;
    if (df == 0.) // smooth surface
      Rii = 1. / (ustar * (pow(St, 2.) / (pow(St, 2.) + 1.)));
    else // rough surface
      Rii = 1. / (ustar * (pow(St, 2.) / (pow(St, 2.) + 400.)));
    return Rii;
  }

  //! Compute turbulent impact resistance (Giardina and Buffa, 2018)
  /*!
    \param ustar surface friction velocity (m/s)
    \param tau particle relaxation time (s)
    \param kin_vis kinematic viscosity of air (m2/s)
    \return Rti turbulent impact resistance (s/m)
  */
  template<class T>
  T ComputeGiardinaTurbulentImpactResistance(T ustar, T tau, T kin_vis)
  {
    const T m = 0.1; // (Giardina and Buffa, 2018) precedent value in the code: m = 0.05
    const T n = 3.; // (Giardina and Buffa, 2018) precedent value in the code: m = 0.75
    T tau_p = tau * pow(ustar, 2.) / kin_vis; // Dimensionless relaxation time (-)
    T Rti = 1. / (ustar * m * pow(tau_p, n));
    return Rti;
  }

  //! Compute Brownian diffusion resistance (Giardina and Buffa, 2018)
  /*!
    \param Schmidt number (-)
    \param ustar surface friction velocity (m/s)
    \return Rdb Brownian diffusion resistance (s/m)
  */
  template<class T>
  T ComputeGiardinaBrownianDiffusionResistance(T Sc, T ustar)
  {
    T Rdb = 1. / (ustar * pow(Sc, -2./3.));
    return Rdb;
  }

  //! Compute surface resistance for particles according to Giardina and Buffa (2018)
  /*!
    \param Rii inertial impact resistance (s/m)
    \param Rti turbulent impact resistance (s/m)
    \param Rdb Brownian diffusion resistance (s/m)
    \param Rr rebound coefficient (-)
    \return Rs surface resistance for particles (s/m)
  */
  template<class T>
  T ComputeGiardinaSurfaceResistance(T Rii, T Rti, T Rdb, T Rr)
  {
    T Rs = 1. / (pow(Rdb, -1.) + pow(Rii * Rr, -1.) + pow((Rii + Rti) * Rr, -1.));
    return Rs;
  }

  //! Additional parametrisations for Browian diffusion
  //! Compute Brownian efficiency coefficient (C. Seigneur)
  /*!
    \param Schmidt number (-)
    \return Eb Brownian diffusion efficiency coefficient (-)
  */
  template<class T>
  T ComputeSeigneurBrownianEfficiency(T Sc)
  {
    const T epsilon0 = 3.; // Zhang et al. (2001) model constant (idependent of land-use)
    T Eb = pow(Sc, -2. / 3.) / (5. * epsilon0);
    return Eb;
  }

  //! Compute Brownian efficiency coefficient (Chamberlain, 1967)
  /*!
    \param Schmidt number (-)
    \param ustar surface friction velocity (m/s)

    \return Eb Brownian diffusion efficiency coefficient (-)
  */
  template<class T>
  T ComputeChamberlainBrownianEfficiency(T Sc, T ustar, T Re_star)
  {
    T Eb = 1. / (ustar * pow(Sc, -0.5) * pow(Re_star, -0.05));
    return Eb;
  }

  //! Compute particles deposition velocity (Venkatran and Pleim, 1999)
  /*!
    \param Vs sedimentation velocity (m/s)
    \param Rs surface resistance for gas or particles (s/m)
    \return Vd deposition velocity (m/s)
  */
  template<class T>
  T ComputeVenkatranDepositionVelocity(T Vs, T Rs)
  {
    T Vd = Vs / (1. - exp(-Vs * Rs));
    return Vd;
  }


}  // namespace AtmoData.

#define ATMODATA_FILE_DEPOSITION_CXX
#endif
