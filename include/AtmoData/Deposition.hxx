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


#ifndef ATMODATA_FILE_DEPOSITION_HXX

namespace AtmoData
{

  template <class T>
  T ComputeWesely(T surface_temperature, T solar_radiation,
                  string species, T reactivity, T diffusivity, T Henry,
                  T Ri, T Rlu, T Rac, T RgsS, T RgsO, T RclS, T RclO,
                  T limit = 1.e10, T D_H2O = 0.25);

  template<class T>
  T ComputeBoundaryLayerResistance(T ustar, T Sc);

  template<class T>
  T ComputeGroundResistance(T RgSO2, T RgO3, T alpha, T beta, T Tair);

  template<class T>
  T ComputeMesophyllResistance(T H, T f0);

  template<class T>
  T ComputeZhangCuticularResistance(T RcutSO2_d0, T RcutO3_d0, T alpha, T beta, T Tair, T RH, T ustar, T LAI);

  template<class T>
  T ComputeWeselyCuticularResistance(T H, T f0, T Tair, T LAI);

  template<class T>
  T ComputeWeselyStomatalResistance(T Dm, T DmH2O, T Rsmin, T G, T Ts);

  template<class T>
  T ComputeMeanFreePath(T Tair, T Pair, T dyn_vis);

  template<class T>
  T ComputeCunninghamNumber(T dp, T L);

  template<class T>
  T ComputeRelaxationTime(T dp, T Cu, T dyn_vis, T rho_p);

  template<class T>
  T ComputeReynoldsNumber(T ustar, T z0s, T kin_vis);

  template <class T>
  T ComputeMolecularDiffusivity(T dp, T Cu, T Tair, T dyn_vis);

  template<class T>
  T ComputeSchmidtNumber(T kin_vis, T Dm);

  template<class T>
  T ComputeSedimentationVelocity(T dp, T Cu, T dyn_vis, T rho_p);

  template<class T>
  T ComputeStokesNumber(T ustar, T df, T Vs, T kin_vis);

  template<class T>
  T ComputeZhangImpactionEfficiency(T St, T alpha, T beta);

  template<class T>
  T ComputeZhangInterceptionEfficiency(T dp, T df);

  template<class T>
  T ComputeZhangBrownianEfficiency(T Sc);

  template<class T>
  T ComputeReboundCoefficient(T St);

  template<class T>
  T ComputeZhangSurfaceResistance(T ustar, T Eim, T Ein, T Eb, T Rr);

  template<class T>
  T ComputeGiardinaInertialImpactResistance(T St, T ustar);

  template<class T>
  T ComputeGiardinaTurbulentImpactResistance(T ustar, T tau, T kin_vis);

  template<class T>
  T ComputeGiardinaBrownianDiffusionResistance(T Sc, T ustar);

  template<class T>
  T ComputeGiardinaSurfaceResistance(T Rii, T Rti, T Rdb, T Rr);

  template<class T>
  T ComputeSeigneurBrownianEfficiency(T Sc);

  template<class T>
  T ComputeChamberlainBrownianEfficiency(T Sc, T ustar, T Re_star);

  template<class T>
  T ComputeVenkatranDepositionVelocity(T Vs, T Rs);

}  // namespace AtmoData.

#define ATMODATA_FILE_DEPOSITION_HXX
#endif
