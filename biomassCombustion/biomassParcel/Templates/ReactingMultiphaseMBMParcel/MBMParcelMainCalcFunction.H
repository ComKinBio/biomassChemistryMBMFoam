/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{

  // set the TInf_ as gas temperature
  TInf_ = td.Tc();

  // initialize the time to zero
  tPar = 0.0;

  // initialize total mass of the particle to zero
  mParSum_ = 0.0;

  // set up mesh particle time step
  scalar t_step = TIME_DT_;//cloud.constProps().deltaTime();

  // read the emissivity of the particle
  scalar emissivity = cloud.constProps().epsilon0();

  // Info << "stepFraction is " << this->stepFraction() << nl;
  // Info << "keepParticle is " << td.keepParticle << nl;
  // Info << "switchProcessor is " << td.switchProcessor << nl;



  // * * * * * * * end initialization and start calculation * * * * * * * * * //

  // Info << "code start" << nl;
  tReal_ += dt;
  while (tPar < dt)//t_stop)
  {
    // Info << "tPar is " << tPar << nl;
    // Info << "dt is " << dt << nl;
    /* Initialize matrix system */
    m_zero(A);
    v_zero(b);

    /* Update particle properties including particle volume */
    /* --------------------------------------------------------------------- */
    // update the mass fraction of each layer
    for (label i = 1; i <= 4; i++)
    {
      for (label j = 1; j <= FINE_PARTICLE_MESH_; j++)
      {
        xi_[i-1][j-1] =
          mPar_[i-1][j-1]
          /
          (mPar_[layer_1][j-1]+mPar_[layer_2][j-1]
          +mPar_[layer_3][j-1]+mPar_[layer_4][j-1]);
      }
    }

    // update the volume of each cell
    // due to drying, devolatilization shrinkage
    for (label j = 1; j <= FINE_PARTICLE_MESH_; j++)
    {
      VCell_[j-1] -=
      DRYING_SHRINKAGE_VOLUME_*
      (
        mLast_[layer_1][j-1]-mPar_[layer_1][j-1]
      )/DENSITY_INITIAL_
     +DEVOL_SHRINKAGE_VOLUME_*
      (
        mLast_[layer_2][j-1]-
        (
          mPar_[layer_2][j-1]-
          (mLast_[layer_1][j-1]-mPar_[layer_1][j-1])
          *(1.0-MOISTURE_FRACTION_WB)
        )
      )/density_drywood;

      VChar_[j-1] +=
      (1.0-DEVOL_SHRINKAGE_VOLUME_) *
      (
        mLast_[layer_2][j-1]-
        (
          mPar_[layer_2][j-1]-
          (mLast_[layer_1][j-1]-mPar_[layer_1][j-1])
         *(1.0-MOISTURE_FRACTION_WB)
        )
      )/density_drywood;
    }
    // exit(EXIT_FAILURE);

    if (combustion_layer != 99999999 && combustion_layer != 0)
    {
      VCell_[combustion_layer] -= CHARCONVERSION_SHRINKAGE_VOLUME_
     *Rb_[char_front][combustion_layer]*t_step
     /(mLast_[layer_3][combustion_layer]/VCharLast_[combustion_layer]);

      VChar_[combustion_layer] -=
      Rb_[char_front][combustion_layer]*t_step
     /(mLast_[layer_3][combustion_layer]/VCharLast_[combustion_layer]);
    }

    // updat VCharLast_ and m_char_last
    for (label j = 1; j <= FINE_PARTICLE_MESH_; j++)
    {
      for (label i = 1; i <= 4; i++)
      {
        VCharLast_[j-1] = VChar_[j-1];
        mLast_[i-1][j-1] = mPar_[i-1][j-1];
      }
    }

    volume_sum = 0.5*VCell_[FINE_PARTICLE_MESH_-1];

    for (label j = FINE_PARTICLE_MESH_; j >= 2; j--)
    {
      xc_[j-1] = radius_from_volume(volume_sum);

      volume_sum = volume_sum+0.5*VCell_[j-1]+0.5*VCell_[j-2];
    }
    xc_[0] = radius_from_volume(volume_sum);
    volume_sum = volume_sum+0.5*VCell_[0];
    rb = radius_from_volume(volume_sum);
    // volume_sum = 0.0;

    /* --------------------------------------------------------------------- */
    /* Update reaction rates and heats of reaction based on current Tp_ */
    /* --------------------------------------------------------------------- */
    RProduct_[H2O_index] = 0.0;
    RProduct_[H2_index] = 0.0;
    RProduct_[CO_index] = 0.0;
    RProduct_[CH4_index] = 0.0;
    RProduct_[H2O_index] = 0.0;
    RProduct_[CO2_index] = 0.0;
    R_O2_mole_need = 0.0;

    combustion_flag = 0;
    combustion_layer = 99999999;
    /*Drying, Devolatilization and Char conversion*/
    for (label i = FINE_PARTICLE_MESH_; i >= 1; i--)
    {
      /* Drying */
      // Rb_[drying_front][i-1] = mPar_[layer_1][i-1]*1.0e27*exp(-25000.0/Tp_->ve[i-1]); /* J/s/J/kg = kg/s */
      Rb_[drying_front][i-1] = mPar_[layer_1][i-1] /* J/s/J/kg = kg/s */
       *MOISTURE_FRACTION_WB
       *DryKinA_*Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i-1]));

      if
      (
        mPar_[layer_1][i-1]-Rb_[drying_front][i-1]/MOISTURE_FRACTION_WB*t_step
        <=
        MY_EXTREMELY_SMALL
      )
      {
        Rb_[drying_front][i-1] =
          mPar_[layer_1][i-1]*MOISTURE_FRACTION_WB/t_step;
      }
      QReaction_[evaporation][i-1] = /* (J/kg)*(kg/s) = W */
        deltaHvap(Tp_->ve[i-1])*Rb_[drying_front][i-1];

      QEquilibrate_[drying_front][i-1] =
        Rb_[drying_front][i-1]
       *gas_cp(index_H2O, 0.5*(Tp_->ve[i-1]+TOuterBoundary_))
       *(TOuterBoundary_-Tp_->ve[i-1]); /* (W) */

      RProduct_[H2O_index] = RProduct_[H2O_index]+Rb_[drying_front][i-1];                                                                                         /* (kg/s) */

      // calculate the pyrolysis rate
      // This part of the code now only valid for single component reaction
      RDevolReac_[gassing][i-1] = DevKinA1_
       *Foam::exp(-DevKinE1_/(8.3145*Tp_->ve[i-1]))*mPar_[layer_2][i-1];
      RDevolReac_[tarring][i-1] = DevKinA2_
       *Foam::exp(-DevKinE2_/(8.3145*Tp_->ve[i-1]))*mPar_[layer_2][i-1];
      RDevolReac_[charring][i-1] = DevKinA3_
       *Foam::exp(-DevKinE3_/(8.3145*Tp_->ve[i-1]))*mPar_[layer_2][i-1];

      // calculate the total mass loss rate due to pyrolysis
      Rb_[devol_front][i-1] = /* (kg/s) */
        RDevolReac_[gassing][i-1]
       +RDevolReac_[tarring][i-1]
       +RDevolReac_[charring][i-1];

      // calculate the char production rate due to pyrolysis
      Rb_[char_produced][i-1] = RDevolReac_[charring][i-1];/* (kg/s) */

      /* heat of devolatilization is equal to ZERO */
      QReaction_[devolatilization][i-1] = 0.0*Rb_[devol_front][i-1];

      // referring to Eq.(A.19)
      QEquilibrate_[devol_front][i-1] = /* (W) */
        (RDevolReac_[gassing][i-1]+RDevolReac_[tarring][i-1])
       *cp_g*(TOuterBoundary_-Tp_->ve[i-1]);

      // calculate the produced char at this time step
      charFromDevol_[i-1] += Rb_[char_produced][i-1]*t_step;

      // calculate the total mass loss due to pyrolysis at this time step
      totalDevol_[i-1] += Rb_[devol_front][i-1]*t_step;

      // R_H2_Su_m[i-1] = RProduct_[H2_index];
      R_O2_mole_need_list[i-1] = R_O2_mole_need; /* (kg/s) */
      // RProduct_[H2_index] = RProduct_[H2_index]+Y_H2/Y_G*(RDevolReac_[gassing][i-1]+RDevolReac_[tarring][i-1]);

      // calculate how much mole of O2 needed to consume all the produced gas
      //assume volatile gas is fully consumed by O2
      R_O2_mole_need += /* (mol/s) */
        (Rb_[devol_front][i-1]*C_FRACTION_-RDevolReac_[charring][i-1])/M_C
       + Rb_[devol_front][i-1]*H_FRACTION_/M_H/4.
       - Rb_[devol_front][i-1]*O_FRACTION_/M_O/2.;

      if (combustion_flag == 0)
      {
        //get the char front location
        charFrontRadius_ =
          charFrontLocation(rb, (i-1), VCell_, VChar_, xc_, xi_[layer_3][0]);

        /*make sure only one char layer is converted, char will not start
        converting before dry wood is fully converted (in that layer)*/
        if (i == 1 && xi_[layer_3][i-1] > MY_EXTREMELY_SMALL /*&& xi_[layer_2][i-1]<=0.*/)
        {

          D_g = 2.03e-10*Foam::pow((TInf_+TOuterBoundary_)*0.5, 2.0); /* Thunman et al., Comb. Flame 2001 */

          hmi = D_g/(2.0*rb)*
          (
            2.0+
            1.1*Foam::pow((mu/(rho_g*D_g)), 1.0/3.0)
           *Foam::pow(rho_g*rb*vel_mag/mu, 0.6)
          );

          O2_diff_suface = surface_area(rb)*hmi*conc_O2; //O2 transport to the surface(mol/s)=mPar_^2*mPar_/3*mol/m3

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }

          conc_O2_correct = (1.-alpha)*conc_O2;

          omegaC_[i-1] = 2.0*
            (1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1])))/
            (2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1])));

          //Henrik 2002
          betaR_[1-1] = 1.715*Tp_->ve[i-1]*Foam::exp(-9000.0/Tp_->ve[i-1]);
          //             betaR_[1-1] = 4000.*exp(-80000./(Tp_->ve[i-1]*8.3145)); //Sven 2011
          //                     betaR_[2-1] = 3.42*Tp_->ve[i-1]*exp(-15600.0/Tp_->ve[i-1]); //Henrik 2002
          //                     betaR_[3-1] = betaR_[2-1];
          //                 betaR_[4-1] = betaR_[2-1]/1000.0;

          hmi = D_g/(2.0*rb)*
            (
              2.0+
              1.1*Foam::pow((mu/(rho_g*D_g)), 1.0/3.0)
             *Foam::pow(rho_g*rb*vel_mag/mu, 0.6)
            );

          D_gAB = 2.03e-10
           *Foam::pow((Tp_->ve[i-1]+TOuterBoundary_)*0.5, 2.0);

          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          //             eta = 1.0; /* change this one!?! */

          beta_d = hmi*hmia/(hmi+hmia);

          // Refer to Eq.(3) in the paper
          RCharReac_[1-1][i-1] =
            M_C*omegaC_[i-1]
           *conc_O2_correct
           *(betaR_[1-1]*beta_d/(betaR_[1-1]+beta_d))
           *surface_area(charFrontRadius_);

          // Char conversion reation rate
          /* rate of consumption of char (kg/s) */
          Rb_[char_front][i-1] = RCharReac_[1-1][i-1] /*+ RCharReac_[2-1][i-1]+RCharReac_[3-1][i-1]+RCharReac_[4-1][i-1]*/;


          // Heat of reaction due to char conversion
          QReaction_[char_combustion][i-1] = /* kg/s*J/kg = W */
            heatRatio(Tp_->ve[i-1])*(RCharReac_[1-1][i-1]/M_C)*
            (
              (2.0*(omegaC_[i-1]-1.0)/omegaC_[i-1])
             *LHV[index_CO]*M_CO
             +
              ((2.0-omegaC_[i-1])/omegaC_[i-1])
             *LHV[index_CO2]*M_CO2-LHV[index_C]*M_C
            );

          QEquilibrate_[char_front][i-1] = 0.0;

          if (Rb_[char_front][i-1] != 0.)
          {
            combustion_layer = i-1;
            combustion_flag = 1;
          }
        }
        else if
        (
          i != 1 && xi_[layer_3][i-1] > MY_EXTREMELY_SMALL
                 && xi_[layer_3][i-2] <= MY_EXTREMELY_SMALL
                 && xi_[layer_2][i-2] <= MY_EXTREMELY_SMALL
                 && xi_[layer_1][i-2] <= MY_EXTREMELY_SMALL
        )
        {
          D_g = 2.03e-10*Foam::pow((TInf_+TOuterBoundary_)*0.5, 2.0); /* Thunman et al., Comb. Flame 2001 */

          hmi = D_g/(2.0*rb) *
          (
            2.0+
            1.1*Foam::pow((mu/(rho_g*D_g)), 1.0/3.0)
               *Foam::pow(rho_g*rb*vel_mag/mu, 0.6)
          );

          O2_diff_suface = surface_area(rb)*hmi*conc_O2; //O2 transport to the surface(mol/s)=mPar_^2*mPar_/3*mol/m3

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }
          conc_O2_correct = (1.-alpha)*conc_O2;

          omegaC_[i-1] = 2.0 *
            (1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1])))/
            (2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1])));

          betaR_[1-1] =
            1.715*Tp_->ve[i-1]*Foam::exp(-9000.0/Tp_->ve[i-1]); //Henrik 2002

          hmi = D_g/(2.0*rb) *
          (
            2.0+1.1*Foam::pow((mu/(rho_g*D_g)), 1.0/3.0)
                    *Foam::pow(rho_g*rb*vel_mag/mu, 0.6)
          );

          D_gAB =
            2.03e-10*Foam::pow((Tp_->ve[i-1]+TOuterBoundary_)*0.5, 2.0);

          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          //                     eta = 1.0; /* change this one!?! */

          beta_d = hmi*hmia/(hmi+hmia);

          RCharReac_[1-1][i-1] = /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */ /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */
            M_C*omegaC_[i-1]*conc_O2_correct*
            (betaR_[1-1]*beta_d/(betaR_[1-1]+beta_d))*
            surface_area(charFrontRadius_);

          Rb_[char_front][i-1] = RCharReac_[1-1][i-1];/* rate of consumption of char (kg/s) */
          QReaction_[char_combustion][i-1] = heatRatio(Tp_->ve[i-1])*
            (RCharReac_[1-1][i-1]/M_C)*
            (
              (2.0*(omegaC_[i-1]-1.0)/omegaC_[i-1]) *
              LHV[index_CO]*M_CO
             +
              ((2.0-omegaC_[i-1])/omegaC_[i-1]) *
              LHV[index_CO2]*M_CO2
             -
              LHV[index_C]*M_C
            ); /* kg/s*J/kg = W */

          QEquilibrate_[char_front][i-1] = 0.0;

          if (Rb_[char_front][i-1] != 0.)
          {
            combustion_layer = i-1;
            combustion_flag = 1;
          }
        }
        else
        {
          // inner char has no conversion
          Rb_[char_front][i-1] = 0.0;
          QReaction_[char_combustion][i-1] = 0.0;
          QEquilibrate_[char_front][i-1] = 0.0;
        }
      }
      else
      {
        // inner char has no conversion
        Rb_[char_front][i-1] = 0.0;
        QReaction_[char_combustion][i-1] = 0.0;
        QEquilibrate_[char_front][i-1] = 0.0;
      }
    }

    /* Source term information storage part 1*/
    R_drySourceimins1High = 0.;
    R_drySourceiHigh = 0.;
    R_devoSourceimins1High = 0.;
    R_devoSourceimins1HighChar = 0.;
    R_devoSourceiHigh = 0.;
    R_devoSourceiHighChar = 0.;

    R_drySourceimins1Low = 0.;
    R_drySourceiLow = 0.;
    R_devoSourceimins1Low = 0.;
    R_devoSourceimins1LowChar = 0.;
    R_devoSourceiLow = 0.;
    R_devoSourceiLowChar = 0.;

    for (label i = FINE_PARTICLE_MESH_; i >= 1; i--)
    {
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sStar_~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (i == FINE_PARTICLE_MESH_)
      {
        sStar_[i-1] = 0.0;
      }
      else if (i == 1)
      {
        sStar_[i-1] = QReaction_[devolatilization][i]
                     -QReaction_[evaporation][i]
                     -QReaction_[char_combustion][i]
                     -QEquilibrate_[char_front][i]
                     -QEquilibrate_[devol_front][i]
                     -QEquilibrate_[drying_front][i]
                     +QReaction_[devolatilization][i-1]
                     -QReaction_[evaporation][i-1]
                     -QReaction_[char_combustion][i-1]
                     -QEquilibrate_[char_front][i-1]
                     -QEquilibrate_[devol_front][i-1]
                     -QEquilibrate_[drying_front][i-1];
      }
      else
      {
        sStar_[i-1] = QReaction_[devolatilization][i]
                     -QReaction_[evaporation][i]
                     -QReaction_[char_combustion][i]
                     -QEquilibrate_[char_front][i]
                     -QEquilibrate_[devol_front][i]
                     -QEquilibrate_[drying_front][i];
      }

      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~dSdT_star~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (i == FINE_PARTICLE_MESH_)
      {
        dSdT_star[i-1] = 0.0;
      }
      else if (i == 1)
      {
        //Dry Devo high~~~~~~~~~~~~~~~~~~~~~
        R_drySourceimins1High =  DryKinA_*
          Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i-1]+1.0))*mPar_[layer_1][i-1];

        scalarField dMassRateHigh(3, 0.0);

        dMassRateHigh[gassing] = DevKinA1_
          *Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i-1]+1)))*mPar_[layer_2][i-1];
        dMassRateHigh[tarring] = DevKinA2_
          *Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i-1]+1)))*mPar_[layer_2][i-1];
        dMassRateHigh[charring] = DevKinA3_
          *Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i-1]+1)))*mPar_[layer_2][i-1];

        R_devoSourceimins1High =
          (
            DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i-1]+1.0)))
           +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i-1]+1.0)))
          )*mPar_[layer_2][i-1];

        R_devoSourceimins1HighChar = DevKinA3_
         *Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i-1]+1.0)))
         *mPar_[layer_2][i-1];

        R_drySourceiHigh = DryKinA_
         *Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i]+1.0))
         *mPar_[layer_1][i];

        R_devoSourceiHigh =
          (
            DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i]+1.0)))
           +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i]+1.0)))
          )*mPar_[layer_2][i];

        R_devoSourceiHighChar = DevKinA3_
         *Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i]+1.0)))
         *mPar_[layer_2][i];

        high =
          -
          (
            deltaHvap(Tp_->ve[i-1]+1.0)
           +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i-1]+1.0)+TOuterBoundary_))
            *(TOuterBoundary_-(Tp_->ve[i-1]+1.0))
          )*R_drySourceimins1High
          -
          R_devoSourceimins1High*cp_g*(TOuterBoundary_-(Tp_->ve[i-1]+1.0))
          -
          (
            deltaHvap(Tp_->ve[i]+1.0)
           +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i]+1.0)+TOuterBoundary_))
            *(TOuterBoundary_-(Tp_->ve[i]+1.0))
          )*R_drySourceiHigh
          -
          R_devoSourceiHigh*cp_g*(TOuterBoundary_-(Tp_->ve[i]+1.0));

        //Dry Devo low~~~~~~~~~~~~~~~~~~~~~
        R_drySourceimins1Low = DryKinA_*
          Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i-1]-1.0))*mPar_[layer_1][i-1];

        R_devoSourceimins1Low =
          (
            DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i-1]-1.0)))
           +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i-1]-1.0)))
          )*mPar_[layer_2][i-1];

        R_devoSourceimins1LowChar = DevKinA3_*
          Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i-1]-1.0)))*mPar_[layer_2][i-1];

        R_drySourceiLow = DryKinA_*
          Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i]-1.0))*mPar_[layer_1][i];

        R_devoSourceiLow =
          (
            DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i]-1.0)))
           +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i]-1.0)))
          )*mPar_[layer_2][i];

        R_devoSourceiLowChar =
          DevKinA3_
            *Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i]-1.0)))*mPar_[layer_2][i];

        low =
        -
          (
            deltaHvap(Tp_->ve[i-1]-1.0)
           +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i-1]-1.0)+TOuterBoundary_))
            *(TOuterBoundary_-(Tp_->ve[i-1]-1.0))
          )*R_drySourceimins1Low-R_devoSourceimins1Low*cp_g
          *(TOuterBoundary_-(Tp_->ve[i-1]-1.0))
        -
          (
            deltaHvap(Tp_->ve[i]-1.0)
           +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i]-1.0)+TOuterBoundary_))
            *(TOuterBoundary_-(Tp_->ve[i]-1.0))
          )*R_drySourceiLow-R_devoSourceiLow*cp_g
          *(TOuterBoundary_-(Tp_->ve[i]-1.0));

        //Char conversion~~~~~~~~~~
        if (i == combustion_layer+1) // combustion at i-1, NO combustion at i
        {

          charFrontRadius_ =
            charFrontLocation(rb, (i-1), VCell_, VChar_, xc_, xi_[layer_3][0]);
          //HIGH~~~~~~~~~~~~~~~~~~~~~~~~~~~

          //alpha = ((R_H2_Su_m[i]+(R_devoSourceimins1High+R_devoSourceiHigh)*Y_H2/Y_G)/M_H2/2.)/(O2_diff_suface);
          R_O2_mole_need =
            R_O2_mole_need_list[i]
            +
            (
              (
                R_devoSourceimins1High
               +R_devoSourceiHigh
               +R_devoSourceimins1HighChar
               +R_devoSourceiHighChar
              )*C_FRACTION_
              -
              (R_devoSourceimins1HighChar+R_devoSourceiHighChar)
            )/M_C
            +
            (
              R_devoSourceimins1High
             +R_devoSourceiHigh
             +R_devoSourceimins1HighChar
             +R_devoSourceiHighChar
            )*(H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }

          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*
            Foam::pow((Tp_->ve[i-1]+1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1]+1.0)))
            /(2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1]+1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*(1.715*(Tp_->ve[i-1]+1.0)
            *Foam::exp(-9000.0/(Tp_->ve[i-1]+1.0))*beta_d/(1.715*(Tp_->ve[i-1]+1.0)
            *Foam::exp(-9000.0/(Tp_->ve[i-1]+1.0))+beta_d))
            *surface_area(charFrontRadius_);
          /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */ /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */

          high -=
            heatRatio(Tp_->ve[i-1]+1.0)*
            (
              (RCharReac_1S/M_C)*
              (
                (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
               +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
               -LHV[index_C]*M_C
              )
            );

          R_O2_mole_need = R_O2_mole_need_list[i]+
          (
            (
              R_devoSourceimins1Low
             +R_devoSourceiLow
             +R_devoSourceimins1LowChar
             +R_devoSourceiLowChar
            )*C_FRACTION_
            -(R_devoSourceimins1LowChar+R_devoSourceiLowChar)
          )/M_C
          +
          (
            R_devoSourceimins1Low
           +R_devoSourceiLow
           +R_devoSourceimins1LowChar
           +R_devoSourceiLowChar
          )*(H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }

          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*
            Foam::pow((Tp_->ve[i-1]-1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1]-1.0)))/
            (2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i-1]-1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*
          (
            1.715*(Tp_->ve[i-1]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i-1]-1.0))
            *beta_d/
            (
              1.715*(Tp_->ve[i-1]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i-1]-1.0))+beta_d
            )
          )*surface_area(charFrontRadius_); /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */ /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */

          low -= heatRatio(Tp_->ve[i-1]-1.0)*
            (
              (RCharReac_1S/M_C)*
                (
                  (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
                 +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
                 -LHV[index_C]*M_C
                )
            );
        }
        else if (i == combustion_layer) // no combustion at i-1, combustion at i
        {
          charFrontRadius_ =
            charFrontLocation(rb, (i), VCell_, VChar_, xc_, xi_[layer_3][0]);

          //HIGH~~~~~~~~~~~~~~~~~~~~~~~~~~~
          R_O2_mole_need =
            R_O2_mole_need_list[i]
            +
            (
              (
                R_devoSourceimins1High
               +R_devoSourceiHigh
               +R_devoSourceimins1HighChar
               +R_devoSourceiHighChar
              )*C_FRACTION_
              -
              (
                R_devoSourceimins1HighChar
               +R_devoSourceiHighChar
              )
            )/M_C
            +
            (
              R_devoSourceimins1High
             +R_devoSourceiHigh
             +R_devoSourceimins1HighChar
             +R_devoSourceiHighChar
            )*(H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }

          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*Foam::pow((Tp_->ve[i]+1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]+1.0)))
            /(2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]+1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*
            (
              1.715*(Tp_->ve[i]+1.0)*Foam::exp(-9000.0/(Tp_->ve[i]+1.0))*beta_d/
              (1.715*(Tp_->ve[i]+1.0)*Foam::exp(-9000.0/(Tp_->ve[i]+1.0))+beta_d)
            )*surface_area(charFrontRadius_); /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */ /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */
          ;

          high -= heatRatio(Tp_->ve[i]+1.0)*
            (
              (RCharReac_1S/M_C)*
              (
                (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
               +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
               -LHV[index_C]*M_C
              )
            );

          //LOW~~~~~~~~~~~~~~~~~~~~~~~~~~~
          //                     alpha = ((R_H2_Su_m[i]+(R_devoSourceimins1Low+R_devoSourceiLow)*Y_H2/Y_G)/M_H2/2.)/(O2_diff_suface);
          R_O2_mole_need =
            R_O2_mole_need_list[i]
            +
            (
              (
                R_devoSourceimins1Low
               +R_devoSourceiLow
               +R_devoSourceimins1LowChar
               +R_devoSourceiLowChar
              )*C_FRACTION_
              -
              (R_devoSourceimins1LowChar+R_devoSourceiLowChar)
            )/M_C
            +
            (
              R_devoSourceimins1Low
             +R_devoSourceiLow
             +R_devoSourceimins1LowChar
             +R_devoSourceiLowChar
            )*(H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);

          if (alpha > 1.0)
          {
            alpha = 1.0;
          }

          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*Foam::pow((Tp_->ve[i]-1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]-1.0)))
            /(2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]-1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*
          (
            1.715*(Tp_->ve[i]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i]-1.0))*beta_d/
            (1.715*(Tp_->ve[i]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i]-1.0))+beta_d)
          )*surface_area(charFrontRadius_); /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */

          low -= heatRatio(Tp_->ve[i]-1.0)*
            (
              (RCharReac_1S/M_C)*
              (
                (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
               +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
               -LHV[index_C]*M_C
              )
            );
        }
        else //no combustion
        {
          //do nothing
        }
      }
      else
      {
        //Dry Devo high~~~~~~~~~~~~~~~~~~~~~
        R_drySourceiHigh = DryKinA_*
          Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i]+1.0))*mPar_[layer_1][i];

        R_devoSourceiHigh =
        (
          DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i]+1.0)))
         +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i]+1.0)))
        )*mPar_[layer_2][i];

        R_devoSourceiHighChar = DevKinA3_*
          Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i]+1.0)))*mPar_[layer_2][i];

        high =
          -
          (
            deltaHvap(Tp_->ve[i]+1.0)
           +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i]+1.0)+TOuterBoundary_))
            *(TOuterBoundary_-(Tp_->ve[i]+1.0))
          )*R_drySourceiHigh
          -
          R_devoSourceiHigh*cp_g*(TOuterBoundary_-(Tp_->ve[i]+1.0));

        //Dry Devo low~~~~~~~~~~~~~~~~~~~~~
        R_drySourceiLow = DryKinA_*
          Foam::exp(-DryKinE_/(8.3145*Tp_->ve[i]-1.0))*mPar_[layer_1][i];

        R_devoSourceiLow =
          (
            DevKinA1_*Foam::exp(-DevKinE1_/(8.3145*(Tp_->ve[i]-1.0)))
           +DevKinA2_*Foam::exp(-DevKinE2_/(8.3145*(Tp_->ve[i]-1.0)))
          )*mPar_[layer_2][i];

        R_devoSourceiLowChar = DevKinA3_*
          Foam::exp(-DevKinE3_/(8.3145*(Tp_->ve[i]-1.0)))*mPar_[layer_2][i];

        low =
          -
          (
            deltaHvap(Tp_->ve[i]-1.0)
            +
            gas_cp(index_H2O, 0.5*((Tp_->ve[i]-1.0)+TOuterBoundary_))*
            (TOuterBoundary_-(Tp_->ve[i]-1.0))
          )*R_drySourceiLow
          -
          R_devoSourceiLow*cp_g*(TOuterBoundary_-(Tp_->ve[i]-1.0));

        //Char conversion~~~~~~~~~~
        if (i == combustion_layer) // no combustion at i-1, combustion at i
        {
          charFrontRadius_ =
            charFrontLocation(rb, (i), VCell_, VChar_, xc_, xi_[layer_3][0]);
          //HIGH~~~~~~~~~~~~~~~~~~~~~~~~~~~
          //                     alpha = ((R_H2_Su_m[i]+R_devoSourceiHigh*Y_H2/Y_G)/M_H2/2.)/(O2_diff_suface);
          R_O2_mole_need = R_O2_mole_need_list[i]
            +
            (
              (R_devoSourceiHigh+R_devoSourceiHighChar)*C_FRACTION_
              -
              (R_devoSourceiHighChar)
            )/M_C
            +
            (R_devoSourceiHigh+R_devoSourceiHighChar)*
            (H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);
          if (alpha > 1.0)
          {
            alpha = 1.0;
          }
          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*Foam::pow((Tp_->ve[i]+1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]+1.0)))/
            (2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]+1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*
            (
              1.715*(Tp_->ve[i]+1.0)*Foam::exp(-9000.0/(Tp_->ve[i]+1.0))*beta_d/
              (1.715*(Tp_->ve[i]+1.0)*Foam::exp(-9000.0/(Tp_->ve[i]+1.0))+beta_d)
            )*surface_area(charFrontRadius_); /* kg/mol*mol/m3*mPar_/s*m2 = kg/s */

          high -= heatRatio(Tp_->ve[i]+1.0)*
            (
              (RCharReac_1S/M_C)*
              (
                (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
               +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
               -LHV[index_C]*M_C
              )
            );

          //LOW~~~~~~~~~~~~~~~~~~~~~~~~~~~
          //alpha = ((R_H2_Su_m[i]+R_devoSourceiLow*Y_H2/Y_G)/M_H2/2.)/(O2_diff_suface);
          R_O2_mole_need =
            R_O2_mole_need_list[i]
            +
            (
              (R_devoSourceiLow+R_devoSourceiLowChar)*C_FRACTION_
              -
              (R_devoSourceiLowChar)
            )/M_C
            +
            (R_devoSourceiLow+R_devoSourceiLowChar)*
            (H_FRACTION_/M_H/4.-O_FRACTION_/M_O/2.);

          alpha = (R_O2_mole_need)/(O2_diff_suface);
          if (alpha > 1.0)
          {
            alpha = 1.0;
          }
          conc_O2_correct = (1.-alpha)*conc_O2; /* molar concentration of O2 in the gas phase */

          D_gAB = 2.03e-10*Foam::pow((Tp_->ve[i]-1.0+TOuterBoundary_)*0.5, 2.0);
          hmia = D_gAB*Foam::pow(0.65, 2.0)/(rb-charFrontRadius_);
          beta_d = hmi*hmia/(hmi+hmia);

          omegaC_S = 2.0*(1.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]-1.0)))
            /(2.0+4.3*Foam::exp(-3390.0/(Tp_->ve[i]-1.0)));

          RCharReac_1S = M_C*omegaC_S*conc_O2_correct*
            (
              1.715*(Tp_->ve[i]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i]-1.0))*beta_d/
             (1.715*(Tp_->ve[i]-1.0)*Foam::exp(-9000.0/(Tp_->ve[i]-1.0))+beta_d)
            )*surface_area(charFrontRadius_);

          low -= heatRatio(Tp_->ve[i]-1.0)*
            (
              (RCharReac_1S/M_C)*
                (
                  (2.0*(omegaC_S-1.0)/omegaC_S)*LHV[index_CO]*M_CO
                 +((2.0-omegaC_S)/omegaC_S)*LHV[index_CO2]*M_CO2
                 -LHV[index_C]*M_C
                )
            );
        }
        else
        {
          //do nothing
        }
      }

      dSdT_star[i-1] = (high-low)/2.0;
    }

    /* --------------------------------------------------------------------- */

    /* Update source term components Su_ and Sp_ */
    /* --------------------------------------------------------------------- */
    /* Picard's method:
      * Su_ = S*-(dS/dT)* * Tp_*
      * Sp_ = (dS/dT)*
      * where*denoted current iteration/time level */
    for (label i = 1; i <= FINE_PARTICLE_MESH_; i++)
    {
      Su_[i-1] = sStar_[i-1]-dSdT_star[i-1]*Tp_->ve[i-1];
      Sp_[i-1] = dSdT_star[i-1];
    }
    /* --------------------------------------------------------------------- */

    /* Update boundary temperature */
    /* --------------------------------------------------------------------- */
    h = k_g/(2.0*rb)*
      (
        2.0
        +
        1.1*Foam::pow((cp_g*mu/k_g), 1.0/3.0)*
        Foam::pow((rho_g*2.0*rb*vel_mag/mu), 0.6)
      );

    a = view_factor*emissivity*Foam::constant::physicoChemical::sigma.value();

    heatHKp = h+
      (
        xi_[layer_1][0]*k_p(layer_1, Tp_->ve[0])
        +
        xi_[layer_2][0]*k_p(layer_2, Tp_->ve[0])
        +
        xi_[layer_3][0]*k_p(layer_3, Tp_->ve[0])
        +
        xi_[layer_4][0]*k_p(layer_4, Tp_->ve[0])
      )/(rb-xc_[0]);

    e =
      -(
        view_factor*emissivity*Foam::constant::physicoChemical::sigma.value()
        *Foam::pow(T_RADIATION_, 4.0)
        +
        h*TInf_
        +
        (
          xi_[layer_1][0]*k_p(layer_1, Tp_->ve[0])
          +
          xi_[layer_2][0]*k_p(layer_2, Tp_->ve[0])
          +
          xi_[layer_3][0]*k_p(layer_3, Tp_->ve[0])
          +
          xi_[layer_4][0]*k_p(layer_4, Tp_->ve[0])
        )*Tp_->ve[0]/(rb-xc_[0])
       );

    q = heatHKp/a;
    delta_0 = 12.0*a*e;
    delta_1 = 27.0*a*Foam::pow(heatHKp, 2.0);

    Q =
      Foam::pow
      (
        0.5*
        (
          delta_1
          +
          Foam::sqrt(Foam::pow(delta_1, 2.0)-4.0*Foam::pow(delta_0, 3.0))
        ),
        (1.0/3.0)
      );

    S = 0.5*Foam::sqrt(1.0/(3.0*a)*(Q+delta_0/Q));
    TOuterBoundary_ = -S+0.5*Foam::sqrt(-4.0*Foam::pow(S, 2.0)+q/S);

    /* --------------------------------------------------------------------- */

    /* Boundary node solution (outer node) */
    /* --------------------------------------------------------------------- */
    //east
    Tf = 0.5*(Tp_->ve[0]+Tp_->ve[1]);
    xif_[layer_1] = 0.5*(xi_[layer_1][0]+xi_[layer_1][1]);
    xif_[layer_2] = 0.5*(xi_[layer_2][0]+xi_[layer_2][1]);
    xif_[layer_3] = 0.5*(xi_[layer_3][0]+xi_[layer_3][1]);
    xif_[layer_4] = 0.5*(xi_[layer_4][0]+xi_[layer_4][1]);
    Af = surface_area(xc_[1]+0.5*(xc_[0]-xc_[1]));

    ae = Af*
      (
        xif_[layer_1]*k_p(layer_1, Tf)
       +xif_[layer_2]*k_p(layer_2, Tf)
       +xif_[layer_3]*k_p(layer_3, Tf)
       +xif_[layer_4]*k_p(layer_4, Tf)
      )/(xc_[0]-xc_[1]);

    ap0 =
      (
        mPar_[layer_1][0]*c_p(layer_1, Tp_->ve[0])
       +mPar_[layer_2][0]*c_p(layer_2, Tp_->ve[0])
       +mPar_[layer_3][0]*c_p(layer_3, Tp_->ve[0])
       +mPar_[layer_4][0]*c_p(layer_4, Tp_->ve[0])
      )/t_step;

    //west, Sp__b = -kA/dx  Su__b = (kA/dx)*T_boundary
    Af = surface_area(rb);

    Sp__b = -Af*
      (
        xi_[layer_1][0]*k_p(layer_1, Tp_->ve[0])
       +xi_[layer_2][0]*k_p(layer_2, Tp_->ve[0])
       +xi_[layer_3][0]*k_p(layer_3, Tp_->ve[0])
       +xi_[layer_4][0]*k_p(layer_4, Tp_->ve[0])
      )/(rb-xc_[0]);

    ap = ap0 /*+ aw*/+ae-Sp__b-Sp_[0]; //no west point
    A->me[0][0] = -ap;
    A->me[0][1] = ae;

    Su__b = Af*
      (
        TOuterBoundary_*
          (
            xi_[layer_1][0]*k_p(layer_1, Tp_->ve[0])
           +xi_[layer_2][0]*k_p(layer_2, Tp_->ve[0])
           +xi_[layer_3][0]*k_p(layer_3, Tp_->ve[0])
           +xi_[layer_4][0]*k_p(layer_4, Tp_->ve[0])
          )/(rb-xc_[0])
      );

    b->ve[0] = -ap0*Tp_->ve[0]-Su__b-Su_[0]; //-ap*Tp0-Su_, Su_ = Af*.... from heat transfer+Su_[0] from reaction
    /* --------------------------------------------------------------------- */
    /* Loop over interior nodes */
    /* --------------------------------------------------------------------- */
    label j = 2;
    for (label i = 2; i <= (FINE_PARTICLE_MESH_-1); i++)
    {
      Tf = 0.5*(Tp_->ve[j-2]+Tp_->ve[j-1]);
      xif_[layer_1] = 0.5*(xi_[layer_1][j-2]+xi_[layer_1][j-1]);
      xif_[layer_2] = 0.5*(xi_[layer_2][j-2]+xi_[layer_2][j-1]);
      xif_[layer_3] = 0.5*(xi_[layer_3][j-2]+xi_[layer_3][j-1]);
      xif_[layer_4] = 0.5*(xi_[layer_4][j-2]+xi_[layer_4][j-1]);

      Af = surface_area(0.5*(xc_[j-2]+xc_[j-1]));

      aw = Af*
        (
          xif_[layer_1]*k_p(layer_1, Tf)
         +xif_[layer_2]*k_p(layer_2, Tf)
         +xif_[layer_3]*k_p(layer_3, Tf)
         +xif_[layer_4]*k_p(layer_4, Tf)
        )/(xc_[j-2]-xc_[j-1]);

      Tf = 0.5*(Tp_->ve[j-1]+Tp_->ve[j]);
      xif_[layer_1] = 0.5*(xi_[layer_1][j-1]+xi_[layer_1][j]);
      xif_[layer_2] = 0.5*(xi_[layer_2][j-1]+xi_[layer_2][j]);
      xif_[layer_3] = 0.5*(xi_[layer_3][j-1]+xi_[layer_3][j]);
      xif_[layer_4] = 0.5*(xi_[layer_4][j-1]+xi_[layer_4][j]);

      Af = surface_area(0.5*(xc_[j-1]+xc_[j]));

      ae = Af*
        (
          xif_[layer_1]*k_p(layer_1, Tf)
         +xif_[layer_2]*k_p(layer_2, Tf)
         +xif_[layer_3]*k_p(layer_3, Tf)
         +xif_[layer_4]*k_p(layer_4, Tf)
        )/(xc_[j-1]-xc_[j]);

      ap0 =
        (
          mPar_[layer_1][j-1]*c_p(layer_1, Tp_->ve[j-1])
         +mPar_[layer_2][j-1]*c_p(layer_2, Tp_->ve[j-1])
         +mPar_[layer_3][j-1]*c_p(layer_3, Tp_->ve[j-1])
         +mPar_[layer_4][j-1]*c_p(layer_4, Tp_->ve[j-1])
        )/t_step;

      ap = ap0+aw+ae-Sp_[j-1];

      A->me[i-1][j-2] = aw;
      A->me[i-1][j-1] = -ap;
      A->me[i-1][j] = ae;
      b->ve[i-1] = -ap0*Tp_->ve[j-1]-Su_[j-1];
      j = j+1;
    }
    /* --------------------------------------------------------------------- */

    /* Boundary node solution (inner node) */
    /* --------------------------------------------------------------------- */
    Tf = 0.5*(Tp_->ve[FINE_PARTICLE_MESH_-2]+Tp_->ve[FINE_PARTICLE_MESH_-1]);
    xif_[layer_1] =
      0.5*(xi_[layer_1][FINE_PARTICLE_MESH_-2]
          +xi_[layer_1][FINE_PARTICLE_MESH_-1]);

    xif_[layer_2] =
      0.5*(xi_[layer_2][FINE_PARTICLE_MESH_-2]
          +xi_[layer_2][FINE_PARTICLE_MESH_-1]);

    xif_[layer_3] =
      0.5*(xi_[layer_3][FINE_PARTICLE_MESH_-2]
          +xi_[layer_3][FINE_PARTICLE_MESH_-1]);

    xif_[layer_4] =
      0.5*(xi_[layer_4][FINE_PARTICLE_MESH_-2]
          +xi_[layer_4][FINE_PARTICLE_MESH_-1]);

    Af = surface_area(0.5*(xc_[FINE_PARTICLE_MESH_-2]+xc_[FINE_PARTICLE_MESH_-1]));
    aw = Af*
      (
        xif_[layer_1]*k_p(layer_1, Tf)
       +xif_[layer_2]*k_p(layer_2, Tf)
       +xif_[layer_3]*k_p(layer_3, Tf)
       +xif_[layer_4]*k_p(layer_4, Tf)
      )/(xc_[FINE_PARTICLE_MESH_-2]-xc_[FINE_PARTICLE_MESH_-1]);

    ap0 =
      (
        mPar_[layer_1][FINE_PARTICLE_MESH_-1]
        *c_p(layer_1, Tp_->ve[FINE_PARTICLE_MESH_-1])
       +mPar_[layer_2][FINE_PARTICLE_MESH_-1]
        *c_p(layer_2, Tp_->ve[FINE_PARTICLE_MESH_-1])
       +mPar_[layer_3][FINE_PARTICLE_MESH_-1]
        *c_p(layer_3, Tp_->ve[FINE_PARTICLE_MESH_-1])
       +mPar_[layer_4][FINE_PARTICLE_MESH_-1]
        *c_p(layer_4, Tp_->ve[FINE_PARTICLE_MESH_-1])
      )/t_step;

    ap = ap0+aw /*+ ae*/-Sp_[FINE_PARTICLE_MESH_-1];
    A->me[FINE_PARTICLE_MESH_-1][FINE_PARTICLE_MESH_-2] = aw;
    A->me[FINE_PARTICLE_MESH_-1][FINE_PARTICLE_MESH_-1] = -ap;
    b->ve[FINE_PARTICLE_MESH_-1] =
      -ap0*Tp_->ve[FINE_PARTICLE_MESH_-1]-Su_[FINE_PARTICLE_MESH_-1];
    /* --------------------------------------------------------------------- */

    for (j = 1; j <= FINE_PARTICLE_MESH_; j++)
    {
      if (j == 1)
      {
        ThermalTimeScale[j-1] = Foam::pow(rb-xc_[j-1], 2.0)/
          (
            (
              xi_[layer_1][j-1]*k_p(layer_1, Tp_->ve[j-1])
             +xi_[layer_2][j-1]*k_p(layer_2, Tp_->ve[j-1])
             +xi_[layer_3][j-1]*k_p(layer_3, Tp_->ve[j-1])
             +xi_[layer_4][j-1]*k_p(layer_4, Tp_->ve[j-1])
            )/
            (
              (
                mPar_[layer_1][j-1]*c_p(layer_1, Tp_->ve[j-1])
               +mPar_[layer_2][j-1]*c_p(layer_2, Tp_->ve[j-1])
               +mPar_[layer_3][j-1]*c_p(layer_3, Tp_->ve[j-1])
               +mPar_[layer_4][j-1]*c_p(layer_4, Tp_->ve[j-1])
              )/VCell_[j-1]
            )
          );
      }
      else
      {
        ThermalTimeScale[j-1] = Foam::pow(xc_[j-2]-xc_[j-1], 2.0)/
          (
            (
              xi_[layer_1][j-1]*k_p(layer_1, Tp_->ve[j-1])
             +xi_[layer_2][j-1]*k_p(layer_2, Tp_->ve[j-1])
             +xi_[layer_3][j-1]*k_p(layer_3, Tp_->ve[j-1])
             +xi_[layer_4][j-1]*k_p(layer_4, Tp_->ve[j-1])
            )/
            (
              (
                mPar_[layer_1][j-1]*c_p(layer_1, Tp_->ve[j-1])
               +mPar_[layer_2][j-1]*c_p(layer_2, Tp_->ve[j-1])
               +mPar_[layer_3][j-1]*c_p(layer_3, Tp_->ve[j-1])
               +mPar_[layer_4][j-1]*c_p(layer_4, Tp_->ve[j-1])
              )/VCell_[j-1]
            )
          );
      }
    }

    /* Solve matrix system */
    /* --------------------------------------------------------------------- */
    // LU = m_get(A->m, A->n);
    // LU = m_copy(A, LU);
    // pivot = px_get(A->m);
    // LUfactor(LU, pivot);
    // Tp_ = LUsolve(LU, pivot, b, Tp_);
    // PX_FREE(pivot);
    // M_FREE(LU);
    EigenCGFunc(A, Tp_, b);

        /* Output to file (and possibly to screen) */
    /* --------------------------------------------------------------------- */
    if (iter_count % output_every == 0)
    {
      if (VERBOSITY_)
      {
        printf("Outer boundary temperature: %f\n", TOuterBoundary_);
        // Info << "Outer boundary temperature:" << nl;
        // Info << "Outer boundary temperature: " << nl
        //      << TOuterBoundary_ << nl;
      }
    }

    if ((iter_count % output_every) == 0)
    {
      if (VERBOSITY_)
      {
        printf("Tp_ at time = %e s:\n", tPar);
        for (label i = 1; i <= FINE_PARTICLE_MESH_; i++)
        {
          printf("%i:%f, ", i, Tp_->ve[i-1]);
        }
        printf("\n");
        printf("mPar_ at time = %e s:\n", tPar);
        for (j = 1; j <= FINE_PARTICLE_MESH_; j++)
        {
          printf(" | ");
          for (label i = 1; i <= 4; i++)
          {
            printf("%e ", mPar_[i-1][j-1]);
          }
        }
        printf("\n");
      }
    }

    // tPar = tPar+dt;
    iter_count = iter_count+1;
    /* --------------------------------------------------------------------- */

    //- Update time
    tPar = tPar + t_step;
    // Info << "tPar+t_step " << tPar << nl;
    // cumTime_ = cumTime_ + t_step;

    /* Update the masses */
    // Info << "before update mass loop" << nl;
    /* --------------------------------------------------------------------- */
    for (label i = 1; i <= FINE_PARTICLE_MESH_; i++)
    {
      // Info << "update mass loop" << nl;
      if (tPar == t_step)
      {
        // Info << "before 0 " << nl;
        charYield_[i-1] = charFromDevol_[i-1]/(totalDevol_[i-1]+VSMALL);
        // Info << "after 0 " << nl;
      }
      else
      {
        // Info << "before 0 1" << nl;
        // Info << "totalDevol_[i-1] " << totalDevol_[i-1] << nl;
        charYield_[i-1] = charFromDevol_[i-1]/(totalDevol_[i-1]+VSMALL);
        dmdt[layer_1][i-1] = -Rb_[drying_front][i-1]/MOISTURE_FRACTION_WB;
        // Info << "before 0 1" << nl;
      }

      // Info << "MOISTURE_FRACTION_WB is " << MOISTURE_FRACTION_WB << nl;

      // Info << "before 1 " << nl;

      if (mPar_[layer_2][i-1]-Rb_[devol_front][i-1]*t_step <= 0.)
      {
        // Info << "t_step is " << t_step << nl;
        Rb_[devol_front][i-1] = mPar_[layer_2][i-1]/t_step;
      }

      // Info << "after 1 " << nl;

      dmdt[layer_2][i-1] =
        Rb_[drying_front][i-1]/MOISTURE_FRACTION_DB_-Rb_[devol_front][i-1];

      // Info << "before 2 " << nl;

      if (combustion_flag == 1)
      {
        if
        (
          mPar_[layer_3][combustion_layer]-Rb_[char_front][combustion_layer]*t_step
          <=
          0.
        )
        {
          Rb_[char_front][combustion_layer] = mPar_[layer_3][combustion_layer]/t_step;
        }
      }

      // Info << "after 2 " << nl;

      dmdt[layer_3][i-1] =
        Rb_[char_produced][i-1]-Rb_[char_front][i-1];

      if (charYield_[i-1] > 0.0)
      {
        dmdt[layer_4][i-1] =
            Rb_[char_front][i-1]*ASH_MASS_FRACTION_/charYield_[i-1];
      }
      else
      {
        dmdt[layer_4][i-1] = 0.0;
      }

      mPar_[layer_1][i-1] += dmdt[layer_1][i-1]*t_step;
      mPar_[layer_2][i-1] += dmdt[layer_2][i-1]*t_step;
      mPar_[layer_3][i-1] += dmdt[layer_3][i-1]*t_step;
      mPar_[layer_4][i-1] += dmdt[layer_4][i-1]*t_step;

      for (j = 1; j <= 4; j++)
      {
        if (mPar_[j-1][i-1] < 0.0)
        {
          mPar_[j-1][i-1] = 0.0;
        }
      }
    }

    if (tPar + t_step >= dt)
    {
        t_step = dt - tPar;
    }

    if (t_step == 0) {
      break;
    }

    /* --------------------------------------------------------------------- */
  }
  // To sum up all the mPars
  forAll(mPar_,i)
  {
    mParSum_ += sum(mPar_[i]);
  }
}

// ************************************************************************* //
