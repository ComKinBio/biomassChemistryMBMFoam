template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReadLHVM(scalar YC)
{
    /* Lower heating values (LHVs) */
    LHV[index_CO2] = 0.0;                                         /* per definition */
    LHV[index_H2] = 120.0e6;                                      /* J/kg D&D p. 26 */
    LHV[index_CO] = 10.25e6;                                      /* J/kg D&D p. 26 */
    LHV[index_H2O] = 0.0;                                         /* per definition */
    LHV[index_CH4] = 50.0e6;                                      /* J/kg D&D p. 26 */
    LHV[index_C2H4] = 47.7e6;                                     /* J/kg D&D p. 26 */
    LHV[index_C10H8] = 39.48e6;                                   /* J/kg Alberto's code */
    LHV[index_C] = 32.79e6;                                       /* J/kg Alberto's code */
    LHV[index_C2H6] = 47.6e6;                                     /* J/kg D&D p. 26 */
    LHV[index_C3H8] = 46.3e6;                                     /* J/kg D&D p. 26 */
    LHV[index_C6H6O] = 32.45e6;                                   /* J/kg http://www.nist.gov/data/PDFfiles/jpcrd6.pdf OR http://pac.iupac.org/publications/pac/pdf/1961/pdf/0201x0125.pdf */
    LHV[index_volatiles] = (18.93e6 - 32.57e6 * YC) / (1.0 - YC); /* J/kg Alberto's code */
    /* Molar masses */
    M[index_CO2] = 44.01e-3;      /* kg/mole */
    M[index_H2] = 2.0 * 1.008e-3; /* kg/mole */
    M[index_CO] = 28.010e-3;      /* kg/mole */
    M[index_H2O] = 18.01528e-3;   /* kg/mole */
    M[index_CH4] = 16.04e-3;      /* kg/mole */
    M[index_C2H4] = 28.05e-3;     /* kg/mole */
    M[index_C10H8] = 128.17e-3;   /* kg/mole */
    M[index_C] = 12.011e-3;       /* kg/mole */
    M[index_C2H6] = 30.07e-3;     /* kg/mole */
    M[index_C3H8] = 44.1e-3;      /* kg/mole */
    M[index_C6H6O] = 94.11e-3;    /* kg/mole */
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::readMeshInput
(
   const dictionary &particleInfoDict_
)
{

  #define MESHINPUT particleInfoDict_.subDict("meshInput").lookup

  MESHINPUT("VERBOSITY") >> VERBOSITY_;
  MESHINPUT("SAVE_EVERY") >> SAVE_EVERY_;
  MESHINPUT("TIME_DT") >> TIME_DT_;
  MESHINPUT("MOISTURE_FRACTION_DB") >> MOISTURE_FRACTION_DB_;
  MESHINPUT("PARTICLE_RADIUS") >> PARTICLE_RADIUS_;
  MESHINPUT("DRYING_SHRINKAGE_VOLUME") >> DRYING_SHRINKAGE_VOLUME_;
  MESHINPUT("T_RADIATION") >> T_RADIATION_;
  MESHINPUT("T_REACTOR") >> T_REACTOR_;
  MESHINPUT("T_P_INIT") >> T_P_INIT_;
  MESHINPUT("MASS_FRAC_O2") >> MASS_FRAC_O2_;
  MESHINPUT("DENSITY_INITIAL") >> DENSITY_INITIAL_;
  MESHINPUT("FINE_PARTICLE_MESH") >> FINE_PARTICLE_MESH_;
  MESHINPUT("FINE_PARTICLE_MESH_DOUBLE") >> FINE_PARTICLE_MESH_DOUBLE_;
  MESHINPUT("K_P_WET") >> K_P_WET_;
  MESHINPUT("K_P_DRY") >> K_P_DRY_;
  MESHINPUT("K_P_CHAR") >> K_P_CHAR_;
  MESHINPUT("ASH_MASS_FRACTION") >> ASH_MASS_FRACTION_;
  MESHINPUT("C_FRACTION") >> C_FRACTION_;
  MESHINPUT("H_FRACTION") >> H_FRACTION_;
  MESHINPUT("O_FRACTION") >> O_FRACTION_;
  MESHINPUT("CHARCONVERSION_SHRINKAGE_VOLUME") >> CHARCONVERSION_SHRINKAGE_VOLUME_;
  MESHINPUT("DEVOL_SHRINKAGE_VOLUME") >> DEVOL_SHRINKAGE_VOLUME_;
  MESHINPUT("P_ATM") >> P_ATM_;
  MESHINPUT("MASS_FRAC_H2O") >> MASS_FRAC_H2O_;
  MESHINPUT("MASS_FRAC_H2") >> MASS_FRAC_H2_;
  MESHINPUT("MASS_FRAC_CO2") >> MASS_FRAC_CO2_;
  MESHINPUT("MY_MAXIMUM_NEWTON_ITER") >> MY_MAXIMUM_NEWTON_ITER_;
  MESHINPUT("N_SPECIES") >> N_SPECIES_;
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::readDevKin
(
   const dictionary &particleInfoDict_
)
{
  #define SIGKIN particleInfoDict_.subDict("sigKin").lookup

  /* pyrolysis kinetics parameters */

  SIGKIN("DevKinA1") >> DevKinA1_;
  SIGKIN("DevKinE1") >> DevKinE1_;
  SIGKIN("DevKinA2") >> DevKinA2_;
  SIGKIN("DevKinE2") >> DevKinE2_;
  SIGKIN("DevKinA3") >> DevKinA3_;
  SIGKIN("DevKinE3") >> DevKinE3_;

  /* drying kinetics parameters */
  SIGKIN("DryKinA") >> DryKinA_;
  SIGKIN("DryKinE") >> DryKinE_;
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::readCoeff(const Time &time)
{
    // dictionary particleInfoDict_("constant/particleInfoDict");
    // Info << "particleInfoDict_ is " << particleInfoDict_ << nl;
    const dictionary& particleInfoDict_
    (
        IFstream("constant/particleInfoDict")()
    );

    Info << particleInfoDict_ << nl;
    readMeshInput(particleInfoDict_);
    readDevKin(particleInfoDict_);

    MOISTURE_FRACTION_WB = MOISTURE_FRACTION_DB_/(1.0+MOISTURE_FRACTION_DB_);
    ASH_MASS_FRACTION_DB = ASH_MASS_FRACTION_/(1.0-MOISTURE_FRACTION_WB);

}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::initFields()
{
    // initialize some properties of the particle
    material = 1;
    shape = sphere;
    r_init = PARTICLE_RADIUS_;
    Info << "PARTICLE_RADIUS_ " << PARTICLE_RADIUS_ << nl;
    L1 = 2.0*r_init;
    L2 = 2.0*r_init;
    L3 = 2.0*r_init;

    // initial mass of the biomass particle
    m_init = rhoPar(layer_1, T_P_INIT_)*volume_of_element(0.0, r_init);

    // initial volume of the biomass particle
    particle_volume0 = volume_of_element(0.0, r_init);

    LHV.resize(N_SPECIES_, 0.0);
    M.resize(N_SPECIES_-1, 0.0);
    ReadLHVM(0.18); /* assumes a char fraction of 0.18 (only relevant for volatile LHV) */

    #define LISTINIT2D(NAME,SIZE1,SIZE2) NAME.resize(SIZE1); \
    forAll(NAME,i)                                         \
    {                                                      \
        NAME[i].resize(SIZE2, 0.0);                             \
    }

    LISTINIT2D(xi_,4,FINE_PARTICLE_MESH_);
    LISTINIT2D(QReaction_,3,FINE_PARTICLE_MESH_);
    LISTINIT2D(Rb_,4,FINE_PARTICLE_MESH_);
    LISTINIT2D(RCharReac_,4,FINE_PARTICLE_MESH_);
    LISTINIT2D(mPar_,4,FINE_PARTICLE_MESH_);
    LISTINIT2D(mLast_,4,FINE_PARTICLE_MESH_);
    LISTINIT2D(QEquilibrate_,3,FINE_PARTICLE_MESH_);
    LISTINIT2D(RDevolReac_,3,FINE_PARTICLE_MESH_);
    LISTINIT2D(dmdt,4,FINE_PARTICLE_MESH_); //pass

    VCell_.resize(FINE_PARTICLE_MESH_, 0.0);
    VCharLast_.resize(FINE_PARTICLE_MESH_, 0.0);
    VChar_.resize(FINE_PARTICLE_MESH_, 0.0);
    xc_.resize(FINE_PARTICLE_MESH_, 0.0);
    Sp_.resize(FINE_PARTICLE_MESH_, 0.0);
    Su_.resize(FINE_PARTICLE_MESH_);
    xif_.resize(FINE_PARTICLE_MESH_); //
    omegaC_.resize(FINE_PARTICLE_MESH_, 0.0);
    sStar_.resize(FINE_PARTICLE_MESH_, 0.0); //
    dSdT_star.resize(FINE_PARTICLE_MESH_, 0.0); //
    charFromDevol_.resize(FINE_PARTICLE_MESH_, 0.0); //
    totalDevol_.resize(FINE_PARTICLE_MESH_, 0.0); //
    RProduct_.resize(7, 0.0); // TO DO as different model
    charYield_.resize(FINE_PARTICLE_MESH_, 0.0); //
    betaR_.resize(4, 0.0);
    R_O2_mole_need_list.resize(FINE_PARTICLE_MESH_, 0.0); //
    ThermalTimeScale.resize(FINE_PARTICLE_MESH_, 0.0); //

    /* Assemble matrix system */
    A = m_get(FINE_PARTICLE_MESH_, FINE_PARTICLE_MESH_);
    Tp_ = v_get(FINE_PARTICLE_MESH_);
    b = v_get(FINE_PARTICLE_MESH_);

    /* 	ap*Tp_ = aw*Tw+ae*Te+ap0*Tp0+Su_
    *       ap = ap0+aw+aE-Sp_
    *       ap0 = rho*cp*deltax/deltat
    *       aw = kw/deltax_wp
    *       ae = ke/deltax_pe
    *       b = S*deltax = Su_+Sp_*Tp_		*/
    /* --- formulated as a matrix problem:
    *       aw*Tw-ap*Tp_+ae*Te =-ap0*Tp0-Su_	*/

    /* Set some parameters */
    output_every = SAVE_EVERY_; //((label) (1.0/dt));
    // t_stop = TIME_TOTAL_;

    /* Initialize variables */
    // particle time
    // tPar = 0.0;
    // set the particle temperature to be zero
    v_zero(Tp_); // take old values or not? No!!

    // initialize the iteration number to be 0
    // this was used for print out result before
    iter_count = 0;

    // read the particle radius
    // this value needs to be updated!!
    rb = PARTICLE_RADIUS_; //

    /*Print out KP*/
    // Info << "K_P_WET_: " << K_P_WET_ << nl;
    // Info << "K_P_DRY_: " << K_P_DRY_ << nl;
    // Info << "K_P_CHAR_: " << K_P_CHAR_ << nl;

    /* Initialize fields */
    if (VERBOSITY_)
    {
    printf("Init Tp_:\n");
    }

    for (unsigned int i = 0; i < Tp_->max_dim; i++)
    {
        // set the initial temperature to the Tp_ vector
        // in IO.C file, this value should not be set again, but read
        Tp_->ve[i] = T_P_INIT_;
        if (VERBOSITY_)
        {
            // printf("%i: %f ", i, Tp_->ve[i]);
            Info << i << ": " << Tp_->ve[i] << " ";
        }
    }

    // Set the initial boundary temperature
    TOuterBoundary_ = T_P_INIT_;
    if (VERBOSITY_)
    {
        printf("\n");
    }

    /* Particle is meshed from the outside in */
    if (VERBOSITY_)
    {
        printf("Particle mesh:\n");
    }

    // calculate the cell Volume
    scalar cellVol = volume_of_element(0.0, rb)/FINE_PARTICLE_MESH_DOUBLE_;

    // volume_sum at the beginning is equal to one cell volume
    volume_sum = 0.5*cellVol;

    // Set the cell center location from the outmost cell using volume
    // cell number 1 is the one close to the particle surface
    for (label j = FINE_PARTICLE_MESH_; j >= 1; j--)
    {
        // xc_ means the cell center location
        Info << "volume_sum " << volume_sum << nl;
        xc_[j-1] = radius_from_volume(volume_sum);
        volume_sum = volume_sum+cellVol;
        if (VERBOSITY_)
        {
            printf("%i: %e-", j, xc_[j-1]);
        }
        charFromDevol_[j-1] = 0.0;
        totalDevol_[j-1] = 0.0;
    }

    if (VERBOSITY_)
    {
        printf("\n");
    }

    // initialize the total mass to be zero at the beginning
    // in IO.C file, this value should not be set again, but read
    total_mass = 0.0;
    // four layers
    for (label i = 0; i <= 3; i++)
    {
        // each cell in layer i
        for (label j = 0; j < FINE_PARTICLE_MESH_; j++)
        {
            // if it is layer 1
            if (i == layer_1)
            {
                // initial composition of wet wood is set as 1 for all the cells
                xi_[i][j] = 1.0;
            }
            else
            {
                // other layer (dry wood, char, and ash) does not exit yet
                xi_[i][j] = 0.0;
            }

            // mass of each layer is initialize
            // mass = xi_ (fraction) * rho * volume
            mPar_[i][j] =
                xi_[i][j]*rhoPar(i, Tp_->ve[i])*cellVol;

            printf("x %e ", xi_[i][j]);
            printf("mPar_ %e ", mPar_[i][j]);
            // initialize the total mass before any reactions
            total_mass = total_mass+mPar_[i][j];
        }
    }
    printf("\n");

    /* Set gas phase data at the reactor temperature */
    TInf_ = T_REACTOR_;

    // gas density
    rho_g = (3.550951e2*Foam::pow(TInf_, -1.000866));

    // viscosity?
    mu =
        1.716e-5*Foam::pow((TInf_/273.11), 3./2.)
        *(273.11+110.56)/(TInf_+110.56);

    // gas heat capacity
    cp_g = 1.0e3*
    (
        (-0.000000000153767)*Foam::pow(TInf_, 3.0)
        +0.000000363209494*Foam::pow(TInf_, 2.0)
        +(-0.000051486068111)*TInf_+1.009174922600606
    ); /* polyfit from table data for N2 on the interval 500-1300 K */

    k_g = (-0.000000000003493)*Foam::pow(TInf_, 3.0)
            +0.000000003319264*Foam::pow(TInf_, 2.0)
            +0.000060059499759*TInf_+0.008533051948052;
    /* polyfit from table data for N2 on the interval 500-1300 K */

    if (VERBOSITY_)
    {
        printf("Gas phase data to use at Tp_ = %f\n", TInf_);
        printf("  rho  = %e\n", rho_g);
        printf("  mu   = %e\n", mu);
        printf("  cp_g = %e\n", cp_g);
        printf("  k_g  = %e\n", k_g);
    }

    // velocity?
    vel_mag = 0.5; /* (mPar_/s) */

    /* molar concentration of O2 in the gas phase */
    conc_O2 = MASS_FRAC_O2_*(M_N2/M_O2)*P_ATM_/(8.3145*TInf_);


    //  initialize volume of each components
    for (label j = FINE_PARTICLE_MESH_; j >= 1; j--)
    {
        VCell_[j-1] = cellVol;
        VChar_[j-1] = 0.0;
        VCharLast_[j-1] = 0.0;

        for (label i = 1; i <= 4; i++)
        {
            // mLast_ is the mass from the last time step
            mLast_[i-1][j-1] = mPar_[i-1][j-1];
        }
        // The thermal time scale should not be initialized again in IO
        if (j == 1)
        {
            ThermalTimeScale[j-1] =
            Foam::pow(rb-xc_[0], 2.0)
            /(k_p(layer_1, Tp_->ve[j-1])
            /(c_p(layer_1, Tp_->ve[j-1])
            *rhoPar(layer_1, Tp_->ve[j-1])));
        }
        else
        {
            ThermalTimeScale[j-1] =
            Foam::pow(xc_[j-2]-xc_[j-1], 2.0)
            /(k_p(layer_1, Tp_->ve[j-1])
            /(c_p(layer_1, Tp_->ve[j-1])
            *rhoPar(layer_1, Tp_->ve[j-1])));
        }
    }

    density_drywood =
        DENSITY_INITIAL_
        *(1.0-MOISTURE_FRACTION_WB)
        /(1.0-DRYING_SHRINKAGE_VOLUME_);

    combustion_layer = 99999999;

    Info << "A->max_m " << A->max_m << nl;
    Info << "A->max_n " << A->max_n << nl;
    Info << "A->max_size " << A->max_size << nl;

    Info << "xi_" << nl << xi_ << nl;

};


template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::linInterp
(
    scalar X, label n_elements_in_list
)
{
    label i = 0;
    scalar my_value;

    if (X <= X_list[0])
    {
        my_value = Y_list[0];
    }
    else if (X >= X_list[n_elements_in_list - 1])
    {
        my_value = Y_list[n_elements_in_list - 1];
    }
    else
    {
        while (!((X >= X_list[i - 1]) && (X < X_list[i])))
        {
            i = i + 1;
        }
        my_value = Y_list[i - 1]
                    +
                    (Y_list[i] - Y_list[i - 1]) * (X - X_list[i - 1])
                    / (X_list[i] - X_list[i - 1]);
    }

    return my_value;
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::EigenCGFunc
(
    MAT* A, VEC* T, VEC* b
)
{
  // Eigen solver
  std::vector<Eigen::Triplet<double>> tripletList;
  int estimation_of_entries = 4+3*(FINE_PARTICLE_MESH_-2);
  tripletList.reserve(estimation_of_entries);

  Eigen::SparseMatrix<double> S(FINE_PARTICLE_MESH_,FINE_PARTICLE_MESH_);
  for (size_t i = 0; i < FINE_PARTICLE_MESH_; i++) 
  {
    for (size_t j = 0; j < FINE_PARTICLE_MESH_; j++) 
    {
      if (A->me[i][j]!=0) 
      {
        tripletList.push_back(Eigen::Triplet<double>(i,j,A->me[i][j]));
      }
    }
  }

  S.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::VectorXd bEigen(FINE_PARTICLE_MESH_), xEigen(FINE_PARTICLE_MESH_);
  for (size_t j = 0; j < FINE_PARTICLE_MESH_; j++) {
    bEigen[j] = b->ve[j];
    xEigen[j] = T->ve[j];
  }

  

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> Eigsolver;
  // Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> Eigsolver;
  Eigsolver.compute(S);
  xEigen = Eigsolver.solve(bEigen);

//   VEC *TReturn = v_get(FINE_PARTICLE_MESH_);
  for (size_t i = 0; i < FINE_PARTICLE_MESH_; i++) {
    T->ve[i] = xEigen[i];
  }

//   return T;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::c_p
(
  label layer, scalar T
)
{
    /* returns c_p in J/kg,K for the particle layer */
    /* N.B.: currently, all materials use the same expressions! */
    scalar heat_capacity = 0.0;
    scalar c_p_wood;
    scalar c_p_water = 4.309e3; /* c_p for water in J/kg,K */
    scalar A;                   /* correction factor in J/kg,K */
    scalar Y_mois = 1.0-1.0/(1.0+MOISTURE_FRACTION_DB_);
    //     scalar hcgas;

    if (layer == layer_1)
    {
        /* moist wood */
        c_p_wood = 4.206*T-37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        A = 1.0e3*
            ((0.02355*T-1.320*Y_mois/(1.0-Y_mois)-6.191)*Y_mois/(1.0-Y_mois));
        heat_capacity = c_p_wood*(1.0-Y_mois)+c_p_water*Y_mois+A;
    }
    else if (layer == layer_2)
    {
        /* dry wood */
        heat_capacity = 4.206*T-37.7; /* eq. (24.i) Thunman et al., Energy & Fuels 2001 */
        // heat_capacity = 1000.*(0.1031+0.003867*T); /* Properties of Wood for Combustion Analysis */
    }
    else if (layer == layer_3)
    {
        /* char */
        heat_capacity = -334.0+4410.0e-3*T-3160.0e-6*Foam::pow(T, 2.0)
            +1010.0e-9*Foam::pow(T, 3.0)-119.0e-12*Foam::pow(T, 4.0); /* Table 5 Thunman et al., Energy & Fuels 2001 */
        //heat_capacity = 1000.*(1.39+0.00036*T); /* Properties of Wood for Combustion Analysis */
    }
    else if (layer == layer_4)
    {
        heat_capacity = 754.0+0.586*(T-273.15); /* Sven */
    }
    else
    {
        printf("ERROR in c_p function!!!\n");
    }

    return heat_capacity;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::k_p
(
  label layer, scalar T
)
{
    /* returns heat conductivity (W/m,K) for each layer as a function of temperature */
    scalar heat_conductivity = 0.0, kgas = 0.0;

    if (layer == layer_1)
    {
        heat_conductivity = K_P_WET_; //
    }
    else if (layer == layer_2)
    {
        heat_conductivity = K_P_DRY_; //
    }
    else if (layer == layer_3)
    {
        //         heat_conductivity = 0.052; //K. W. Ragland, D. J. Aerts, Properties of Wood for Combustion Analysis Wheast cites "R. C. (1976). Handbook of Chemistry and Physics, CRC Press, Cleveland, Ohio, p. E-5."
        heat_conductivity = K_P_CHAR_;
    }
    else if (layer == layer_4)
    {
        kgas = (-0.000000000003493)*Foam::pow(T, 3.0)
            +0.000000003319264*Foam::pow(T, 2.0)
            +0.000060059499759*T+0.008533051948052;

        heat_conductivity = 1.03*(1.0-0.65)+0.65*kgas; /* Sven */
    }
    else
    {
        printf("ERROR in heat conductivity function!!\n");
    }

    return heat_conductivity;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::rhoPar
(
  label layer, scalar T
)
{
    /* returns density (kg/m3) for each layer as a function of temperature */
    scalar density = 0.0;
    scalar dgas;

    if (material == 1) /* in this version of the code there is only one material implemented in this function */
    {
        if (layer == layer_1)
        {
            density = DENSITY_INITIAL_;
        }
        else if (layer == layer_2)
        {
            density = (1.0-MOISTURE_FRACTION_WB)*DENSITY_INITIAL_;
        }
        else if (layer == layer_3)
        {
            density = 350.0;
        }
        else if (layer == layer_4)
        {
            dgas = 101325.0*(2.0*14.00672*1.0e-3)/(8.3145*T);
            density = 2000.0*(1.0-0.65)+0.65*dgas; /* Sven */
        }
        else
        {
            printf("ERROR in density function!!\n");
        }
    }
    else
    {
        printf("ERROR in density function!!\n");
    }

    return density;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::gas_cp
(
  label index, scalar T
)
{
    /* returns heat capacity of the gas in J/mol,K */
    label i;
    scalar A, B, C, D, E;
    scalar temp = 0.0;
    scalar cp = 0.0;

    scalar T_list_C10H8[19] =
    {50.0, 100.0, 150.0, 200.0, 273.15, 298.15, 300.0, 400.0, 500.0, 600.0,
    700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0};

    scalar Cp_list_C10H8[19] =
    {36.18, 47.50, 63.89, 84.99, 120.52, 133.02, 133.94, 181.16, 220.70,
    252.37, 277.77, 298.43, 315.50, 329.77, 341.8, 352.0, 360.8, 368.2, 374.7};

    scalar T_list_C2H6[31] =
    {100.0, 200.0, 298.15, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0,
    1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600.0, 1700.0, 1800.0,
    1900.0, 2000.0, 2100.0, 2200.0, 2300.0, 2400.0, 2500.0, 2600.0, 2700.0,
    2800.0, 2900.0, 3000.0};

    scalar Cp_list_C2H6[31] =
    {35.70, 42.30, 52.49, 52.71, 65.46, 77.94, 89.19, 99.14, 107.94, 115.71,
    122.55, 128.55, 133.80, 138.39, 142.40, 145.90, 148.98, 151.67, 154.04,
    156.14, 158.00, 159.65, 161.12, 162.43, 163.61, 164.67, 165.63, 166.49,
    167.28, 168.00, 168.65};

    scalar T_list_C3H8[19] =
    {50.0, 100.0, 150.0, 200.0, 273.15, 298.15, 300.0, 400.0, 500.0, 600.0,
    700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0};

    scalar Cp_list_C3H8[19] =
    {34.06, 41.30, 48.79, 56.07, 68.74, 73.60, 73.93, 94.01, 112.59, 128.70,
    142.67, 154.77, 165.35, 174.60, 182.67, 189.74, 195.85, 201.21, 205.89};

    scalar T_list_C6H6O[19] =
    {50.0, 100.0, 150.0, 200.0, 273.15, 298.15, 300.0, 400.0, 500.0, 600.0,
    700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0};

    scalar Cp_list_C6H6O[19] =
    {33.91, 41.38, 54.19, 69.65, 94.61, 103.22, 103.86, 135.79, 161.91, 182.48,
    198.84, 212.14, 223.19, 232.49, 240.41, 247.20, 253.06, 258.12, 262.52};
    //     scalar BensonBuss_H = 0.85;
    //     scalar BensonBuss_C = 3.75;
    //     scalar BensonBuss_O = 3.40;

    if (index == index_CO2)
    {
        if (T < 1200.0)
        {
            /* valid to 1200K */
            A = 24.99735;
            B = 55.18696;
            C = -33.69137;
            D = 7.948387;
            E = -0.136638;
        }
        else
        {
            /* valid from 1200 to 6000 K */
            A = 58.16639;
            B = 2.720074;
            C = -0.492289;
            D = 0.038844;
            E = -6.447293;
        }
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_H2)
    {
        if (T < 1000.0)
        {
            /* valid from 298 to 1000K */
            A = 33.066178;
            B = -11.363417;
            C = 11.432816;
            D = -2.772874;
            E = -0.158558;
        }
        else
        {
            /* valid from 1000 to 2500K */
            A = 18.563083;
            B = 12.257357;
            C = -2.859786;
            D = 0.268238;
            E = 1.977990;
        }
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_CO)
    {
        /* valid to 1300K */
        A = 25.56759;
        B = 6.096130;
        C = 4.054656;
        D = -2.671301;
        E = 0.131021;
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_H2O)
    {
        /* valid from 500 to 1700 K */
        A = 30.09200;
        B = 6.832514;
        C = 6.793435;
        D = -2.534480;
        E = 0.082139;
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_CH4)
    {
        /* valid to 1300K */
        A = -0.703029;
        B = 108.4773;
        C = -42.52157;
        D = 5.862788;
        E = 0.678565;
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_C2H4)
    {
        if (T < 1200.0)
        {
            /* valid from 298 to 1200K */
            A = -6.387880;
            B = 184.4019;
            C = -112.9718;
            D = 28.49593;
            E = 0.315540;
        }
        else
        {
            /* valid from 1200 K to 6000 K */
            A = 106.5104;
            B = 13.73260;
            C = -2.628481;
            D = 0.174595;
            E = -26.14469;
        }
        /* Shomate Equation from NIST Webbook */
        temp = T/1000.0;
        cp = A+B*temp+C*Foam::pow(temp, 2.0)+D*Foam::pow(temp, 3.0)+
            E/Foam::pow(temp, 2.0);
    }
    else if (index == index_C10H8)
    {
        for (i = 1; i <= 19; i++)
        {
            X_list[i-1] = T_list_C10H8[i-1];
            Y_list[i-1] = Cp_list_C10H8[i-1];
        }
        cp = linInterp(T, 19);
    }
    else if (index == index_C)
    {
        /* specific heat of graphite (solid)-valid from 250 to 3000 K (taken from NIST Webbook)
        *       %     cp = 0.538657+9.11129e-6*T-90.2725*T^(-1)-43449.3*T^(-2)+1.59309e7*T^(-3)-1.43688e9*T^(-4); % cal/g*K
        *       %     cp = cp*12.0107;    % (*g/mol) --> cal/mol,K
        *       %     cp = cp*4.184;      % (*J/cal) --> J/mol,K */
        cp =
            (
                0.538657+9.11129e-6*T-90.2725*Foam::pow(T, -1.0)
                -43449.3*Foam::pow(T, -2.0)+1.59309e7*Foam::pow(T, -3.0)
                -1.43688e9*Foam::pow(T, -4.0)
            )*12.0107*4.184; /* J/mol*K */
    }
    else if (index == index_C2H6)
    {
        for (i = 1; i <= 31; i++)
        {
            X_list[i-1] = T_list_C2H6[i-1];
            Y_list[i-1] = Cp_list_C2H6[i-1];
        }
        cp = linInterp(T, 31);
    }
    else if (index == index_C3H8)
    {
        for (i = 1; i <= 19; i++)
        {
            X_list[i-1] = T_list_C3H8[i-1];
            Y_list[i-1] = Cp_list_C3H8[i-1];
        }
        cp = linInterp(T, 19);
    }
    else if (index == index_C6H6O)
    {
        for (i = 1; i <= 19; i++)
        {
            X_list[i-1] = T_list_C6H6O[i-1];
            Y_list[i-1] = Cp_list_C6H6O[i-1];
        }
        cp = linInterp(T, 19);
    }
    else
    {
        printf("ERROR in gas_cp function!!! (1)\n");
    }

    return cp;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::deltaHvap(scalar T)
{
    /* returns heat of vaporization of water in J/kg for T in Kelvin */
    /* on the interval 273-433K */
    label i;
    scalar T_list[9] =
        {273.15, 298.15, 313.15, 333.15, 353.15,
        373.15, 393.15, 413.15, 433.15};

    scalar dH_list[9] =
        {45.054e3, 43.99e3, 43.35e3, 42.482e3,
        41.585e3, 40.657e3, 39.684e3, 38.643e3, 37.518e3};

    scalar dH;

    for (i = 1; i <= 9; i++)
    {
        X_list[i-1] = T_list[i-1];
        Y_list[i-1] = dH_list[i-1];
    }

    dH = linInterp(T, 9);
    dH = dH/18.0e-3; /* convert J/mol to J/kg */

    return dH;
    //     return 2264000.705;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::deltaHdevol
(
  label i, scalar T, scalar char_fraction
)
{
    /* returns heat of devolatilization in J/mole for T in Kelvin */
    scalar dH = 0.0;

    return dH;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::volume_of_element
(
  scalar r1, scalar r2
)
{
    /* returns volume in m3 of element */
    /* 1 = cylinder, 2 = sphere, 3 = parallelepiped */
    scalar vol = 0.0;
    scalar eta = 0.0;

    if (shape == cylinder)
    {
        eta = L1-L3;
        if (L3 == L1)
        {
            vol = 2*M_PI*(Foam::pow(r2, 3.0)-Foam::pow(r1, 3.0));
        }
        else if (L3 > L1)
        {
            vol = 2*M_PI*(Foam::pow(r2, 3.0)-Foam::pow(r2, 2.0)*eta/2)
                 -2*M_PI*(Foam::pow(r1, 3.0)-Foam::pow(r1, 2.0)*eta/2);
        }
        else if (L3 < L1)
        {
            vol =
                2*M_PI*
                    (
                        Foam::pow(r2, 3.0)+Foam::pow(r2, 2.0)*eta
                        +0.25*Foam::pow(eta, 2.0)
                    )
                 -2*M_PI*
                    (
                        Foam::pow(r1, 3.0)+Foam::pow(r1, 2.0)*eta
                        +0.25*Foam::pow(eta, 2.0)
                    );
        }
    }
    else if (shape == sphere)
    {
        vol =
            4.0*M_PI*(Foam::pow(r2, 3.0))/3.0-4.0*M_PI*(Foam::pow(r1, 3.0))/3.0;
    }
    else if (shape == parallelepiped)
    {
        vol = 8*Foam::pow(r2, 3.0)
              +4*Foam::pow(r2, 2.0)*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1))
              -
              (
                8*Foam::pow(r1, 3.0)
                +4*Foam::pow(r1, 2.0)*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1))
              );
    }
    else
    {
        printf("ERROR in volume of element function!!!\n");
    }

    return vol;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::surface_area
(
  scalar r
)
{
    /* returns surface area in m2 for the surface */
    /* 1 = cylinder, 2 = sphere, 3 = parallellepiped */
    scalar surf = 0.0;
    scalar eta = 0.0;

    if (shape == cylinder)
    {
        eta = L1-L3;
        if (L3 == L1)
        {
            surf = 6.0*M_PI*Foam::pow(r, 2.0);
        }
        else if (L3 > L1)
        {
            surf = 2.0*M_PI*(3.0*Foam::pow(r, 2.0)-r*eta);
        }
        else
        {
            surf = 2.0*M_PI*
            (
                3.0*Foam::pow(r, 2.0)+2.0*r*eta+0.25*Foam::pow(eta, 2.0)
            );
        }
    }
    else if (shape == sphere)
    {
        surf = 4.0*M_PI*(Foam::pow(r, 2.0));
    }
    else if (shape == parallelepiped)
    {
        surf = 24.0*Foam::pow(r, 2.0)
            +8.0*r*((L2-L1)+(L3-L1)+2.0*(L2-L1)*(L3-L1));
    }
    else
    {
        printf("ERROR in surface area function!!!\n");
    }

    return surf;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::radius_from_volume
(
  scalar V
)
{
    scalar r0, r1;
    scalar f = 0.0, f_prime = 0.0;
    scalar eta = L1-L3;
    label i = 0;

    r0 = 1.0*L1; /* initial guess */
    r1 = 0.5*L1;

    while
    (
        (Foam::mag((r0-r1)/r0) > MY_QUITE_SMALL)&&(i <= MY_MAXIMUM_NEWTON_ITER_)
    )
    {
        i = i+1;
        r0 = r1;

        if (shape == cylinder)
        {
            if (eta == 0)
            {
                f = (2.0*M_PI*Foam::pow(r0, 3.0)-V);
                f_prime = 6.0*M_PI*Foam::pow(r0, 2.0);
            }
            else if (L3 > L1)
            {
                f =
                    (
                        2.0*M_PI*
                            (Foam::pow(r0, 3.0)-Foam::pow(r0, 2.0)*eta/2.0)-V
                    );
                f_prime = 2.0*M_PI*(3.0*Foam::pow(r0, 2.0)-r0*eta);
            }
            else
            {
                f =
                    (
                        2.0*M_PI*
                            (
                                Foam::pow(r0, 3.0)+Foam::pow(r0, 2.0)*eta
                                +0.25*Foam::pow(eta, 2.0)
                            )
                        -V
                    );
                f_prime = 2.0*M_PI*(3.0*Foam::pow(r0, 2.0)+2.0*r0*eta);
            }
        }
        else if (shape == sphere)
        {
            f = (4.0*M_PI*Foam::pow(r0, 3.0)/3.0-V);
            f_prime = 4.0*M_PI*Foam::pow(r0, 2.0);
        }
        else if (shape == parallelepiped)
        {
            f =
                (
                    8.0*Foam::pow(r0, 3.0)
                    +4.0*Foam::pow(r0, 2.0)*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1))
                    -V
                );
            f_prime = 8.0*3.0*Foam::pow(r0, 2.0)
                +4.0*2.0*r0*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1));
        }
        else
        {
            printf("ERROR in radius from volume function!!\n");
        }

        r1 = r0-f/f_prime;
    }


    return r1;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::
radius_from_volume_element(scalar V, scalar r_inner)
{
    scalar r0, r1;
    scalar f = 0.0, f_prime = 0.0;
    scalar eta = L1-L3;
    label i = 0;

    r0 = 1.0*L1; /* initial guess */
    r1 = 0.5*L1;

    while (
        (Foam::mag((r0-r1)/r0) > MY_QUITE_SMALL)&&(i <= MY_MAXIMUM_NEWTON_ITER_))
    {
        i = i+1;
        r0 = r1;

        if (shape == cylinder)
        {
            if (eta == 0)
            {
                f = (2.0*M_PI*Foam::pow(r0, 3.0)-V);
                f_prime = 6.0*M_PI*Foam::pow(r0, 2.0);
            }
            else if (L3 > L1)
            {
                f =
                    (
                        2.0*M_PI*(Foam::pow(r0, 3.0)-Foam::pow(r0, 2.0)*eta/2.0)
                        -V
                    );
                f_prime = 2.0*M_PI*(3.0*Foam::pow(r0, 2.0)-r0*eta);
            }
            else
            {
                f =
                    (
                        2.0*M_PI*
                            (
                                Foam::pow(r0, 3.0)+Foam::pow(r0, 2.0)*eta
                                +0.25*Foam::pow(eta, 2.0)
                            )
                            -V
                    );
                f_prime = 2.0*M_PI*(3.0*Foam::pow(r0, 2.0)+2.0*r0*eta);
            }
            printf
            (
                "ERROR in radius from volume element function -- NOT PREPARED\n"
            );
        }
        else if (shape == sphere)
        {
            f =
                (
                    4.0*M_PI*Foam::pow(r0, 3.0)/3.0
                   -4.0*M_PI*Foam::pow(r_inner, 3.0)/3.0
                   -V
                );
            f_prime = 4.0*M_PI*Foam::pow(r0, 2.0);
        }
        else if (shape == parallelepiped)
        {
            f =
                (
                    8.0*Foam::pow(r0, 3.0)
                   +4.0*Foam::pow(r0, 2.0)*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1))
                   -V
                );
            f_prime = 8.0*3.0*Foam::pow(r0, 2.0)
                +4.0*2.0*r0*((L2-L1)+(L3-L1)+2*(L2-L1)*(L3-L1));
            printf
            (
                "ERROR in radius from volume element function -- NOT PREPARED\n"
            );
        }
        else
        {
            printf("ERROR in radius from volume element function!!\n");
        }

        r1 = r0-f/f_prime;
    }

    return r1;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::
temperature_gradient(scalar Ti, scalar Tj, scalar ri, scalar rj)
{
    scalar dTdr = 0.0;

    /* currently only for spheres in this version of the code */
    if (shape == sphere)
    {
        /* dTdr = (Ti-Tj)/(ri*(ri/rj-1.0)); */
        dTdr = (Ti-Tj)/(ri-rj);
    }
    else
    {
        printf("Error in temperature gradient function!\n");
    }

    return dTdr;
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::update_surface_temperatures
(
    scalar nextT_p3, scalar k_p3, scalar A, scalar r_b3, scalar r_c3,
    scalar Q_particle_conduction, scalar TInf_, scalar h, scalar emissivity
)
{
    /* currently only for spheres */
    label i = 0;
    scalar f = 0.0, f_prime = 0.0;
    scalar T0, T1;
    // scalar eta = L1-L3;

    T0 = TInf_;
    T1 = (TInf_+nextT_p3);

    while
    (
        (Foam::mag((T0-T1)/T0) > MY_QUITE_SMALL)&&(i <= MY_MAXIMUM_NEWTON_ITER_)
    )
    {
        i = i+1;
        T0 = T1;

        if (shape == sphere)
        {
            f =
                (h*A*(TInf_-T0)
                +
                view_factor*emissivity
                *Foam::constant::physicoChemical::sigma.value()
                *A*(Foam::pow(TInf_, 4.0)-Foam::pow(T0, 4.0))
                -
                k_p3*A*(T0-nextT_p3)/(r_b3-r_c3));
            f_prime =
                (
                    -h*A
                    +
                    view_factor*emissivity
                    *Foam::constant::physicoChemical::sigma.value()
                    *A*(-4.0*Foam::pow(T0, 3.0))
                    -
                    k_p3*A/(r_b3-r_c3)
                );
        }
        else
        {
            printf("ERROR in update surface temperatures function!!\n");
        }

        T1 = T0-f/f_prime;
    }

    T_send[0] = T1;
}

template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::heatRatio
(
  scalar temperature
)
{
    // scalar COCO2Ratio = 4.3*Foam::exp(-3390./temperature);
    //     return (COCO2Ratio+0.3)/(COCO2Ratio+1.);
    return 1.; //strange
}


template <class ParcelType>
Foam::scalar Foam::ReactingMultiphaseMBMParcel<ParcelType>::charFrontLocation
(
    scalar rb,
    label i,
    scalarField VCell_,
    scalarField VChar_,
    scalarField xc_,
    scalar xi_
)
{
    scalar rCharFront;

    if (i == FINE_PARTICLE_MESH_-1)
    {
        rCharFront = radius_from_volume(VChar_[i]);
    }
    else if (i == 0)
    {
        rCharFront = radius_from_volume
            (volume_of_element(0., xc_[i+1])+0.5*VCell_[i+1]+xi_*VChar_[i]);
    }
    else
    {
        rCharFront = radius_from_volume
            (volume_of_element(0., xc_[i+1])+0.5*VCell_[i+1]+VChar_[i]);
    }

    return rCharFront;
}

template <class ParcelType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::
update_boundary_temperatures
(
    scalar k_p1, scalar k_p2, scalar k_p3, scalar k_p4,
    scalar A[4], scalar rb[4], scalar rp[4],
    scalar T_p[4], scalar T_b[4],
    scalar Q_evap, scalar Q_devol, scalar Q_char_comb
)
{
    //     scalar T0 = 0.0;
    scalar T1 = 0.0;
    scalar T2 = 0.0;
    scalar T3 = 0.0;
    //     scalar eta = L1-L3;
    scalar nextT_b1;
    scalar nextT_b2;
    scalar nextT_b3;
    scalar CA, CB;

    /* currently only for spheres */

    if (shape == sphere)
    {
        /* new Tb3: */
        CA = A[boundary_3]*k_p4/(rp[layer_4]-rb[boundary_3]);
        CB = A[boundary_3]*k_p3/(rb[boundary_3]-rp[layer_3]);
        T3 = (CA*T_p[layer_4]+CB*T_p[layer_3]-Q_char_comb)/(CA+CB);
        /* new Tb2: */
        CA = A[boundary_2]*k_p3/(rp[layer_3]-rb[boundary_2]);
        CB = A[boundary_2]*k_p2/(rb[boundary_2]-rp[layer_2]);
        T2 = (CA*T_p[layer_3]+CB*T_p[layer_2]-Q_devol)/(CA+CB);
        /* new Tb1: */
        CA = A[boundary_1]*k_p2/(rp[layer_2]-rb[boundary_1]);
        CB = A[boundary_1]*k_p1/(rb[boundary_1]-rp[layer_1]);
        T1 = (CA*T_p[layer_2]+CB*T_p[layer_1]-Q_evap)/(CA+CB);
    }
    else
    {
        printf("ERROR in update boundary temperatures function!!\n");
    }

    nextT_b3 = T3;
    nextT_b2 = T2;
    nextT_b1 = T1;

    /* return temperatures in array T_send */
    T_send[0] = nextT_b1;
    T_send[1] = nextT_b2;
    T_send[2] = nextT_b3;
}
