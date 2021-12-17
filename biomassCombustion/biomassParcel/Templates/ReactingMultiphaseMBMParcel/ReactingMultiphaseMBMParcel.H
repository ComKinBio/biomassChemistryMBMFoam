/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::ReactingMultiphaseMBMParcel

Description
    Multiphase variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    ReactingMultiphaseMBMParcelI.H
    ReactingMultiphaseMBMParcel.C
    ReactingMultiphaseMBMParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseMBMParcel_H
#define ReactingMultiphaseMBMParcel_H

// include Eigen solver
#include <eigen-3.4.0/Eigen/Dense>
#include <eigen-3.4.0/Eigen/Sparse>
#include <eigen-3.4.0/Eigen/IterativeLinearSolvers>

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

#ifdef __cplusplus
extern "C"
{
  #endif

  #include "matrix2.h"

  #ifdef __cplusplus
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingMultiphaseMBMParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingMultiphaseMBMParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                 Class ReactingMultiphaseMBMParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseMBMParcel
:
    public ParcelType
{
    // Private Data

        /* model parameters */
        #define MY_QUITE_SMALL 1.0e-3
        #define MY_VERY_SMALL 1.0e-6
        #define MY_EXTREMELY_SMALL 1.0e-20
        #define MY_QUITE_LARGE 1.0e3
        #define MY_VERY_LARGE 1.0e6
        // #define T_P_INIT 300.0
        // #define emissivity 0.85

        /* molar masses */
        const scalar M_C = 12.011e-3;           /* kg/mol */
        const scalar M_O = 15.999e-3;           /* kg/mol */
        const scalar M_H = 1.008e-3;            /* kg/mol */
        const scalar M_H2O = 18.01528e-3;       /* kg/mol */
        const scalar M_O2 = (2.0 * M_O);        /* kg/mol */
        const scalar M_H2 = (2.0 * M_H);        /* kg/mol */
        const scalar M_N2 = (2.0 * 14.007e-3);  /* kg/mol */
        const scalar M_CO2 = (M_C + M_O2);      /* kg/mol */
        const scalar M_CO = (M_C + M_O);        /* kg/mol */
        const scalar M_CH4 = (M_C + 4.0 * M_H); /* kg/mol */

        /* model indices */
        enum
        {
          layer_1,
          layer_2,
          layer_3,
          layer_4
        };
        enum
        {
          boundary_1,
          boundary_2,
          boundary_3,
          boundary_4
        };
        enum
        {
          drying_front, //drf
          devol_front, //dvf
          char_front, //cf
          char_produced //cp
        };
        enum
        {
          charring,
          tarring,
          gassing,
        };
        enum
        {
          evaporation, //ev
          devolatilization, //dv
          char_combustion //cc
        };
        enum
        {
          cylinder,
          sphere,
          parallelepiped
        };

        /* property indices */
        enum
        {
          index_CO2,
          index_H2,
          index_CO,
          index_H2O,
          index_CH4,
          index_C2H4,
          index_C10H8,
          index_C,
          index_C2H6,
          index_C3H8,
          index_C6H6O,
          index_volatiles
        };

        enum
        {
          H2O_index,
          O2_index,
          CH4_index,
          H2_index,
          CO_index,
          CO2_index,
          N2_index
        };

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        /* static interpolation lists */
        scalar X_list[50];
        scalar Y_list[50];

        /* static message passing */
        scalar T_send[3];

        scalar view_factor = 1.;
        scalar m_init;
        scalar r_init;
        scalar L1, L2, L3;
        scalar particle_volume0;

        /* Redux statics */
        label material, shape;

        /* read some parameters from file */
        label VERBOSITY_, SAVE_EVERY_, FINE_PARTICLE_MESH_, N_SPECIES_;

        scalar TIME_DT_, PARTICLE_RADIUS_, T_REACTOR_,
        T_RADIATION_,T_P_INIT_, MASS_FRAC_O2_, DENSITY_INITIAL_, MOISTURE_FRACTION_DB_,
        MOISTURE_FRACTION_WB, FINE_PARTICLE_MESH_DOUBLE_,
        DRYING_SHRINKAGE_VOLUME_, ASH_MASS_FRACTION_DB, K_P_WET_, K_P_DRY_,
        K_P_CHAR_,ASH_MASS_FRACTION_, C_FRACTION_, H_FRACTION_, O_FRACTION_,
        CHARCONVERSION_SHRINKAGE_VOLUME_,DEVOL_SHRINKAGE_VOLUME_,
        P_ATM_, MASS_FRAC_H2O_, MASS_FRAC_H2_, MASS_FRAC_CO2_,
        MY_MAXIMUM_NEWTON_ITER_; /*DRY_DENSITY*/;

        /* drying kinetics parameters */
        scalar DryKinA_, DryKinE_;

        /* pyrolysis kinetics parameters */
        scalar DevKinA1_, DevKinE1_, DevKinA2_, DevKinE2_, DevKinA3_, DevKinE3_;

        // count interatrion number for print out results
        label iter_count, output_every;

        scalar tPar; //,t_stop;
        scalar aw, ae, ap, ap0, Su__b, Sp__b;
        scalar TOuterBoundary_;
        scalar Tf, Af;
        scalar rb, TInf_;
        scalar S, Q, q, a, heatHKp, e, delta_0, delta_1;
        scalar rho_g, mu, cp_g, k_g, vel_mag, h;
        scalar RCharReac_1S;

        scalar total_mass;

        scalar R_drySourceimins1High, R_drySourceiHigh, R_devoSourceimins1High,
        R_devoSourceiHigh, R_devoSourceimins1HighChar, R_devoSourceiHighChar,
        R_drySourceimins1Low, R_drySourceiLow, R_devoSourceimins1Low,
        R_devoSourceiLow, R_devoSourceimins1LowChar, R_devoSourceiLowChar,
        density_drywood;

        scalar D_g, D_gAB;

        scalar hmi = 0.0, hmia, beta_d, alpha = 0.0, conc_O2, conc_O2_correct,
        O2_diff_suface = 0.0, omegaC_S;

        scalar volume_sum;

        scalar mParSum_;

        scalar high = 0.0, low = 0.0;

        scalar R_O2_mole_need = 0.;

        scalar charFrontRadius_;

        scalar tReal_ = 0.0;

        label combustion_flag, combustion_layer;

        List<scalarField> xi_,QReaction_,Rb_,RCharReac_;
        List<scalarField> mPar_,mLast_;
        List<scalarField> QEquilibrate_,RDevolReac_;
        List<scalarField> dmdt;

        scalarField VCell_;
        scalarField VCharLast_,VChar_,xc_,Sp_,Su_,xif_,omegaC_;
        scalarField sStar_,dSdT_star;
        scalarField charFromDevol_,totalDevol_;
        scalarField RProduct_;
        scalarField charYield_;
        scalarField betaR_;
        scalarField R_O2_mole_need_list;
        scalarField ThermalTimeScale;

        /* static library variables */
        scalarField LHV,M;

        // meschach members
        VEC *b, *Tp_;
        MAT *A, *LU;
        PERM *pivot;


public:

    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data

    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;


private:

    // Private Member Functions
    // initialize
    void ReadLHVM(scalar YC); //YC is the char Fraction

    void initFields();

    void readMeshInput(const dictionary &particleInfoDict_);

    void readDevKin(const dictionary &particleInfoDict_);

    // call both readMeshInput and readDevKin
    void readCoeff(const Time &time);

    // Helper functions(Functions added from the mesh code)
    scalar c_p(label layer, scalar T);

    scalar k_p(label layer, scalar T);

    scalar rhoPar(label layer, scalar T);

    scalar gas_cp(label index, scalar T);

    scalar deltaHvap(scalar T);

    scalar deltaHdevol(label i, scalar T, scalar char_fraction);

    scalar volume_of_element(scalar r1, scalar r2);

    scalar surface_area(scalar r);

    scalar radius_from_volume(scalar V);

    scalar radius_from_volume_element(scalar V, scalar r_inner);

    scalar temperature_gradient(scalar Ti, scalar Tj, scalar ri, scalar rj);

    void update_surface_temperatures
    (
        scalar nextT_p3, scalar k_p3, scalar A, scalar r_b3, scalar r_c3,
        scalar Q_particle_conduction, scalar TInf_, scalar h,
        scalar emissivity
    );

    scalar heatRatio(scalar temperature);

    scalar charFrontLocation
    (
        scalar rb,
        label i,
        scalarField V_cell,
        scalarField VChar_,
        scalarField xc_,
        scalar xi_
    );

    void update_boundary_temperatures
    (
        scalar k_p1, scalar k_p2, scalar k_p3, scalar k_p4,
        scalar A[4], scalar rb[4], scalar rp[4],
        scalar T_p[4], scalar T_b[4],
        scalar Q_evap, scalar Q_devol, scalar Q_char_comb
    );

    //Do linear interpolation
    scalar linInterp(scalar X, label n_elements_in_list);

    // Eigen CG solver
    void EigenCGFunc(MAT* A, VEC* T, VEC* b);

protected:

    // Protected data

        // Parcel properties


    // Protected Member Functions



public:

    // Static Data Members

        //- Runtime type information
        TypeName("ReactingMultiphaseMBMParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " TOuterBoundary"
          + " charFrontRadius"
          + " TInf"
          + " mParSum"
          + " VCell"
          + " Tp"
          + " VChar"
          + " VCharLast"
          + " xc"
          + " RProduct"
          + " omegaC"
          + " betaR"
          + " Su"
          + " Sp"
          + " totalDevol"
          + " charFromDevol"
          + " xi"
          + " mPar"
          + " mLast"
          + " RCharReac"
          + " QEquilibrate"
          + " Rb"
          + " QReaction"
          + " RDevolReac"
        );

    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline ReactingMultiphaseMBMParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingMultiphaseMBMParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline ReactingMultiphaseMBMParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalarField& Y0,
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingMultiphaseMBMParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingMultiphaseMBMParcel(const ReactingMultiphaseMBMParcel& p);

        //- Construct as a copy
        ReactingMultiphaseMBMParcel
        (
            const ReactingMultiphaseMBMParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingMultiphaseMBMParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ReactingMultiphaseMBMParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingMultiphaseMBMParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<ReactingMultiphaseMBMParcel<ParcelType>>
                (
                    new ReactingMultiphaseMBMParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline scalar TOuterBoundary() const;
            inline scalar charFrontRadius() const;
            inline scalar TInf() const;
            inline scalar mParSum() const;

            inline scalarField VCell() const;
            inline scalarField Tp() const;
            inline scalarField VChar() const;
            inline scalarField VCharLast() const;
            inline scalarField xc() const;
            inline scalarField RProduct() const;
            inline scalarField omegaC() const;
            inline scalarField betaR() const;
            inline scalarField Su() const;
            inline scalarField Sp() const;
            inline scalarField totalDevol() const;
            inline scalarField charFromDevol() const;

            inline List<scalarField> xi() const;
            inline List<scalarField> mPar() const;
            inline List<scalarField> mLast() const;
            inline List<scalarField> RCharReac() const;
            inline List<scalarField> QEquilibrate() const;
            inline List<scalarField> Rb() const;
            inline List<scalarField> QReaction() const;
            inline List<scalarField> RDevolReac() const;

            // some initial constant value
            inline scalar MOISTURE_FRACTION_DB() const;



        // Edit

            //- Return access to front tempearture
            inline scalar& TOuterBoundary();
            inline scalar& charFrontRadius();
            inline scalar& TInf();
            inline scalar& mParSum();

            inline scalarField& VCell();
            inline scalarField& Tp();
            inline scalarField& VChar();
            inline scalarField& VCharLast();
            inline scalarField& xc();
            inline scalarField& RProduct();
            inline scalarField& omegaC();
            inline scalarField& betaR();
            inline scalarField& Su();
            inline scalarField& Sp();
            inline scalarField& totalDevol();
            inline scalarField& charFromDevol();

            inline List<scalarField>& xi();
            inline List<scalarField>& mPar();
            inline List<scalarField>& mLast();
            inline List<scalarField>& RCharReac();
            inline List<scalarField>& QEquilibrate();
            inline List<scalarField>& Rb();
            inline List<scalarField>& QReaction();
            inline List<scalarField>& RDevolReac();

            // some initial constant value
            inline scalar& MOISTURE_FRACTION_DB();


        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingMultiphaseMBMParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseMBMParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseMBMParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
