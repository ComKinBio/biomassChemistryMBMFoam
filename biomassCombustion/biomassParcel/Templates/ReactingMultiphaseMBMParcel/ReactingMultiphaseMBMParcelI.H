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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties()
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict)
{}


template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti)
{
    readCoeff(mesh.time());
    initFields();
}


template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli)
{
    readCoeff(mesh.time());
    initFields();
}


template<class ParcelType>
inline Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
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
)
:
    ParcelType
    (
        mesh,
        coordinates,
        celli,
        tetFacei,
        tetPti,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    )
{
    readCoeff(mesh.time());
    initFields();
}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //
template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseMBMParcel<ParcelType>::TOuterBoundary() const
{
      return TOuterBoundary_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseMBMParcel<ParcelType>::charFrontRadius() const
{
      return charFrontRadius_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseMBMParcel<ParcelType>::TInf() const
{
      return TInf_;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mParSum() const
{
      return mParSum_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VCell() const
{
      return VCell_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Tp() const
{
    scalarField TpSave;
    for (unsigned int j = 0; j < Tp_->max_dim; j++)
    {
        TpSave.append(Tp_->ve[j]);
    }

    return TpSave;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VChar() const
{
      return VChar_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VCharLast() const
{
      return VCharLast_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::xc() const
{
      return xc_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RProduct() const
{
      return RProduct_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::omegaC() const
{
      return omegaC_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::betaR() const
{
      return betaR_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Su() const
{
      return Su_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Sp() const
{
      return Sp_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::totalDevol() const
{
      return totalDevol_;
}

template<class ParcelType>
inline Foam::scalarField
Foam::ReactingMultiphaseMBMParcel<ParcelType>::charFromDevol() const
{
      return charFromDevol_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::xi() const
{
      return xi_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mPar() const
{
      return mPar_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mLast() const
{
      return mLast_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RCharReac() const
{
      return RCharReac_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::QEquilibrate() const
{
      return QEquilibrate_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Rb() const
{
      return Rb_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::QReaction() const
{
      return QReaction_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RDevolReac() const
{
      return RDevolReac_;
}

// some initial constant value
template<class ParcelType>
inline Foam::scalar
Foam::ReactingMultiphaseMBMParcel<ParcelType>::MOISTURE_FRACTION_DB() const
{
      return MOISTURE_FRACTION_DB_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Edit
template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::TOuterBoundary()
{
      return TOuterBoundary_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::charFrontRadius()
{
      return charFrontRadius_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::TInf()
{
      return TInf_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mParSum()
{
      return mParSum_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VCell()
{
    return VCell_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Tp()
{
    scalarField TpSave;
    for (unsigned int j = 0; j < Tp_->max_dim; j++)
    {
        TpSave.append(Tp_->ve[j]);
    }

    return TpSave;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VChar()
{
      return VChar_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::VCharLast()
{
      return VCharLast_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::xc()
{
      return xc_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RProduct()
{
      return RProduct_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::omegaC()
{
      return omegaC_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::betaR()
{
      return betaR_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Su()
{
      return Su_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Sp()
{
      return Sp_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::totalDevol()
{
      return totalDevol_;
}

template<class ParcelType>
inline Foam::scalarField&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::charFromDevol()
{
      return charFromDevol_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::xi()
{
    return xi_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mPar()
{
    return mPar_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::mLast()
{
    return mLast_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RCharReac()
{
      return RCharReac_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::QEquilibrate()
{
      return QEquilibrate_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::Rb()
{
      return Rb_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::QReaction()
{
      return QReaction_;
}

template<class ParcelType>
inline Foam::List<Foam::scalarField>&
Foam::ReactingMultiphaseMBMParcel<ParcelType>::RDevolReac()
{
      return RDevolReac_;
}
// ************************************************************************* //
