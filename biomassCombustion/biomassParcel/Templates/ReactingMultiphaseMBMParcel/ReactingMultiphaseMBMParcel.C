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

#include "ReactingMultiphaseMBMParcel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


// * * * * * * * * * * * *  Helper Member Functions * * * * * * * * * * * * * //
#include "externalFunctions.H"
// * * * * * * * * * * * * Main clac() Function * * * * * * * * * * * * * * * //
#include "MBMParcelMainCalcFunction.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
(
    const ReactingMultiphaseMBMParcel<ParcelType>& p
)
:
    ParcelType(p)
{}


template<class ParcelType>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
(
    const ReactingMultiphaseMBMParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseMBMParcelIO.C"

// ************************************************************************* //
