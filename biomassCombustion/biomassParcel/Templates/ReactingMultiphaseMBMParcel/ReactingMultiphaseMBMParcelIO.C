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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ReactingMultiphaseMBMParcel<ParcelType>::propertyList_ =
    Foam::ReactingMultiphaseMBMParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ReactingMultiphaseMBMParcel<ParcelType>::sizeofFields_
(
    4*sizeof(scalar)+12*sizeof(scalarField) + 8*sizeof(List<scalarField>)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields)
{
    readCoeff(mesh.time());
    initFields();

    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
          Info << "is.format() == IOstream::ASCII " << nl;
            scalarField TpIO;
            // Info << "is.name()" << is.name() << nl;
            is >>  TOuterBoundary_ >> charFrontRadius_ >> TInf_ >> mParSum_
               >> VCell_ >> TpIO >> VChar_ >> VCharLast_ >> xc_ >> RProduct_
               >> omegaC_ >> betaR_ >> Su_ >> Sp_ >> totalDevol_
               >> charFromDevol_
               >> xi_ >> mPar_ >> mLast_ >> RCharReac_ >> QEquilibrate_ >> Rb_
               >> QReaction_ >> RDevolReac_;

            for (unsigned int i = 0; i < Tp_->max_dim; i++)
            {
              // read Tp from TpIO in to VEC* Tp_
              Tp_->ve[i] = TpIO[i];
            }
        }
        else
        {
            is.read(reinterpret_cast<char*>(&TOuterBoundary_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "ReactingMultiphaseMBMParcel<ParcelType>::ReactingMultiphaseMBMParcel"
        "("
            "const polyMesh&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    bool valid = c.size();

    ParcelType::readFields(c, compModel);

    IOField<scalarField> VCell(c.fieldIOobject("VCell", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, VCell);

    IOField<scalarField> xi0(c.fieldIOobject("xi0", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, xi0);

    IOField<scalarField> xi1(c.fieldIOobject("xi1", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, xi1);

    IOField<scalarField> xi2(c.fieldIOobject("xi2", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, xi2);

    IOField<scalarField> xi3(c.fieldIOobject("xi3", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, xi3);

    IOField<scalarField> mPar0(c.fieldIOobject("mPar0", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mPar0);

    IOField<scalarField> mPar1(c.fieldIOobject("mPar1", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mPar1);

    IOField<scalarField> mPar2(c.fieldIOobject("mPar2", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mPar2);

    IOField<scalarField> mPar3(c.fieldIOobject("mPar3", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mPar3);

    IOField<scalarField> mLast0(c.fieldIOobject("mLast0", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mLast0);

    IOField<scalarField> mLast1(c.fieldIOobject("mLast1", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mLast1);

    IOField<scalarField> mLast2(c.fieldIOobject("mLast2", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mLast2);

    IOField<scalarField> mLast3(c.fieldIOobject("mLast3", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mLast3);

    IOField<scalarField> Tp(c.fieldIOobject("Tp", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, Tp);

    IOField<scalarField> VChar(c.fieldIOobject("VChar", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, VChar);

    IOField<scalarField> RbDrf(c.fieldIOobject("RbDrf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RbDrf);

    IOField<scalarField> RbDvf(c.fieldIOobject("RbDvf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RbDvf);

    IOField<scalarField> RbCf(c.fieldIOobject("RbCf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RbCf);

    IOField<scalarField> RbCp(c.fieldIOobject("RbCp", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RbCp);

    IOField<scalarField> QEquilibrateDrf(c.fieldIOobject("QEquilibrateDrf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QEquilibrateDrf);

    IOField<scalarField> QEquilibrateDvf(c.fieldIOobject("QEquilibrateDvf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QEquilibrateDvf);

    IOField<scalarField> QEquilibrateCf(c.fieldIOobject("QEquilibrateCf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QEquilibrateCf);

    IOField<scalarField> QReactionEv(c.fieldIOobject("QReactionEv", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QReactionEv);

    IOField<scalarField> QReactionDv(c.fieldIOobject("QReactionDv", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QReactionDv);

    IOField<scalarField> QReactionCC(c.fieldIOobject("QReactionCC", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, QReactionCC);

    IOField<scalarField> VCharLast(c.fieldIOobject("VCharLast", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, VCharLast);

    IOField<scalarField> xc(c.fieldIOobject("xc", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, xc);

    IOField<scalarField> RProduct(c.fieldIOobject("RProduct", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RProduct);

    IOField<scalarField> RDevolReac0(c.fieldIOobject("RDevolReac0", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RDevolReac0);

    IOField<scalarField> RDevolReac1(c.fieldIOobject("RDevolReac1", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RDevolReac1);

    IOField<scalarField> RDevolReac2(c.fieldIOobject("RDevolReac2", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RDevolReac2);

    IOField<scalarField> RCharReac0(c.fieldIOobject("RCharReac0", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RCharReac0);

    IOField<scalarField> RCharReac1(c.fieldIOobject("RCharReac1", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RCharReac1);

    IOField<scalarField> RCharReac2(c.fieldIOobject("RCharReac2", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RCharReac2);

    IOField<scalarField> RCharReac3(c.fieldIOobject("RCharReac3", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, RCharReac3);

    IOField<scalarField> omegaC(c.fieldIOobject("omegaC", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, omegaC);

    IOField<scalarField> betaR(c.fieldIOobject("betaR", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, betaR);

    IOField<scalarField> Su(c.fieldIOobject("Su", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, Su);

    IOField<scalarField> Sp(c.fieldIOobject("Sp", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, Sp);

    IOField<scalarField> totalDevol(c.fieldIOobject("totalDevol", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, totalDevol);

    IOField<scalarField> charFromDevol(c.fieldIOobject("charFromDevol", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, charFromDevol);


    IOField<scalar> TOuterBoundary(c.fieldIOobject("TOuterBoundary", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, TOuterBoundary);

    IOField<scalar> charFrontRadius(c.fieldIOobject("charFrontRadius", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, charFrontRadius);

    IOField<scalar> TInf(c.fieldIOobject("TInf", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, TInf);

    IOField<scalar> mParSum(c.fieldIOobject("mParSum", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, mParSum);

    IOField<scalar> tReal(c.fieldIOobject("tReal", IOobject::READ_IF_PRESENT), valid);
    c.checkFieldIOobject(c, tReal);

    label i = 0;
    forAllIter(typename Cloud<ReactingMultiphaseMBMParcel<ParcelType> >, c, iter)
    {
        ReactingMultiphaseMBMParcel<ParcelType>& p = iter();

        p.VCell_ = VCell[i];

        p.xi_[layer_1] = xi0[i];
        p.xi_[layer_2] = xi1[i];
        p.xi_[layer_3] = xi2[i];
        p.xi_[layer_4] = xi3[i];

        p.mPar_[layer_1] = mPar0[i];
        p.mPar_[layer_2] = mPar1[i];
        p.mPar_[layer_3] = mPar2[i];
        p.mPar_[layer_4] = mPar3[i];

        p.mLast_[layer_1] = mLast0[i];
        p.mLast_[layer_2] = mLast1[i];
        p.mLast_[layer_3] = mLast2[i];
        p.mLast_[layer_4] = mLast3[i];

        for (unsigned int j = 0; j < p.Tp_->max_dim; j++)
        {
            p.Tp_->ve[j] = Tp[i][j];
        }

        p.VChar_ = VChar[i];

        p.Rb_[drying_front] = RbDrf[i];
        p.Rb_[devol_front] = RbDvf[i];
        p.Rb_[char_front] = RbCf[i];
        p.Rb_[char_produced] = RbCp[i];

        p.VCharLast_ = VCharLast[i];
        p.xc_ = xc[i];
        p.RProduct_ = RProduct[i];

        p.QEquilibrate_[drying_front] = QEquilibrateDrf[i];
        p.QEquilibrate_[devol_front] = QEquilibrateDvf[i];
        p.QEquilibrate_[char_front] = QEquilibrateCf[i];

        p.QReaction_[evaporation] = QReactionEv[i];
        p.QReaction_[devolatilization] = QReactionDv[i];
        p.QReaction_[char_combustion] = QReactionCC[i];

        p.RDevolReac_[0] = RDevolReac0[i];
        p.RDevolReac_[1] = RDevolReac1[i];
        p.RDevolReac_[2] = RDevolReac2[i];

        p.TOuterBoundary_ = TOuterBoundary[i];

        p.RCharReac_[0] = RCharReac0[i];
        p.RCharReac_[1] = RCharReac1[i];
        p.RCharReac_[2] = RCharReac2[i];
        p.RCharReac_[3] = RCharReac3[i];

        p.omegaC_ = omegaC[i];

        p.betaR_ = betaR[i];

        p.Su_ = Su[i];
        p.Sp_ = Sp[i];

        p.charFrontRadius_ = charFrontRadius[i];
        p.TInf_ = TInf[i];

        p.mParSum_ = mParSum[i];

        p.totalDevol_ = totalDevol[i];

        p.charFromDevol_ = charFromDevol[i];

        p.tReal_ = tReal[i];

        i++;
    }

}


template<class ParcelType>
template<class CloudType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ReactingMultiphaseMBMParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();
    {
        IOField<scalarField> VCell(c.fieldIOobject("VCell", IOobject::NO_READ), np);
        IOField<scalarField> xi0(c.fieldIOobject("xi0", IOobject::NO_READ), np);
        IOField<scalarField> xi1(c.fieldIOobject("xi1", IOobject::NO_READ), np);
        IOField<scalarField> xi2(c.fieldIOobject("xi2", IOobject::NO_READ), np);
        IOField<scalarField> xi3(c.fieldIOobject("xi3", IOobject::NO_READ), np);
        IOField<scalarField> mPar0(c.fieldIOobject("mPar0", IOobject::NO_READ), np);
        IOField<scalarField> mPar1(c.fieldIOobject("mPar1", IOobject::NO_READ), np);
        IOField<scalarField> mPar2(c.fieldIOobject("mPar2", IOobject::NO_READ), np);
        IOField<scalarField> mPar3(c.fieldIOobject("mPar3", IOobject::NO_READ), np);
        IOField<scalarField> mLast0(c.fieldIOobject("mLast0", IOobject::NO_READ), np);
        IOField<scalarField> mLast1(c.fieldIOobject("mLast1", IOobject::NO_READ), np);
        IOField<scalarField> mLast2(c.fieldIOobject("mLast2", IOobject::NO_READ), np);
        IOField<scalarField> mLast3(c.fieldIOobject("mLast3", IOobject::NO_READ), np);
        IOField<scalarField> Tp(c.fieldIOobject("Tp", IOobject::NO_READ), np);
        IOField<scalarField> VChar(c.fieldIOobject("VChar", IOobject::NO_READ), np);
        IOField<scalarField> RbDrf(c.fieldIOobject("RbDrf", IOobject::NO_READ), np);
        IOField<scalarField> RbDvf(c.fieldIOobject("RbDvf", IOobject::NO_READ), np);
        IOField<scalarField> RbCf(c.fieldIOobject("RbCf", IOobject::NO_READ), np);
        IOField<scalarField> RbCp(c.fieldIOobject("RbCp", IOobject::NO_READ), np);
        IOField<scalarField> QEquilibrateDrf(c.fieldIOobject("QEquilibrateDrf", IOobject::NO_READ), np);
        IOField<scalarField> QEquilibrateDvf(c.fieldIOobject("QEquilibrateDvf", IOobject::NO_READ), np);
        IOField<scalarField> QEquilibrateCf(c.fieldIOobject("QEquilibrateCf", IOobject::NO_READ), np);
        IOField<scalarField> QReactionEv(c.fieldIOobject("QReactionEv", IOobject::NO_READ), np);
        IOField<scalarField> QReactionDv(c.fieldIOobject("QReactionDv", IOobject::NO_READ), np);
        IOField<scalarField> QReactionCC(c.fieldIOobject("QReactionCC", IOobject::NO_READ), np);
        IOField<scalarField> RDevolReac0(c.fieldIOobject("RDevolReac0", IOobject::NO_READ), np);
        IOField<scalarField> RDevolReac1(c.fieldIOobject("RDevolReac1", IOobject::NO_READ), np);
        IOField<scalarField> RDevolReac2(c.fieldIOobject("RDevolReac2", IOobject::NO_READ), np);
        IOField<scalarField> RCharReac0(c.fieldIOobject("RCharReac0", IOobject::NO_READ), np);
        IOField<scalarField> RCharReac1(c.fieldIOobject("RCharReac1", IOobject::NO_READ), np);
        IOField<scalarField> RCharReac2(c.fieldIOobject("RCharReac2", IOobject::NO_READ), np);
        IOField<scalarField> RCharReac3(c.fieldIOobject("RCharReac3", IOobject::NO_READ), np);
        IOField<scalarField> omegaC(c.fieldIOobject("omegaC", IOobject::NO_READ), np);
        IOField<scalarField> betaR(c.fieldIOobject("betaR", IOobject::NO_READ), np);
        IOField<scalarField> Su(c.fieldIOobject("Su", IOobject::NO_READ), np);
        IOField<scalarField> Sp(c.fieldIOobject("Sp", IOobject::NO_READ), np);
        IOField<scalarField> VCharLast(c.fieldIOobject("VCharLast", IOobject::NO_READ), np);
        IOField<scalarField> xc(c.fieldIOobject("xc", IOobject::NO_READ), np);
        IOField<scalarField> RProduct(c.fieldIOobject("RProduct", IOobject::NO_READ), np);
        IOField<scalarField> totalDevol(c.fieldIOobject("totalDevol", IOobject::NO_READ), np);
        IOField<scalarField> charFromDevol(c.fieldIOobject("charFromDevol", IOobject::NO_READ), np);

        IOField<scalar> TOuterBoundary(c.fieldIOobject("TOuterBoundary", IOobject::NO_READ), np);
        IOField<scalar> charFrontRadius(c.fieldIOobject("charFrontRadius", IOobject::NO_READ), np);
        IOField<scalar> TInf(c.fieldIOobject("TInf", IOobject::NO_READ), np);
        IOField<scalar> mParSum(c.fieldIOobject("mParSum", IOobject::NO_READ), np);
        IOField<scalar> tReal(c.fieldIOobject("tReal", IOobject::NO_READ), np);

        label i = 0;
        forAllConstIter(typename Cloud<ReactingMultiphaseMBMParcel<ParcelType> >, c, iter)
        {
            const ReactingMultiphaseMBMParcel<ParcelType>& p = iter();
            VCell[i] = p.VCell_;
            xi0[i] = p.xi_[layer_1];
            xi1[i] = p.xi_[layer_2];
            xi2[i] = p.xi_[layer_3];
            xi3[i] = p.xi_[layer_4];
            mPar0[i] = p.mPar_[layer_1];
            mPar1[i] = p.mPar_[layer_2];
            mPar2[i] = p.mPar_[layer_3];
            mPar3[i] = p.mPar_[layer_4];
            mLast0[i] = p.mLast_[layer_1];
            mLast1[i] = p.mLast_[layer_2];
            mLast2[i] = p.mLast_[layer_3];
            mLast3[i] = p.mLast_[layer_4];

            scalarField TpSave;
            for (unsigned int j = 0; j < p.Tp_->max_dim; j++)
            {
                TpSave.append(p.Tp_->ve[j]);
            }

            Tp[i] = TpSave;
            VChar[i] = p.VChar_;

            RbDrf[i] = p.Rb_[drying_front];
            RbDvf[i] = p.Rb_[devol_front];
            RbCf[i] = p.Rb_[char_front];
            RbCp[i] = p.Rb_[char_produced];

            QEquilibrateDrf[i] = p.QEquilibrate_[drying_front];
            QEquilibrateDvf[i] = p.QEquilibrate_[devol_front];
            QEquilibrateCf[i] = p.QEquilibrate_[char_front];
            //
            QReactionEv[i] = p.QReaction_[evaporation];
            QReactionDv[i] = p.QReaction_[devolatilization];
            QReactionCC[i] = p.QReaction_[char_combustion];


            VCharLast[i] = p.VCharLast_;
            xc[i] = p.xc_;
            RProduct[i] = p.RProduct_;

            RDevolReac0[i] = p.RDevolReac_[0];
            RDevolReac1[i] = p.RDevolReac_[1];
            RDevolReac2[i] = p.RDevolReac_[2];

            TOuterBoundary[i] = p.TOuterBoundary_;

            RCharReac0[i] = p.RCharReac_[0];
            RCharReac1[i] = p.RCharReac_[1];
            RCharReac2[i] = p.RCharReac_[2];
            RCharReac3[i] = p.RCharReac_[3];

            omegaC[i] = p.omegaC_;

            betaR[i] = p.betaR_;

            Su[i] = p.Su_;
            Sp[i] = p.Sp_;

            charFrontRadius[i] = p.charFrontRadius_;
            TInf[i] = p.TInf_;

            mParSum[i] = p.mParSum_;

            totalDevol[i] = p.totalDevol_;

            charFromDevol[i] = p.charFromDevol_;

            tReal[i] = p.tReal_;

            i++;
        }

        VCell.write(np > 0);
        xi0.write(np > 0);
        xi1.write(np > 0);
        xi2.write(np > 0);
        xi3.write(np > 0);
        mPar0.write(np > 0);
        mPar1.write(np > 0);
        mPar2.write(np > 0);
        mPar3.write(np > 0);
        mLast0.write(np > 0);
        mLast1.write(np > 0);
        mLast2.write(np > 0);
        mLast3.write(np > 0);
        Tp.write(np > 0);
        VChar.write(np > 0);
        RbDrf.write(np > 0);
        RbDvf.write(np > 0);
        RbCf.write(np > 0);
        RbCp.write(np > 0);
        QEquilibrateDrf.write(np > 0);
        QEquilibrateDvf.write(np > 0);
        QEquilibrateCf.write(np > 0);
        QReactionEv.write(np > 0);
        QReactionDv.write(np > 0);
        QReactionCC.write(np > 0);
        VCharLast.write(np > 0);
        xc.write(np > 0);
        RProduct.write(np > 0);
        RDevolReac0.write(np > 0);
        RDevolReac1.write(np > 0);
        RDevolReac2.write(np > 0);
        TOuterBoundary.write(np > 0);
        RCharReac0.write(np > 0);
        RCharReac1.write(np > 0);
        RCharReac2.write(np > 0);
        RCharReac3.write(np > 0);
        omegaC.write(np > 0);
        betaR.write(np > 0);
        Su.write(np > 0);
        Sp.write(np > 0);
        charFrontRadius.write(np > 0);
        TInf.write(np > 0);
        mParSum.write(np > 0);
        totalDevol.write(np > 0);
        charFromDevol.write(np > 0);
        tReal.write(np > 0);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ReactingMultiphaseMBMParcel<ParcelType>& p
)
{

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.TOuterBoundary()
            << token::SPACE << p.charFrontRadius()
            << token::SPACE << p.TInf()
            << token::SPACE << p.mParSum()
            << token::SPACE << p.VCell()
            << token::SPACE << p.Tp()
            << token::SPACE << p.VChar()
            << token::SPACE << p.VCharLast()
            << token::SPACE << p.xc()
            << token::SPACE << p.RProduct()
            << token::SPACE << p.omegaC()
            << token::SPACE << p.betaR()
            << token::SPACE << p.Su()
            << token::SPACE << p.Sp()
            << token::SPACE << p.totalDevol()
            << token::SPACE << p.charFromDevol()
            << token::SPACE << p.xi()
            << token::SPACE << p.mPar()
            << token::SPACE << p.mLast()
            << token::SPACE << p.RCharReac()
            << token::SPACE << p.QEquilibrate()
            << token::SPACE << p.Rb()
            << token::SPACE << p.QReaction()
            << token::SPACE << p.RDevolReac();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.TOuterBoundary_),
            ReactingMultiphaseMBMParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const ReactingMultiphaseMBMParcel<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
