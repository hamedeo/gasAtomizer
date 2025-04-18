/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "SolidifyingSprayParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::SolidifyingSprayParcel<ParcelType>::propertyList_ =
    Foam::SolidifyingSprayParcel<ParcelType>::propertyList();


template<class ParcelType>
const std::size_t Foam::SolidifyingSprayParcel<ParcelType>::sizeofFields
(
    sizeof(SolidifyingSprayParcel<ParcelType>) - sizeof(ParcelType)
    // 0    // added (should this line be added?)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    d0_(0.0),
    position0_(Zero),
    sigma_(0.0),
    mu_(0.0),
    liquidCore_(0.0),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(GREAT),
    user_(0.0),
    YGas_(0), // added
    YLiquid_(0), // added
    YSolid_(0) // added
{
    Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> d0_
                >> position0_
                >> sigma_
                >> mu_
                >> liquidCore_
                >> KHindex_
                >> y_
                >> yDot_
                >> tc_
                >> ms_
                >> injector_
                >> tMom_
                >> user_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &d0_);
            readRawScalar(is, position0_.data(), vector::nComponents);
            readRawScalar(is, &sigma_);
            readRawScalar(is, &mu_);
            readRawScalar(is, &liquidCore_);
            readRawScalar(is, &KHindex_);
            readRawScalar(is, &y_);
            readRawScalar(is, &yDot_);
            readRawScalar(is, &tc_);
            readRawScalar(is, &ms_);
            readRawScalar(is, &injector_);
            readRawScalar(is, &tMom_);
            readRawScalar(is, &user_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d0_), sizeofFields);
        }
        DynamicList<scalar> Yg;   // added
        DynamicList<scalar> Yl;   // added
        DynamicList<scalar> Ys;   // added

        is >> Yg >> Yl >> Ys;   // added

        YGas_.transfer(Yg);   // added
        YLiquid_.transfer(Yl);   // added
        YSolid_.transfer(Ys);   // added
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::readFields(CloudType& c)
{
    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SolidifyingSprayParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    const bool valid = c.size();

    ParcelType::readFields(c, compModel);

    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, d0);

    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, position0);

    IOField<scalar> sigma(c.fieldIOobject("sigma", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, sigma);

    IOField<scalar> mu(c.fieldIOobject("mu", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, mu);

    IOField<scalar> liquidCore
    (
        c.fieldIOobject("liquidCore", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, liquidCore);

    IOField<scalar> KHindex
    (
        c.fieldIOobject("KHindex", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, KHindex);

    IOField<scalar> y
    (
        c.fieldIOobject("y", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, y);

    IOField<scalar> yDot
    (
        c.fieldIOobject("yDot", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, yDot);

    IOField<scalar> tc
    (
        c.fieldIOobject("tc", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, tc);

    IOField<scalar> ms
    (
        c.fieldIOobject("ms", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, ms);

    IOField<scalar> injector
    (
        c.fieldIOobject("injector", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, injector);

    IOField<scalar> tMom
    (
        c.fieldIOobject("tMom", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, tMom);

    IOField<scalar> user
    (
        c.fieldIOobject("user", IOobject::MUST_READ),
        valid
    );
    c.checkFieldIOobject(c, user);

    label i = 0;
    for (SolidifyingSprayParcel<ParcelType>& p : c)
    {
        p.d0_ = d0[i];
        p.position0_ = position0[i];
        p.sigma_ = sigma[i];
        p.mu_ = mu[i];
        p.liquidCore_ = liquidCore[i];
        p.KHindex_ = KHindex[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        p.tc_ = tc[i];
        p.ms_ = ms[i];
        p.injector_ = injector[i];
        p.tMom_ = tMom[i];
        p.user_ = user[i];

        ++i;
    }

    // added <<
    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    Pout << __FILE__ << ": " << __LINE__ << " => idGas: " << idGas << endl;
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();

    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    // Set storage for each Y... for each parcel
    for (SolidifyingSprayParcel<ParcelType>& p : c)
    {
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }

    // Pout << __FILE__ << ": " << __LINE__ << " => idGas: " << idGas << endl;
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        // Pout << __FILE__ << ": " << __LINE__ << $p.YGas_ << endl;
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (SolidifyingSprayParcel<ParcelType>& p : c)
        {
            p.YGas_[j] = YGas[i]/(max(p.Y()[GAS], SMALL));
            ++i;
        }
    }
    // Populate YLiquid for each parcel
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (SolidifyingSprayParcel<ParcelType>& p : c)
        {
            p.YLiquid_[j] = YLiquid[i]/(max(p.Y()[LIQ], SMALL));
            ++i;
        }
    }
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            ),
            valid
        );

        label i = 0;
        for (SolidifyingSprayParcel<ParcelType>& p : c)
        {
            p.YSolid_[j] = YSolid[i]/(max(p.Y()[SLD], SMALL));
            ++i;
        }
    }
    // added >>
}


template<class ParcelType>
template<class CloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SolidifyingSprayParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    const label np = c.size();
    const bool valid = np;


    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::NO_READ), np);
    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::NO_READ),
        np
    );
    IOField<scalar> sigma(c.fieldIOobject("sigma", IOobject::NO_READ), np);
    IOField<scalar> mu(c.fieldIOobject("mu", IOobject::NO_READ), np);
    IOField<scalar> liquidCore
    (
        c.fieldIOobject("liquidCore", IOobject::NO_READ),
        np
    );
    IOField<scalar> KHindex(c.fieldIOobject("KHindex", IOobject::NO_READ), np);
    IOField<scalar> y(c.fieldIOobject("y", IOobject::NO_READ), np);
    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::NO_READ), np);
    IOField<scalar> tc(c.fieldIOobject("tc", IOobject::NO_READ), np);
    IOField<scalar> ms(c.fieldIOobject("ms", IOobject::NO_READ), np);
    IOField<scalar> injector
    (
        c.fieldIOobject("injector", IOobject::NO_READ),
        np
    );
    IOField<scalar> tMom(c.fieldIOobject("tMom", IOobject::NO_READ), np);
    IOField<scalar> user(c.fieldIOobject("user", IOobject::NO_READ), np);

    label i = 0;
    for (const SolidifyingSprayParcel<ParcelType>& p : c)
    {
        d0[i] = p.d0_;
        position0[i] = p.position0_;
        sigma[i] = p.sigma_;
        mu[i] = p.mu_;
        liquidCore[i] = p.liquidCore_;
        KHindex[i] = p.KHindex_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        tc[i] = p.tc_;
        ms[i] = p.ms_;
        injector[i] = p.injector_;
        tMom[i] = p.tMom_;
        user[i] = p.user_;
        ++i;
    }

    d0.write(valid);
    position0.write(valid);
    sigma.write(valid);
    mu.write(valid);
    liquidCore.write(valid);
    KHindex.write(valid);
    y.write(valid);
    yDot.write(valid);
    tc.write(valid);
    ms.write(valid);
    injector.write(valid);
    tMom.write(valid);
    user.write(valid);

    // added <<
    Pout << __FILE__ << ": " << __LINE__ << " is reached : " << __FUNCTION__ << endl;
    // Write the composition fractions
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YGas[i] = p0.YGas()[j]*max(p0.Y()[GAS], SMALL);
                ++i;
            }

            YGas.write(np > 0);
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YLiquid[i] = p0.YLiquid()[j]*max(p0.Y()[LIQ], SMALL);
                ++i;
            }

            YLiquid.write(np > 0);
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YSolid[i] = p0.YSolid()[j]*max(p0.Y()[SLD], SMALL);
                ++i;
            }

            YSolid.write(np > 0);
        }
    }
    // added >>
}


template<class ParcelType>
void Foam::SolidifyingSprayParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("d0", d0_);
    writeProp("position0", position0_);
    writeProp("sigma", sigma_);
    writeProp("mu", mu_);
    writeProp("liquidCore", liquidCore_);
    writeProp("KHindex", KHindex_);
    writeProp("y", y_);
    writeProp("yDot", yDot_);
    writeProp("tc", tc_);
    writeProp("ms", ms_);
    writeProp("injector", injector_);
    writeProp("tMom", tMom_);
    writeProp("user", user_);
    writeProp("YGas", YGas_);   // added
    writeProp("YLiquid", YLiquid_);   // added
    writeProp("YSolid", YSolid_);   // added

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);
}


template<class ParcelType>
template<class CloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SolidifyingSprayParcel<ParcelType>::readObjects
(
    CloudType& c,
    const CompositionType& compModel,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, compModel, obr);

    if (!c.size()) return;

    const auto& d0 = cloud::lookupIOField<scalar>("d0", obr);
    const auto& position0 = cloud::lookupIOField<vector>("position0", obr);
    const auto& sigma = cloud::lookupIOField<scalar>("sigma", obr);
    const auto& mu = cloud::lookupIOField<scalar>("mu", obr);
    const auto& liquidCore = cloud::lookupIOField<scalar>("liquidCore", obr);
    const auto& KHindex = cloud::lookupIOField<scalar>("KHindex", obr);
    const auto& y = cloud::lookupIOField<scalar>("y", obr);
    const auto& yDot = cloud::lookupIOField<scalar>("yDot", obr);
    const auto& tc = cloud::lookupIOField<scalar>("tc", obr);
    const auto& ms = cloud::lookupIOField<scalar>("ms", obr);
    const auto& injector = cloud::lookupIOField<scalar>("injector", obr);
    const auto& tMom = cloud::lookupIOField<scalar>("tMom", obr);
    const auto& user = cloud::lookupIOField<scalar>("user", obr);

    label i = 0;
    for (SolidifyingSprayParcel<ParcelType>& p : c)
    {
        p.d0_ = d0[i];
        p.position0_ = position0[i];
        p.sigma_ = sigma[i];
        p.mu_ = mu[i];
        p.liquidCore_ = liquidCore[i];
        p.KHindex_ = KHindex[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        p.tc_ = tc[i];
        p.ms_ = ms[i];
        p.injector_ = injector[i];
        p.tMom_ = tMom[i];
        p.user_ = user[i];

        ++i;
    }

    // added <<
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    ParcelType::readObjects(c, obr);

    const label np = c.size();

    // The composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            const word fieldName = "Y" + gasNames[j] + stateLabels[idGas];
            const auto& YGas = cloud::lookupIOField<scalar>(fieldName, obr);

            label i = 0;
            for (SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                p0.YGas()[j]*max(p0.Y()[GAS], SMALL) = YGas[i];
                ++i;
            }
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            const word fieldName = "Y" + liquidNames[j] + stateLabels[idLiquid];
            const auto& YLiquid = cloud::lookupIOField<scalar>(fieldName, obr);

            label i = 0;
            for (SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                p0.YLiquid()[j]*max(p0.Y()[LIQ], SMALL) = YLiquid[i];
                ++i;
            }
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            const word fieldName = "Y" + solidNames[j] + stateLabels[idSolid];
            const auto& YSolid = cloud::lookupIOField<scalar>(fieldName, obr);

            label i = 0;
            for (SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                p0.YSolid()[j]*max(p0.Y()[SLD], SMALL) = YSolid[i];
                ++i;
            }
        }
    }
    // added >>
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SolidifyingSprayParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    const CompositionType& compModel,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, compModel, obr);

    const label np = c.size();

    auto& d0 = cloud::createIOField<scalar>("d0", np, obr);
    auto& position0 = cloud::createIOField<vector>("position0", np, obr);
    auto& sigma = cloud::createIOField<scalar>("sigma", np, obr);
    auto& mu = cloud::createIOField<scalar>("mu", np, obr);
    auto& liquidCore = cloud::createIOField<scalar>("liquidCore", np, obr);
    auto& KHindex = cloud::createIOField<scalar>("KHindex", np, obr);
    auto& y = cloud::createIOField<scalar>("y", np, obr);
    auto& yDot= cloud::createIOField<scalar>("yDot", np, obr);
    auto& tc = cloud::createIOField<scalar>("tc", np, obr);
    auto& ms = cloud::createIOField<scalar>("ms", np, obr);
    auto& injector = cloud::createIOField<scalar>("injector", np, obr);
    auto& tMom = cloud::createIOField<scalar>("tMom", np, obr);
    auto& user = cloud::createIOField<scalar>("user", np, obr);

    label i = 0;
    for (const SolidifyingSprayParcel<ParcelType>& p : c)
    {
        d0[i] = p.d0_;
        position0[i] = p.position0_;
        sigma[i] = p.sigma_;
        mu[i] = p.mu_;
        liquidCore[i] = p.liquidCore_;
        KHindex[i] = p.KHindex_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        tc[i] = p.tc_;
        ms[i] = p.ms_;
        injector[i] = p.injector_;
        tMom[i] = p.tMom_;
        user[i] = p.user_;

        ++i;
    }

    // added <<
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    ParcelType::writeObjects(c, obr);

    // const label np = c.size();   // Redeclaration

    // Write the composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            const word fieldName = "Y" + gasNames[j] + stateLabels[idGas];
            auto& YGas = cloud::createIOField<scalar>(fieldName, np, obr);

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YGas[i] = p0.YGas()[j]*max(p0.Y()[GAS], SMALL);
                ++i;
            }
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            const word fieldName = "Y" + liquidNames[j] + stateLabels[idLiquid];
            auto& YLiquid = cloud::createIOField<scalar>(fieldName, np, obr);

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YLiquid[i] = p0.YLiquid()[j]*max(p0.Y()[LIQ], SMALL);
                ++i;
            }
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            const word fieldName = "Y" + solidNames[j] + stateLabels[idSolid];
            auto& YSolid = cloud::createIOField<scalar>(fieldName, np, obr);

            label i = 0;
            for (const SolidifyingSprayParcel<ParcelType>& p0 : c)
            {
                YSolid[i] = p0.YSolid()[j]*max(p0.Y()[SLD], SMALL);
                ++i;
            }
        }
    }
    // added >>
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SolidifyingSprayParcel<ParcelType>& p
)
{
    scalarField YGasLoc(p.YGas());  // added
    Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    scalarField YLiquidLoc(p.YLiquid());  // added
    scalarField YSolidLoc(p.YSolid());  // added

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.d0()
            << token::SPACE << p.position0()
            << token::SPACE << p.sigma()
            << token::SPACE << p.mu()
            << token::SPACE << p.liquidCore()
            << token::SPACE << p.KHindex()
            << token::SPACE << p.y()
            << token::SPACE << p.yDot()
            << token::SPACE << p.tc()
            << token::SPACE << p.ms()
            << token::SPACE << p.injector()
            << token::SPACE << p.tMom()
            << token::SPACE << p.user()
            << token::SPACE << YGasLoc  // added
            << token::SPACE << YLiquidLoc   // added
            << token::SPACE << YSolidLoc;   // added

    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os  << YGasLoc << YLiquidLoc << YSolidLoc;  // added 
        os.write
        (
            reinterpret_cast<const char*>(&p.d0_),
            SolidifyingSprayParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
