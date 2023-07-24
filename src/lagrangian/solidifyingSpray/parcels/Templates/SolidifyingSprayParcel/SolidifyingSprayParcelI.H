/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
inline Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    sigma0_(this->dict_, 0.0),
    mu0_(this->dict_, 0.0)
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    sigma0_(cp.sigma0_),
    mu0_(cp.mu0_)
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    sigma0_(this->dict_, "sigma0"),
    mu0_(this->dict_, "mu0")
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::constantProperties
(
    const label parcelTypeId,
    const scalar rhoMin,
    const scalar rho0,
    const scalar minParcelMass,
    const scalar youngsModulus,
    const scalar poissonsRatio,
    const scalar T0,
    const scalar TMin,
    const scalar TMax,
    const scalar Cp0,
    const scalar epsilon0,
    const scalar f0,
    const scalar Pr,
    const scalar pMin,
    const bool constantVolume,
    const scalar sigma0,
    const scalar mu0
)
:
    ParcelType::constantProperties
    (
        parcelTypeId,
        rhoMin,
        rho0,
        minParcelMass,
        youngsModulus,
        poissonsRatio,
        T0,
        TMin,
        TMax,
        Cp0,
        epsilon0,
        f0,
        Pr,
        pMin,
        constantVolume
    ),
    sigma0_(this->dict_, sigma0),
    mu0_(this->dict_, mu0)
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    d0_(this->d()),
    position0_(this->position()),
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
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    d0_(this->d()),
    position0_(this->position()),
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
{}


template<class ParcelType>
inline Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel
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
    const scalar liquidCore,
    const scalar KHindex,
    const scalar y,
    const scalar yDot,
    const scalar tc,
    const scalar ms,
    const scalar injector,
    const scalar tMom,
    const scalar user,
    const scalarField& YGas0,
    const scalarField& YLiquid0,    // added
    const scalarField& YSolid0,    // added
    // const constantProperties& constProps    // added
    const typename ParcelType::constantProperties& constProps
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
    ),
    d0_(d0),
    position0_(this->position()),
    sigma_(constProps.sigma0()),
    mu_(constProps.mu0()),
    liquidCore_(liquidCore),
    KHindex_(KHindex),
    y_(y),
    yDot_(yDot),
    tc_(tc),
    ms_(ms),
    injector_(injector),
    tMom_(tMom),
    user_(user),
    YGas_(YGas0),   // added
    YLiquid_(YLiquid0),   // added
    YSolid_(YSolid0)   // added
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::sigma0() const
{
    return sigma0_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::SolidifyingSprayParcel<ParcelType>::constantProperties::mu0() const
{
    return mu0_.value();
}


// * * * * * * * * * * SolidifyingSprayParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>  // added
inline const Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>  // added
inline const Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>  // added
inline const Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>  // added
inline Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>  // added
inline Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>  // added
inline Foam::scalarField& Foam::SolidifyingSprayParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::d0() const
{
    return d0_;
}


template<class ParcelType>
inline const Foam::vector& Foam::SolidifyingSprayParcel<ParcelType>::position0() const
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::sigma() const
{
    return sigma_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::mu() const
{
    return mu_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::liquidCore() const
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::KHindex() const
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::y() const
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::yDot() const
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::tc() const
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::ms() const
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::injector() const
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::tMom() const
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::user() const
{
    return user_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::d0()
{
    return d0_;
}


template<class ParcelType>
inline Foam::vector& Foam::SolidifyingSprayParcel<ParcelType>::position0()
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::sigma()
{
    return sigma_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::mu()
{
    return mu_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::liquidCore()
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::KHindex()
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::y()
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::yDot()
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::tc()
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::ms()
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::injector()
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::tMom()
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SolidifyingSprayParcel<ParcelType>::user()
{
    return user_;
}


// ************************************************************************* //
