/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SolidifyingSprayCloud.H"
#include "AtomizationModel.H"
#include "BreakupModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::setModels()
{
    atomizationModel_.reset
    (
        AtomizationModel<SolidifyingSprayCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    breakupModel_.reset
    (
        BreakupModel<SolidifyingSprayCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::cloudReset
(
    SolidifyingSprayCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    atomizationModel_.reset(c.atomizationModel_.ptr());
    breakupModel_.reset(c.breakupModel_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SolidifyingSprayCloud<CloudType>::SolidifyingSprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    solidifyingSprayCloud(),
    cloudCopyPtr_(nullptr),
    // constProps_(this->particleProperties()), // added
    averageParcelMass_(0.0),
    atomizationModel_(nullptr),
    breakupModel_(nullptr)
{
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    if (this->solution().active())
    {
        setModels();
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;

        averageParcelMass_ = this->injectors().averageParcelMass();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }

        Info << "Average parcel mass: " << averageParcelMass_ << endl;
    }

    if (this->solution().resetSourcesOnStartup())
    {
        CloudType::resetSourceTerms();
    }
}


template<class CloudType>
Foam::SolidifyingSprayCloud<CloudType>::SolidifyingSprayCloud
(
    SolidifyingSprayCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    solidifyingSprayCloud(),
    cloudCopyPtr_(nullptr),
    // constProps_(c.constProps_), // added
    averageParcelMass_(c.averageParcelMass_),
    atomizationModel_(c.atomizationModel_->clone()),
    breakupModel_(c.breakupModel_->clone())
{}


template<class CloudType>
Foam::SolidifyingSprayCloud<CloudType>::SolidifyingSprayCloud
(
    const fvMesh& mesh,
    const word& name,
    const SolidifyingSprayCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    solidifyingSprayCloud(),
    cloudCopyPtr_(nullptr),
    // constProps_(),  // added
    averageParcelMass_(0.0),
    atomizationModel_(nullptr),
    breakupModel_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SolidifyingSprayCloud<CloudType>::~SolidifyingSprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    const liquidMixtureProperties& liqMix = this->composition().liquids();
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    label idGas = this->composition().idGas();  // added
    label idLiquid = this->composition().idLiquid();  // added
    label idSolid = this->composition().idSolid();  // added

    const scalarField& Y(parcel.Y());
    scalarField X(liqMix.X(Y));
    const scalar pc = this->p()[parcel.cell()];

    // override rho and Cp from constantProperties
    parcel.Cp() = liqMix.Cp(pc, parcel.T(), X);
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    parcel.rho() = liqMix.rho(pc, parcel.T(), X);
    parcel.sigma() = liqMix.sigma(pc, parcel.T(), X);
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    parcel.mu() = liqMix.mu(pc, parcel.T(), X);
    parcel.YGas() = this->composition().Y0(idGas);  // added
    parcel.YLiquid() = this->composition().Y0(idLiquid);  // added
    parcel.YSolid() = this->composition().Y0(idSolid);  // added

    // If rho0 was given in constProp use it. If not use the composition
    // to set tho
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    if (constProps_.rho0() == -1)   // added
    {
        const scalarField& Ygas = this->composition().Y0(idGas);
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        const scalarField& Yliq = this->composition().Y0(idLiquid);
        const scalarField& Ysol = this->composition().Y0(idSolid);
        const scalar p0 =
            this->composition().thermo().thermo().p()[parcel.cell()];
        const scalar T0 = constProps_.T0();

        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        parcel.rho() = this->composition().rho(Ygas, Yliq, Ysol, T0, p0);
    }

}


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();

    parcel.liquidCore() = atomization().initLiquidCore();

    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;

    if (fullyDescribed) // added
    {
        label idGas = this->composition().idGas();
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        
        Pout << "SSC line 238" << endl,
        
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }

}


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<SolidifyingSprayCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>   // added
void Foam::SolidifyingSprayCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::evolve()
{
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        this->solve(*this, td);
    }
    else
    {
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    }
}


/*template<class CloudType>   // added
void Foam::SolidifyingSprayCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}*/


template<class CloudType>
void Foam::SolidifyingSprayCloud<CloudType>::info()
{
    CloudType::info();
    scalar d32 = 1.0e+6*this->Dij(3, 2);
    scalar d10 = 1.0e+6*this->Dij(1, 0);
    scalar dMax = 1.0e+6*this->Dmax();
    scalar pen = this->penetration(0.95);

    Info << "    D10, D32, Dmax (mu)             = " << d10 << ", " << d32
         << ", " << dMax << nl
         << "    Liquid penetration 95% mass (m) = " << pen << endl;
}

template<class CloudType>   // added
void Foam::SolidifyingSprayCloud<CloudType>::writeFields() const
{
    Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
    if (this->compositionModel_)
    {
        Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
