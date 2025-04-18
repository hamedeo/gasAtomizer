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

Class
    Foam::SolidifyingSprayParcel

Description
    Reacting spray parcel, with added functionality for atomization and breakup

\*---------------------------------------------------------------------------*/

#ifndef SolidifyingSprayParcel_H
#define SolidifyingSprayParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class SolidifyingSprayParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const SolidifyingSprayParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                          Class SolidifyingSprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class SolidifyingSprayParcel
:
    public ParcelType
{
public:

    //- Size in bytes of the fields
    // static const std::size_t sizeofFields;  // added


    // IDs of phases in ReactingParcel phase list (Y)

        static const label GAS;  // added
        static const label LIQ;  // added
        static const label SLD;  // added


    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Particle initial surface tension [N/m]
            demandDrivenEntry<scalar> sigma0_;

            //- Particle initial dynamic viscosity [Pa.s]
            demandDrivenEntry<scalar> mu0_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);

            //- Construct from components
            constantProperties
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
            );


        // Access

            //- Return const access to the initial surface tension
            inline scalar sigma0() const;

            //- Return const access to the initial dynamic viscosity
            inline scalar mu0() const;

    };


    //- Use base tracking data
    typedef typename ParcelType::trackingData trackingData;

// added <<
private:

    // Private Member Functions

        //- Return the solidification switch
        //template<class TrackCloudType>
        bool Solidification
        (
            // TrackCloudType& cloud,
            // trackingData& td,
            // const scalar p,
            const scalar T,
            const scalar Tth
            // scalarField& YLiquid_,
            // scalarField& YSolid_
        );

        //- Return the mixture effective specific heat capacity
        template<class TrackCloudType>
        scalar CpEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackCloudType>
        scalar HsEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackCloudType>
        scalar LEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractions
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );
// added >>

protected:

    // Protected data

        // Spray parcel properties

            //- Initial droplet diameter
            scalar d0_;

            //- Injection position
            vector position0_;

            //- Liquid surface tension [N/m]
            scalar sigma_;

            //- Liquid dynamic viscosity [Pa.s]
            scalar mu_;

            //- Part of liquid core ( >0.5=liquid, <0.5=droplet )
            scalar liquidCore_;

            //- Index for KH Breakup
            scalar KHindex_;

            //- Spherical deviation
            scalar y_;

            //- Rate of change of spherical deviation
            scalar yDot_;

            //- Characteristic time (used in atomization and/or breakup model)
            scalar tc_;

            //- Stripped parcel mass due to breakup
            scalar ms_;

            //- Injected from injector (needed e.g. for calculating distance
            //  from injector)
            scalar injector_;

            //- Momentum relaxation time (needed for calculating parcel acc.)
            scalar tMom_;

            //- Passive scalar (extra variable to be defined by user)
            scalar user_;

            //- Mass fractions of gases []
            scalarField YGas_;  // added

            //- Mass fractions of liquids []
            scalarField YLiquid_;  // added

            //- Mass fractions of solids []
            scalarField YSolid_;  // added

            //- Mass fractions Mixture []
            scalarField YMixture0_; //added

            //- Mass fractions of mixture []
            // scalarField Y_;

    // added <<
    // Protected Member Functions


        //- Return change of volume due to mass exchange
        template<class TrackCloudType>
        scalar updatedDeltaVolume
        (
            TrackCloudType& cloud,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid,
            const label idG,
            const label idL,
            const label idS,
            const scalar p,
            const scalar T
        );
    // added >>


public:

    // Static data members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields;

        //- Runtime type information
        TypeName("SolidifyingSprayParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " d0"
          + " position0"
          + " sigma"
          + " mu"
          + " liquidCore"
          + " KHindex"
          + " y"
          + " yDot"
          + " tc"
          + " ms"
          + " injector"
          + " tMom"
          + " user"
          + " nGas(Y1..YN)" // added
          + " nLiquid(Y1..YN)" // added
          + " nSolid(Y1..YN)" // added
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline SolidifyingSprayParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline SolidifyingSprayParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from components
        inline SolidifyingSprayParcel
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
            const scalarField& YGas0,   // added
            const scalarField& YLiquid0,   // added
            const scalarField& YSolid0,   // added
            // const constantProperties& constProps   // added
            const typename ParcelType::constantProperties& constProps
        );

        //- Construct from Istream
        SolidifyingSprayParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        SolidifyingSprayParcel
        (
            const SolidifyingSprayParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        SolidifyingSprayParcel(const SolidifyingSprayParcel& p);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new SolidifyingSprayParcel<ParcelType>(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>
            (
                new SolidifyingSprayParcel<ParcelType>(*this, mesh)
            );
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

            autoPtr<SolidifyingSprayParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<SolidifyingSprayParcel<ParcelType>>
                (
                    new SolidifyingSprayParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial droplet diameter
            inline scalar d0() const;

            //- Return const access to initial droplet position
            inline const vector& position0() const;

            //- Return const access to the liquid surface tension
            inline scalar sigma() const;

            //- Return const access to the liquid dynamic viscosity
            inline scalar mu() const;

            //- Return const access to liquid core
            inline scalar liquidCore() const;

            //- Return const access to Kelvin-Helmholtz breakup index
            inline scalar KHindex() const;

            //- Return const access to spherical deviation
            inline scalar y() const;

            //- Return const access to rate of change of spherical deviation
            inline scalar yDot() const;

            //- Return const access to atomization characteristic time
            inline scalar tc() const;

            //- Return const access to stripped parcel mass
            inline scalar ms() const;

            //- Return const access to injector id
            inline scalar injector() const;

            //- Return const access to momentum relaxation time
            inline scalar tMom() const;

            //- Return const access to passive user scalar
            inline scalar user() const;

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const; // added

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const; // added

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const; // added

            //- Return const access to the canCombust flag
            inline label canCombust() const; // added


        // Edit

            //- Return access to initial droplet diameter
            inline scalar& d0();

            //- Return access to initial droplet position
            inline vector& position0();

            //- Return access to the liquid surface tension
            inline scalar& sigma();

            //- Return access to the liquid dynamic viscosity
            inline scalar& mu();

            //- Return access to liquid core
            inline scalar& liquidCore();

            //- Return access to Kelvin-Helmholtz breakup index
            inline scalar& KHindex();

            //- Return access to spherical deviation
            inline scalar& y();

            //- Return access to rate of change of spherical deviation
            inline scalar& yDot();

            //- Return access to atomization characteristic time
            inline scalar& tc();

            //- Return access to stripped parcel mass
            inline scalar& ms();

            //- Return access to injector id
            inline scalar& injector();

            //- Return access to momentum relaxation time
            inline scalar& tMom();

            //- Return access to passive user scalar
            inline scalar& user();

            //- Return access to mass fractions of gases
            inline scalarField& YGas(); // added

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid(); // added

            //- Return access to mass fractions of solids
            inline scalarField& YSolid(); // added

            //- Return access to mass fractions of mixture []
            // inline scalarField& Y();

            //- Return access to the canCombust flag
            inline label& canCombust(); // added

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct parcel properties according to atomization model
            template<class TrackCloudType>
            void calcAtomization
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct parcel properties according to breakup model
            template<class TrackCloudType>
            void calcBreakup
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );
            /*
            //- Correct surface values due to emitted species
            template<class TrackCloudType>
            void correctSurfaceValues
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappa
            );*/

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Calculate the chi-factor for flash-boiling for the
            //  atomization model
            template<class TrackCloudType>
            scalar chi
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalarField& X
            ) const;

            //- Solve the TAB equation
            template<class TrackCloudType>
            void solveTABEq
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

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly
            ) const;

            //- Read particle fields as objects from the obr registry
            //  - no composition
            template<class CloudType>
            static void readObjects
            (
                CloudType& c,
                const objectRegistry& obr
            );

            //- Read particle fields as objects from the obr registry
            template<class CloudType, class CompositionType>
            static void readObjects
            (
                CloudType& c,
                const CompositionType& compModel,
                const objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            //  - no composition
            template<class CloudType>
            static void writeObjects
            (
                const CloudType& c,
                objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            template<class CloudType, class CompositionType>
            static void writeObjects
            (
                const CloudType& c,
                const CompositionType& compModel,
                objectRegistry& obr
            );


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const SolidifyingSprayParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SolidifyingSprayParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidifyingSprayParcel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
