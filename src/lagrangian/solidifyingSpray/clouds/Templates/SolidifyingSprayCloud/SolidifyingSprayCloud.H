/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::SolidifyingSprayCloud

Description
    Templated base class for spray cloud

    - sub-models:
      - atomization model
      - break-up model

\*---------------------------------------------------------------------------*/

#ifndef SolidifyingSprayCloud_H
#define SolidifyingSprayCloud_H

#include "solidifyingSprayCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class CloudType> class AtomizationModel;
template<class CloudType> class BreakupModel;
template<class CloudType> class SolidifyingSprayCloud;


/*---------------------------------------------------------------------------*\
                      Class SolidifyingSprayCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SolidifyingSprayCloud
:
    public CloudType,
    public solidifyingSprayCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef SolidifyingSprayCloud<CloudType> solidifyingSprayCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<SolidifyingSprayCloud<CloudType>> cloudCopyPtr_;

        //- Average parcel mass
        scalar averageParcelMass_;


    // Private Member Functions

        //- No copy construct
        SolidifyingSprayCloud(const SolidifyingSprayCloud&) = delete;

        //- No copy assignment
        void operator=(const SolidifyingSprayCloud&) = delete;


protected:

    // Protected data

        // References to the cloud sub-models

            //- Atomization model
            autoPtr<AtomizationModel<SolidifyingSprayCloud<CloudType>>>
                atomizationModel_;

            //- Break-up model
            autoPtr<BreakupModel<SolidifyingSprayCloud<CloudType>>> breakupModel_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(SolidifyingSprayCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        SolidifyingSprayCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        SolidifyingSprayCloud(SolidifyingSprayCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        SolidifyingSprayCloud
        (
            const fvMesh& mesh,
            const word& name,
            const SolidifyingSprayCloud<CloudType>& c
        );


        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SolidifyingSprayCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new SolidifyingSprayCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~SolidifyingSprayCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const SolidifyingSprayCloud& cloudCopy() const;

            //- Return const-access to the average parcel mass
            inline scalar averageParcelMass() const;


        // Check

            //- Penetration for fraction [0-1] of the current total mass
            inline scalar penetration(const scalar fraction) const;


            // Sub-models

                //- Return const-access to the atomization model
                inline const AtomizationModel<SolidifyingSprayCloud<CloudType>>&
                    atomization() const;

                //- Return reference to the atomization model
                inline AtomizationModel<SolidifyingSprayCloud<CloudType>>& atomization();

                //- Return const-access to the breakup model
                inline const BreakupModel<SolidifyingSprayCloud<CloudType>>&
                    breakup() const;

                //- Return reference to the breakup model
                inline BreakupModel<SolidifyingSprayCloud<CloudType>>& breakup();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Evolve the spray (inject, move)
            void evolve();


        // I-O

            //- Print cloud information
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SolidifyingSprayCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SolidifyingSprayCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
