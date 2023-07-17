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
    Foam::KelvinHelmholtz

Group
    grpLagrangianIntermediateBreakupSubModels

Description
    Secondary breakup model which uses the Kelvin-Helmholtz
    instability theory to predict the 'stripped' droplets... and
    the Raleigh-Taylor instability as well.

\*---------------------------------------------------------------------------*/

#include "myBreakupModel.H"

#ifndef KelvinHelmholtz_H
#define KelvinHelmholtz_H


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class KelvinHelmholtz Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KelvinHelmholtz
:
    public myBreakupModel<CloudType>
{
private:

    // Private data

        // model constants
        scalar b0_;
        scalar b1_;
        scalar cTau_;
        scalar cRT_;
        scalar msLimit_;
        scalar weberLimit_;


public:

    //- Runtime type information
    TypeName("KelvinHelmholtz");


    // Constructors

        //- Construct from dictionary
        KelvinHelmholtz(const dictionary&, CloudType&);

        //- Construct copy
        KelvinHelmholtz(const KelvinHelmholtz<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<myBreakupModel<CloudType>> clone() const
        {
            return autoPtr<myBreakupModel<CloudType>>
            (
                new KelvinHelmholtz<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~KelvinHelmholtz();


    // Member Functions

        //- Update the parcel diameter
        virtual bool update
        (
            const scalar dt,
            const vector& g,
            scalar& d,
            scalar& tc,
            scalar& ms,
            scalar& nParticle,
            scalar& KHindex,
            scalar& y,
            scalar& yDot,
            const scalar d0,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const vector& U,
            const scalar rhoc,
            const scalar muc,
            const vector& Urel,
            const scalar Urmag,
            const scalar tMom,
            scalar& dChild,
            scalar& massChild
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KelvinHelmholtz.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
