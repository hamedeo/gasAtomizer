/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::basicSolidifyingSprayParcel

Description
    Definition of spray parcel

SourceFiles
    basicSolidifyingSprayParcel.C
    basicSolidifyingSprayParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef basicSolidifyingSprayParcel_H
#define basicSolidifyingSprayParcel_H

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "ThermoParcel.H"
#include "ReactingParcel.H"
#include "SolidifyingSprayParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef SolidifyingSprayParcel
    <
        ReactingParcel
        <
            ThermoParcel
            <
                KinematicParcel
                <
                    particle
                >
            >
        >
    > basicSolidifyingSprayParcel;

    //- Non-contiguous data for basicSolidifyingSprayParcel
    template<> struct is_contiguous<basicSolidifyingSprayParcel> : std::false_type {};
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
