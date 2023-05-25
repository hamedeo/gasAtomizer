/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "myDistortedSphereDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::myDistortedSphereDragForce<CloudType>::CdRe
(
    const scalar Re
) const
{
    // (AOB:Eq. 35; LMR:Eq. 9)
    if (Re > 1000.0)
    {
            Info << "Re" << Re << endl;
        return 0.424*Re;
    }

    Info << "Re" << Re << endl;

    return 24.0*(1.0 + (1.0/6.0)*pow(Re, 2.0/3.0));

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::myDistortedSphereDragForce<CloudType>::myDistortedSphereDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{
    Info << " Test 1" << endl;
}


template<class CloudType>
Foam::myDistortedSphereDragForce<CloudType>::myDistortedSphereDragForce
(
    const myDistortedSphereDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{
    Info << " Test 2" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::myDistortedSphereDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    // Limit the drop distortion to y=0 (sphere) and y=1 (disk)
    const scalar y = min(max(p.y(), scalar(0)), scalar(1));

    // I have added this lines that can be safely removed #deoRemovable    
    // Log_<< "    Parcel fate: system (number, mass)" << nl   // #deoRemovable
     Info<< "Level of distortion y: " << y << endl;     // #deoRemovable

    //Here is what i added to test but you can remove it
    //Pout<< "y = " << y << endl;
 
    // (LMR:Eq. 10)
    return
        forceSuSp
        (
            Zero,
            mass*0.75*muc*CdRe(Re)*(1.0 + 2.632*y)/(p.rho()*sqr(p.d()))
        );
}


// ************************************************************************* //
