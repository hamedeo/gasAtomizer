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

Class
    Foam::Fe

Description
    water

SourceFiles
    Fe.C

\*---------------------------------------------------------------------------*/

#ifndef Fe_H
#define Fe_H

#include "../liquidProperties/liquidProperties.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc5.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "APIdiffCoefFunc.H"
#include "thermophysicalConstant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Fe Declaration
\*---------------------------------------------------------------------------*/

class Fe
:
    public liquidProperties
{
    // Private data

        NSRDSfunc0 rho_;
        NSRDSfunc1 pv_;
        thermophysicalConstant hl_;
        thermophysicalConstant Cp_;
        thermophysicalConstant h_;
        thermophysicalConstant mu_;
        thermophysicalConstant kappa_;
        thermophysicalConstant sigma_;


public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("Fe");


    // Constructors

        //- Construct null
        Fe();

        //- Construct from components
        Fe
        (
            const liquidProperties& l,
            const NSRDSfunc0& density,
            const NSRDSfunc1& vapourPressure,
            const thermophysicalConstant& heatOfVapourisation,
            const thermophysicalConstant& heatCapacity,
            const thermophysicalConstant& enthalpy,
            const thermophysicalConstant& dynamicViscosity,
            const thermophysicalConstant& thermalConductivity,
            const thermophysicalConstant& surfaceTension
        );

        //- Construct from dictionary
        Fe(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>::NewFrom<Fe>(*this);
        }


    // Member Functions

        //- Liquid density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //- Vapour pressure [Pa]
        inline scalar pv(scalar p, scalar T) const;

        //- Heat of vapourisation [J/kg]
        inline scalar hl(scalar p, scalar T) const;

        //- Liquid heat capacity [J/(kg K)]
        inline scalar Cp(scalar p, scalar T) const;

        //- Liquid Enthalpy [J/(kg)]
        inline scalar h(scalar p, scalar T) const;

        //- Liquid viscosity [Pa s]
        inline scalar mu(scalar p, scalar T) const;

        //- Liquid thermal conductivity  [W/(m K)]
        inline scalar kappa(scalar p, scalar T) const;

        //- Surface tension [N/m]
        inline scalar sigma(scalar p, scalar T) const;

        //- Vapour diffusivity [m2/s]
        inline scalar D(scalar p, scalar T) const;

        //- Vapour diffusivity [m2/s] with specified binary pair
        inline scalar D(scalar p, scalar T, scalar Wb) const;


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const Fe& l);
};


Ostream& operator<<(Ostream& os, const Fe& l);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
