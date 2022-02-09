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

#include "Fe.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Fe, 0);
    addToRunTimeSelectionTable(liquidProperties, Fe,);
    addToRunTimeSelectionTable(liquidProperties, Fe, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Fe::Fe()
:
    liquidProperties
    (
        55.845,	//Fe
        8500,		//Fe
        510e8,		//Au -> Fe not available
        2.6823e-7,	// Vc = (Zc*R*Tc)/pc	
        0.2,		// Near 0.2 for liquid metals, see: Kulinskii, V.L. (2014), The critical compressibility factor value: associative fluids and liquid alkali metals
        1548.6,	//Fe
        0.16,		//Fe
        3135,		//Fe
        6.1709e-30,	//H2O
        0.3449,	//H2O
        4.7813e+4	//H2O
    ),
    rho_(6980, 0, 0, 0, 0, 0),		//Fe
    pv_(11.353, -19574, 0, 0, 0),		//Fe
    hl_(6340765, 0, 0, 0, 0, 0),		//Fe
    Cp_(450, 0, 0, 0, 0, 0),			//Fe
    h_(2220432, 0, 0, 0, 0, 0),		//Fe
    Cpg_(477.1, 0, 0, 0, 0, 0),		//Air
    B_						//H2O
    (
       -0.0012789342214821,
        1.4909797391063,
       -1563696.91923397,
        1.85445462114904e+19,
       -7.68082153760755e+21
    ),
    mu_(1e-2, 0, 0, 0, 0, 0),			//Fe
    mug_(2e-5, 0, 0, 0, 0, 0),		//Guess, based on air
    kappa_(38, 0, 0, 0, 0, 0),		//Fe
    kappag_(6.977e-05, 1.1243, 844.9, -148850),	//H2O
    sigma_(1.78, 0, 0, 0, 0, 0),		//Fe
    D_(15.0, 15.0, 18.015, 28)		//H2O
{}


Foam::Fe::Fe
(
            /*const liquidProperties& l,
            const NSRDSfunc5& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const NSRDSfunc0& heatCapacity,
            const NSRDSfunc0& enthalpy,
            const NSRDSfunc7& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc1& dynamicViscosity,
            const NSRDSfunc2& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const NSRDSfunc2& vapourThermalConductivity,
            const NSRDSfunc6& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity*/
            
            const liquidProperties& l,
            const NSRDSfunc0& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc0& heatOfVapourisation,
            const NSRDSfunc0& heatCapacity,
            const NSRDSfunc0& enthalpy,
            const NSRDSfunc0& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc0& dynamicViscosity,
            const NSRDSfunc0& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const NSRDSfunc2& vapourThermalConductivity,
            const NSRDSfunc0& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::Fe::Fe(const dictionary& dict)
:
    Fe()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Fe::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const Fe& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
