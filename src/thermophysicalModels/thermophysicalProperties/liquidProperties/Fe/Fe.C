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
        55.845,	//W -> molar mass				(Fe)
        9250,		//Tc -> critical temperature			(Fe) Beutl1994
        8750e5,	//Pc -> critical pressure			(Fe) Beutl1994
        1.758e-5,	//Vc -> critical volume			(Fe) Vc = Zc*R*Tc/pc
        0.2,		//Zc -> critical compressibility factor	(Fe) see: Kulinskii, V.L. (2014)
        1548.6,	//Tt -> triple point temperature		(Fe)
        0.16,		//Pt -> triple point pressure			(Fe)
        3135,		//Tb -> normal boiling point			(Fe)
        // 1811		//Tm or Tref -> normal melting point			(Fe)
        6.1709e-30,
        0.3449,
        4.7813e+4
    ),
    rho_				//density	(Fe) Brillo 2005 "Surface tension of ...  iron and their binary alloys" NSRDS0
    (
        8711.946,
        -0.926,
        0,
        0,
        0,
        0
    ),
    pv_				//vapour pressure		(Fe) Beutl1994
    (
        26.4317,
        -48769,
        -1.3217,
        0,
        0
    ),
    hl_(6340765),			//heat of vapourisation	(Fe)
    Cp_(825),				//heat capacity		(Fe)
    h_(2220432),			//enthalpy			(Fe)
    Cpg_
    (
        1851.73466555648,
        1487.53816264224,
        2609.3,
        493.366638912018,
        1167.6
    ),
    B_
    (
       -0.0012789342214821,
        1.4909797391063,
       -1563696.91923397,
        1.85445462114904e+19,
       -7.68082153760755e+21
    ),
    mu_     				//dynamic viscosity		(Fe) Assael 2006 "Reference Data for the Density ... Iron" NSRDS3   
    (
        0,
        0.00190151606,
        -6205.35,
        1
    ),
    mug_(2.6986e-06, 0.498, 1257.7, -19570),
    kappa_(38),	    		//thermal conductivity		(Fe)
    kappag_(6.977e-05, 1.1243, 844.9, -148850),
    sigma_      		//surface tension		(Fe) Reference Data for the Density and Viscosity of Liquid Aluminum and Liquid Iron" NSRDS0
    (
        1.900586,         
        -3.97e-4,
        0,
        0,
        0,
        0
    ),
    // sigma_(1.2)
    D_(15.0, 15.0, 18.015, 28)
{}


Foam::Fe::Fe
(
    const liquidProperties& l,
    const NSRDSfunc0& density,
    const NSRDSfunc1& vapourPressure,
    const thermophysicalConstant& heatOfVapourisation,
    const thermophysicalConstant& heatCapacity,
    const thermophysicalConstant& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc3& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const thermophysicalConstant& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc0& surfaceTension,
    // const thermophysicalConstant& surfaceTension
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
