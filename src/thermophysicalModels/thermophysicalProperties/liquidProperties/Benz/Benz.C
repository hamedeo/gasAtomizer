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

#include "Benz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Benz, 0);
    addToRunTimeSelectionTable(liquidProperties, Benz,);
    addToRunTimeSelectionTable(liquidProperties, Benz, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Benz::Benz()
:
    liquidProperties
    (
        18.015,
        647.13,
        2.2055e+7,
        0.05595,
        0.229,
        273.16,
        6.113e+2,
        373.15,
        6.1709e-30,
        0.3449,
        4.7813e+4
    ),
rho_(									//Density (Benz)
    842,
    ),
    pv_(73.649, -7258.2, -7.3037, 4.1653e-06, 2),
    hl_(647.13, 2889425.47876769, 0.3199, -0.212, 0.25795, 0),
    Cp_
    (
        15341.1046350264,
       -116.019983347211,
        0.451013044684985,
       -0.000783569247849015,
        5.20127671384957e-07,
        0
    ),
    h_
    (
       -17957283.7993676,
        15341.1046350264,
       -58.0099916736053,
        0.150337681561662,
       -0.000195892311962254,
        1.04025534276991e-07
    ),
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
    mu_(2.17e-3),							//Dynamic viscosity (Benz)
    mug_(2.6986e-06, 0.498, 1257.7, -19570),
    kappa_(-0.4267, 0.0056903, -8.0065e-06, 1.815e-09, 0, 0),
    kappag_(6.977e-05, 1.1243, 844.9, -148850),
    sigma_(0.02),							//Surface tension (Benz)
    D_(15.0, 15.0, 18.015, 28)
{}


Foam::Benz::Benz
(
    const liquidProperties& l,
    const thermophysicalConstant& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const thermophysicalConstant& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const thermophysicalConstant& surfaceTension
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

/*
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
        1811		//Tm -> normal melting point			(Fe)
    ),
    rho_				//density			(Fe) Beutl1994
    (
    8497.70,
    -0.82268,
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
    mu_(1e-2),				//dynamic viscosity		(Fe)
    kappa_(38),			//thermal conductivity		(Fe)
    sigma_(1.78)			//surface tension		(Fe)
{}


Foam::Fe::Fe
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
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    mu_(dynamicViscosity),
    kappa_(thermalConductivity),
    sigma_(surfaceTension)
{}
*/


Foam::Benz::Benz(const dictionary& dict)
:
    Benz()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Benz::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const Benz& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
