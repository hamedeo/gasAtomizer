/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
//#include "myBasicSprayCloud.H"
#include "basicSprayCloud.H"

/*#include "makeReactingParcelCloudFunctionObjects.H"

// Kinematic

#include "makeThermoParcelForces.H" // thermo variant
#include "makeThermoParcelTurbulenceForces.H" // add turbulence variant
#include "makeParcelDispersionModels.H"
#include "makeParcelTurbulenceDispersionModels.H" // add turbulence variant
#include "makeSprayParcelInjectionModels.H" // Spray variant
#include "makeParcelPatchInteractionModels.H"
#include "makeSprayParcelStochasticCollisionModels.H" // Spray variant


// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting

#include "makeMyReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"


// Spray
*/
//#include "myDistortedSphereDragForce.H"
/*#include "DistortedSphereDragForce.H"
//#include "makeSprayParcelAtomizationModels.H"
//#include "makeSprayParcelBreakupModels.H"
*/
#include "myMakeSprayParcelBreakupModels.H"
//#include "myMakeSprayParcelDragModels.H"
/*
// MPPIC sub-models

#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Kinematic sub-models

makeThermoParcelForces(myBasicSprayCloud);
makeThermoParcelTurbulenceForces(myBasicSprayCloud);
makeParcelDispersionModels(myBasicSprayCloud);
makeParcelTurbulenceDispersionModels(myBasicSprayCloud);
makeSprayParcelInjectionModels(myBasicSprayCloud);
makeParcelPatchInteractionModels(myBasicSprayCloud);
makeSprayParcelStochasticCollisionModels(myBasicSprayCloud);


// Thermo sub-models
makeParcelHeatTransferModels(myBasicSprayCloud);


// Reacting sub-models

makeMyReactingParcelCompositionModels(myBasicSprayCloud);
makeReactingParcelPhaseChangeModels(myBasicSprayCloud);
makeReactingParcelSurfaceFilmModels(myBasicSprayCloud);

*/
// Spray sub-models
//myMakeSprayParcelBreakupModels(myBasicSprayCloud);
//makeParticleForceModelType(DistortedSphereDragForce, myBasicSprayCloud);
myMakeSprayParcelBreakupModels(basicSprayCloud);
//makeParticleForceModelType(myDistortedSphereDragForce, basicSprayCloud);


//makeSprayParcelAtomizationModels(basicSprayCloud);
//makeSprayParcelBreakupModels(myBasicSprayCloud);

// makeParticleForceModelType(myDistortedSphereDragForce, myBasicSprayCloud);
// myMakeSprayParcelDragModels(myBasicSprayCloud);


// MPPIC sub-models
/*
makeMPPICParcelDampingModels(myBasicSprayCloud);
makeMPPICParcelIsotropyModels(myBasicSprayCloud);
makeMPPICParcelPackingModels(myBasicSprayCloud);
*/

// ************************************************************************* //
