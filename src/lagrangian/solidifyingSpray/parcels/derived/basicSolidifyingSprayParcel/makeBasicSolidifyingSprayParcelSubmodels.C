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

#include "basicSolidifyingSprayCloud.H"

#include "makeReactingParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeThermoParcelTurbulenceForces.H" // add turbulence variant
#include "makeParcelDispersionModels.H"
#include "makeParcelTurbulenceDispersionModels.H" // add turbulence variant
#include "makeSolidifyingSprayParcelInjectionModels.H" // Spray variant
#include "makeParcelPatchInteractionModels.H"
#include "makeSolidifyingSprayParcelStochasticCollisionModels.H" // Spray variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeMyReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Spray
#include "DistortedSphereDragForce.H"
#include "makeSolidifyingSprayParcelAtomizationModels.H"
#include "makeSolidifyingSprayParcelBreakupModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactingParcelCloudFunctionObjects(basicSolidifyingSprayCloud);

// Kinematic sub-models
makeThermoParcelForces(basicSolidifyingSprayCloud);
makeThermoParcelTurbulenceForces(basicSolidifyingSprayCloud);
makeParcelDispersionModels(basicSolidifyingSprayCloud);
makeParcelTurbulenceDispersionModels(basicSolidifyingSprayCloud);
makeSolidifyingSprayParcelInjectionModels(basicSolidifyingSprayCloud);
makeParcelPatchInteractionModels(basicSolidifyingSprayCloud);
makeSolidifyingSprayParcelStochasticCollisionModels(basicSolidifyingSprayCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicSolidifyingSprayCloud);

// Reacting sub-models
makeMyReactingParcelCompositionModels(basicSolidifyingSprayCloud);
makeReactingParcelPhaseChangeModels(basicSolidifyingSprayCloud);
makeReactingParcelSurfaceFilmModels(basicSolidifyingSprayCloud);

// Spray sub-models
makeParticleForceModelType(DistortedSphereDragForce, basicSolidifyingSprayCloud);
makeSolidifyingSprayParcelAtomizationModels(basicSolidifyingSprayCloud);
makeSolidifyingSprayParcelBreakupModels(basicSolidifyingSprayCloud);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicSolidifyingSprayCloud);
makeMPPICParcelIsotropyModels(basicSolidifyingSprayCloud);
makeMPPICParcelPackingModels(basicSolidifyingSprayCloud);

// ************************************************************************* //
