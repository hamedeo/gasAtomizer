/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "densitySpherePostProcess.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(densitySpherePostProcess, 0);
    addToRunTimeSelectionTable(functionObject, densitySpherePostProcess, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::densitySpherePostProcess::output
(
    const word& fieldName,
    const word& outputName,
    const scalar& shapePreservation,
    const scalar& volumeConservation
)
{
    OFstream& file = this->file();
    
    writeCurrentTime(file);
    
    writeTabbed(file, fieldName);
    
    file << token::TAB << shapePreservation;
    file << token::TAB << volumeConservation;
    
    file << endl;
    
    Log << " value of shapePreservation is " << shapePreservation << endl;
    Log << " value of volumeConservation is " << volumeConservation << endl;

    // Write state/results information
    word nameStr('(' + outputName + ')');
    this->setResult("value" + nameStr, shapePreservation);
    this->setResult("value" + nameStr, volumeConservation);
}

void Foam::functionObjects::densitySpherePostProcess::writeFileHeader(Ostream& os)
{
    if (!fieldSet_.updateSelection())
    {
        return;
    }

    if (writtenHeader_)
    {
        writeBreak(file());
    }
    else
    {
        writeHeader(os, "Comparison to initial alpha field");
    }

    writeCommented(os, "Time");
    
    writeTabbed(os, "field");
    writeTabbed(os, "shapePreservation");
    writeTabbed(os, "volumeConservation");

    os  << endl;

    writtenHeader_ = true;
}

Foam::scalar Foam::functionObjects::densitySpherePostProcess::calcShapePreservation
(
            const scalarField& cellVolume,
            const volScalarField& alphaFieldZeroTime,
            const volScalarField& alphaFieldCurrTime
)
{
	// Calculate the volume summations for the shape preservation
        
        // \sum_i V_i |\alpha_i(t) - \alpha_i(t=0)|
        scalar alphaDifferenceSum = 0;
        
        // \sum_i V_i \alpha_i(t=0)
        scalar alphaZeroTimeSum = 0;
        
        // Loop over all cells to make the calculation
        forAll(alphaFieldZeroTime, cellI)
        {
		const scalar alphaDifference = mag(alphaFieldCurrTime[cellI] - alphaFieldZeroTime[cellI]);
		
		alphaDifferenceSum += (cellVolume[cellI] * alphaDifference);
		
		alphaZeroTimeSum += (cellVolume[cellI] * alphaFieldZeroTime[cellI]);
        }
        
        return alphaDifferenceSum / alphaZeroTimeSum;
}

Foam::scalar Foam::functionObjects::densitySpherePostProcess::calcVolumeConservation
(
            const scalarField& cellVolume,
            const volScalarField& alphaFieldZeroTime,
            const volScalarField& alphaFieldCurrTime
)
{
	// Calculate the volume summations for the volume conservation
	
	// \sum_i \alpha_i(t=0) V_i (Field at zero time)
	scalar zeroTimeSum = 0;
	
	// \sum_i \alpha_i(t) V_i (Field at current time)
	scalar currTimeSum = 0;
	
	// Loop over all cells to make the calculation
	forAll(alphaFieldZeroTime, cellI)
	{
		zeroTimeSum += (alphaFieldZeroTime[cellI] * cellVolume[cellI]);
		
		currTimeSum += (alphaFieldCurrTime[cellI] * cellVolume[cellI]);
	}
	
	return (currTimeSum - zeroTimeSum) / zeroTimeSum;
	
}

void Foam::functionObjects::densitySpherePostProcess::calcDensitySpherePostProcess(const word& fieldName)
{
    // If the field exists as a scalarfield, calculate the shape preservation
    if (obr_.foundObject<volScalarField>(fieldName))
    {
        // Get the cell volumes 
        // FIX: check for dynamic meshes
        const scalarField& cellVolume = mesh_.V();
        
        // The volume fraction field at time zero       
        const volScalarField alphaFieldZeroTime
        (
        	IOobject
        	(
        		fieldName,
        		time_.timeName(0, 1e-12),
        		mesh_,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE,
        		false
        	),
        	mesh_
        ); 
        
        // The volume fraction field at the current time        
        const volScalarField alphaFieldCurrTime
        (
        	IOobject
        	(
        		fieldName,
        		time_.timeName(),
        		mesh_,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE,
        		false
        	),
        	mesh_
        ); 
        
        scalar shapePreservation = calcShapePreservation(cellVolume, alphaFieldZeroTime, alphaFieldCurrTime);
        
        scalar volumeConservation = calcVolumeConservation(cellVolume, alphaFieldZeroTime, alphaFieldCurrTime);
        
        output
    	(
        	alphaFieldZeroTime.name(),
        	fieldName,
        	shapePreservation,
        	volumeConservation
    	);
    }
    
    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::densitySpherePostProcess::densitySpherePostProcess
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    fieldSet_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::densitySpherePostProcess::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::densitySpherePostProcess::execute()
{
    return true;
}


bool Foam::functionObjects::densitySpherePostProcess::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;

    for (const word& fieldName : fieldSet_.selectionNames())
    {
        calcDensitySpherePostProcess(fieldName);
        
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
