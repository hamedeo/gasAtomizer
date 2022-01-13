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

#include "volumeConservationVOF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volumeConservationVOF, 0);
    addToRunTimeSelectionTable(functionObject, volumeConservationVOF, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volumeConservationVOF::output
(
    const word& fieldName,
    const word& outputName,
    const scalar& volumeConservation
)
{
    OFstream& file = this->file();
    
    writeCurrentTime(file);
    
    writeTabbed(file, fieldName);
    
    file << token::TAB << volumeConservation;
    
    file << endl;
    
    Log << " value of volumeConservation is " << volumeConservation << endl;

    // Write state/results information
    word nameStr('(' + outputName + ')');
    this->setResult("value" + nameStr, volumeConservation);
}

void Foam::functionObjects::volumeConservationVOF::writeFileHeader(Ostream& os)
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
        writeHeader(os, "Conservation of volume in VOF method");
    }

    writeCommented(os, "Time");
    
    writeTabbed(os, "field");
    writeTabbed(os, "volumeConservation");

    os  << endl;

    writtenHeader_ = true;
}

void Foam::functionObjects::volumeConservationVOF::calcVolumeConservation(const word& fieldName)
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
	
	scalar volumeConservation = (currTimeSum - zeroTimeSum) / zeroTimeSum;
        
        output
    	(
        	alphaFieldZeroTime.name(),
        	fieldName,
        	volumeConservation
    	);
    }
    
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volumeConservationVOF::volumeConservationVOF
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

bool Foam::functionObjects::volumeConservationVOF::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::volumeConservationVOF::execute()
{
    return true;
}


bool Foam::functionObjects::volumeConservationVOF::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;
    
    //Assign field names
    for (const word& fieldName : fieldSet_.selectionNames())
    {
    	calcVolumeConservation(fieldName);

    }

    Log << endl;

    return true;
}


// ************************************************************************* //
