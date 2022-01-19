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

#include "normalizedKineticEnergy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(normalizedKineticEnergy, 0);
    addToRunTimeSelectionTable(functionObject, normalizedKineticEnergy, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::normalizedKineticEnergy::output
(
    const word& fieldName,
    const scalar& normKinEnergy,
    const scalar& Ek_t,
    const scalar& Ek_0
)
{
    OFstream& file = this->file();
    
    writeCurrentTime(file);
    
    writeTabbed(file, fieldName);
    
    file << token::TAB << normKinEnergy;
    file << token::TAB << Ek_t;
    file << token::TAB << Ek_0;
    
    file << endl;
    
    Log << " value of normKinEnergy is " << normKinEnergy << endl;

}

void Foam::functionObjects::normalizedKineticEnergy::writeFileHeader(Ostream& os)
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
        writeHeader(os, "Normalized kinetic energy");
    }

    writeCommented(os, "Time");
    
    writeTabbed(os, "field");
    writeTabbed(os, "normKinEnergy");
    writeTabbed(os, "Ek_t");
    writeTabbed(os, "Ek_0");

    os  << endl;

    writtenHeader_ = true;
}

void Foam::functionObjects::normalizedKineticEnergy::calcNormalizedKineticEnergy(const word& fieldName)
{
    // If the field exists as a scalarfield, calculate the kinetic energy
    if (obr_.foundObject<volVectorField>(fieldName))
    {
        
	// The alpha field at time zero
	const volScalarField alphaFieldZeroTime
	(
		IOobject
		(
			"alpha.liquid",
			time_.timeName(0, 1e-12),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		mesh_
	);
	
	// The alpha field at the current time
	const volScalarField alphaFieldCurrTime
	(
		IOobject
		(
			"alpha.liquid",
			time_.timeName(),
			mesh_,
			IOobject::NO_READ,
			IOobject::NO_WRITE,
			false
		),
		mesh_
	);
	
	// The velocity field at zero time
	const volVectorField velocityFieldZeroTime
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
	
	// The velocity field at the current time
	const volVectorField& velocityFieldCurrTime = lookupObject<volVectorField>(fieldName);
	
	// Get the densities of both phases
        const dictionary& transportPropertiesDict = mesh_.lookupObject<IOdictionary>
        (
        	"transportProperties"
        );
        
        const wordList phaseNamesList
        (
        	transportPropertiesDict.lookup("phases")
        );
        
        const word phase1 = phaseNamesList[0];
        const word phase2 = phaseNamesList[1];
        const dictionary& phaseDict1 = transportPropertiesDict.subDict(phase1);
        const dictionary& phaseDict2 = transportPropertiesDict.subDict(phase2);
        const scalar rho1 = phaseDict1.getScalar("rho");
        const scalar rho2 = phaseDict2.getScalar("rho");
        
        // Calculate the kinetic energy at current and zero time
        scalar Ek_t = 0;
        scalar Ek_0 = 0;
        
        forAll(mesh_.C(), cellI)
        {
        	const scalar rhoCellCurrTime = alphaFieldCurrTime[cellI] * rho1 + (1 - alphaFieldCurrTime[cellI]) * rho2;
        	
        	const scalar rhoCellZeroTime = alphaFieldZeroTime[cellI] * rho1 + (1 - alphaFieldZeroTime[cellI]) * rho2;
        
        	Ek_t += rhoCellCurrTime * mag(velocityFieldCurrTime[cellI]) * mag(velocityFieldCurrTime[cellI]);
        	
        	Ek_0 += rhoCellZeroTime * mag(velocityFieldZeroTime[cellI]) * mag(velocityFieldZeroTime[cellI]);
        }
        
        scalar normKinEnergy = Ek_t / Ek_0;
        
        output
    	(
        	fieldName,
        	normKinEnergy,
        	Ek_t,
        	Ek_0
    	);
    }
    
    
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::normalizedKineticEnergy::normalizedKineticEnergy
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

bool Foam::functionObjects::normalizedKineticEnergy::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    fieldSet_.read(dict);

    return true;
}


bool Foam::functionObjects::normalizedKineticEnergy::execute()
{
    return true;
}


bool Foam::functionObjects::normalizedKineticEnergy::write()
{
    writeFileHeader(file());

    Log << type() << " " << name() <<  " write:" << nl;
    
    //Assign field names
    for (const word& fieldName : fieldSet_.selectionNames())
    {
    	calcNormalizedKineticEnergy(fieldName);

    }

    Log << endl;

    return true;
}


// ************************************************************************* //
