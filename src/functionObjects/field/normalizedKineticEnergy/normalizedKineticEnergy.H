/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::functionObjects::normalizedKineticEnergy

Group
    grpFieldFunctionObjects

Description
    Computes the values and locations of field minima and maxima.
    These are good indicators of calculation performance, e.g. to confirm that
    predicted results are within expected bounds, or how well a case is
    converging.

    Multiple fields can be processed, where for rank > 0 primitives, e.g.
    vectors and tensors, the extrema can be calculated per component, or by
    magnitude.  In addition, spatial location and local processor index are
    included in the output.

    Operands:
    \table
      Operand        | Type | Location
      input          | -    | -
      output file    | dat  | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | -    | -
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    normalizedKineticEnergy1
    {
        // Mandatory entries (unmodifiable)
        type        normalizedKineticEnergy;
        libs        (myFieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        mode        magnitude;
        fields      (<field1> <field2> ... <fieldN>);

        // Optional entries (runtime modifiable)
        location    true;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: normalizedKineticEnergy             | word |  yes  | -
      libs         | Library name: fieldFunctionObjects | word |  yes  | -
      fields       | List of operand fields             | wordList | yes | -
      location     | Write location of the min/max value | bool | no   | true
      mode   | Calculation mode: magnitude or component | word | no | magnitude
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::normalizedKineticEnergy

SourceFiles
    normalizedKineticEnergy.C
    normalizedKineticEnergyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_normalizedKineticEnergy_H
#define functionObjects_normalizedKineticEnergy_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "vector.H"
#include "scalar.H"
#include "volFieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class normalizedKineticEnergy Declaration
\*---------------------------------------------------------------------------*/

class normalizedKineticEnergy
:
    public fvMeshFunctionObject,
    public writeFile
{

protected:

    // Protected Data

        //- Fields to assess
        volFieldSelection fieldSet_;

    // Protected Member Functions

        //- Helper function to write the output        
        void output
        (
            const word& fieldName,
            const scalar& normKinEnergy,
            const scalar& Ek_t,
            const scalar& Ek_0
        );

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

	//- Post processing function
        void calcNormalizedKineticEnergy
        (
	    const word& fieldName
        );

public:

    //- Runtime type information
    TypeName("normalizedKineticEnergy");


    // Constructors

        //- Construct from Time and dictionary
        normalizedKineticEnergy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        normalizedKineticEnergy(const normalizedKineticEnergy&) = delete;

        //- No copy assignment
        void operator=(const normalizedKineticEnergy&) = delete;


    //- Destructor
    virtual ~normalizedKineticEnergy() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the normalizedKineticEnergy
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
