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
    Foam::thermophysicalConstant

Description


\*---------------------------------------------------------------------------*/

#ifndef thermophysicalConstant_H
#define thermophysicalConstant_H

#include "thermophysicalFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermophysicalConstant Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalConstant
:
    public thermophysicalFunction
{
    // Private data

        
        // property value
        scalar value_;

public:

    //- Runtime type information
    TypeName("thermophysicalConstant");


    // Constructors

        //- Construct from components
        thermophysicalConstant
        (
            const scalar value
        );

        //- Construct from dictionary
        thermophysicalConstant(const dictionary& dict);


    // Member Functions

        //- Evaluate the function and return the result
        scalar f(scalar, scalar T) const
        {
            return value_;
        }

        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            os  << value_;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const thermophysicalConstant& f)
        {
            f.writeData(os);
            return os;
        }
        


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
