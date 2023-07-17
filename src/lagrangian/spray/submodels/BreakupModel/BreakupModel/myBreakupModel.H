/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::myBreakupModel

Group
    grpLagrangianIntermediateBreakupSubModels

Description
    Templated break-up model class

SourceFiles
    myBreakupModel.C
    myBreakupModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef myBreakupModel_H
#define myBreakupModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class myBreakupModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class myBreakupModel
:
    public CloudSubModelBase<CloudType>
{

protected:

    // Protected data

        bool solveOscillationEq_;

        scalar y0_;
        scalar yDot0_;

        scalar TABComega_;
        scalar TABCmu_;
        scalar TABtwoWeCrit_;


public:

    //- Runtime type information
    TypeName("breakupModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        myBreakupModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        myBreakupModel(CloudType& owner);

        //- Construct from dictionary
        myBreakupModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type,
            bool solveOscillationEq = false
        );

        //- Construct copy
        myBreakupModel(const myBreakupModel<CloudType>& bum);

        //- Construct and return a clone
        virtual autoPtr<myBreakupModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~myBreakupModel();


    //- Selector
    static autoPtr<myBreakupModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        inline bool solveOscillationEq() const
        {
            return solveOscillationEq_;
        }

        inline const scalar& y0() const
        {
            return y0_;
        }

        inline const scalar& yDot0() const
        {
            return yDot0_;
        }

        inline const scalar& TABComega() const
        {
            return TABComega_;
        }

        inline const scalar& TABCmu() const
        {
            return TABCmu_;
        }

        inline const scalar& TABtwoWeCrit() const
        {
            return TABtwoWeCrit_;
        }


    // Member Functions

        //- Update the parcel properties and return true if a child parcel
        //  should be added
        virtual bool update
        (
            const scalar dt,
            const vector& g,
            scalar& d,
            scalar& tc,
            scalar& ms,
            scalar& nParticle,
            scalar& KHindex,
            scalar& y,
            scalar& yDot,
            const scalar d0,
            const scalar rho,
            const scalar mu,
            const scalar sigma,
            const vector& U,
            const scalar rhoc,
            const scalar muc,
            const vector& Urel,
            const scalar Urmag,
            const scalar tMom,
            scalar& dChild,
            scalar& massChild
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define myMakeBreakupModel(CloudType)                                            \
                                                                               \
    typedef Foam::CloudType::mysprayCloudType mysprayCloudType;                    \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::myBreakupModel<mysprayCloudType>,                                    \
        0                                                                      \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            myBreakupModel<mysprayCloudType>,                                      \
            dictionary                                                         \
        );                                                                     \
    }


#define myMakeBreakupModelType(SS, CloudType)                                    \
                                                                               \
    typedef Foam::CloudType::mysprayCloudType mysprayCloudType;                    \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<mysprayCloudType>, 0);          \
                                                                               \
    Foam::myBreakupModel<mysprayCloudType>::                                       \
        adddictionaryConstructorToTable<Foam::SS<mysprayCloudType>>              \
            add##SS##CloudType##mysprayCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myBreakupModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
