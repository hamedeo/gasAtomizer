/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::SolidifyingSprayCloud<CloudType>&
Foam::SolidifyingSprayCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>   // added
inline const typename CloudType::particleType::constantProperties&
Foam::SolidifyingSprayCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>   // added
inline typename CloudType::particleType::constantProperties&
Foam::SolidifyingSprayCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::AtomizationModel<Foam::SolidifyingSprayCloud<CloudType>>&
Foam::SolidifyingSprayCloud<CloudType>::atomization() const
{
    return atomizationModel_;
}


template<class CloudType>
inline Foam::AtomizationModel<Foam::SolidifyingSprayCloud<CloudType>>&
Foam::SolidifyingSprayCloud<CloudType>::atomization()
{
    return *atomizationModel_;
}


template<class CloudType>
inline const Foam::BreakupModel<Foam::SolidifyingSprayCloud<CloudType>>&
Foam::SolidifyingSprayCloud<CloudType>::breakup() const
{
    return breakupModel_;
}


template<class CloudType>
inline Foam::BreakupModel<Foam::SolidifyingSprayCloud<CloudType>>&
Foam::SolidifyingSprayCloud<CloudType>::breakup()
{
    return *breakupModel_;
}


template<class CloudType>
inline Foam::scalar Foam::SolidifyingSprayCloud<CloudType>::averageParcelMass() const
{
    return averageParcelMass_;
}


template<class CloudType>
inline Foam::scalar Foam::SolidifyingSprayCloud<CloudType>::penetration
(
    const scalar fraction
) const
{
    if ((fraction < 0) || (fraction > 1))
    {
        FatalErrorInFunction
            << "fraction should be in the range 0 < fraction < 1"
            << exit(FatalError);
    }

    scalar distance = 0.0;

    const label nParcel = this->size();
    globalIndex globalParcels(nParcel);
    const label nParcelSum = globalParcels.size();

    if (nParcelSum == 0)
    {
        return distance;
    }

    // lists of parcels mass and distance from initial injection point
    List<List<scalar>> procMass(Pstream::nProcs());
    List<List<scalar>> procDist(Pstream::nProcs());

    List<scalar>& mass = procMass[Pstream::myProcNo()];
    List<scalar>& dist = procDist[Pstream::myProcNo()];

    mass.setSize(nParcel);
    dist.setSize(nParcel);

    label i = 0;
    scalar mSum = 0.0;
    for (const parcelType& p : *this)
    {
        scalar m = p.nParticle()*p.mass();
        scalar d = mag(p.position() - p.position0());
        mSum += m;

        mass[i] = m;
        dist[i] = d;

        ++i;
    }

    // calculate total mass across all processors
    reduce(mSum, sumOp<scalar>());
    Pstream::gatherList(procMass);
    Pstream::gatherList(procDist);

    if (Pstream::master())
    {
        // flatten the mass lists
        List<scalar> allMass(nParcelSum, Zero);
        SortableList<scalar> allDist(nParcelSum, Zero);
        for (const int proci : Pstream::allProcs())
        {
            SubList<scalar>
            (
                allMass,
                globalParcels.localSize(proci),
                globalParcels.offset(proci)
            ) = procMass[proci];

            // flatten the distance list
            SubList<scalar>
            (
                allDist,
                globalParcels.localSize(proci),
                globalParcels.offset(proci)
            ) = procDist[proci];
        }

        // sort allDist distances into ascending order
        // note: allMass masses are left unsorted
        allDist.sort();

        if (nParcelSum > 1)
        {
            const scalar mLimit = fraction*mSum;
            const labelList& indices = allDist.indices();

            if (mLimit > (mSum - allMass[indices.last()]))
            {
                distance = allDist.last();
            }
            else
            {
                // assuming that 'fraction' is generally closer to 1 than 0,
                // loop through in reverse distance order
                const scalar mThreshold = (1.0 - fraction)*mSum;
                scalar mCurrent = 0.0;
                label i0 = 0;

                forAllReverse(indices, i)
                {
                    label indI = indices[i];

                    mCurrent += allMass[indI];

                    if (mCurrent > mThreshold)
                    {
                        i0 = i;
                        break;
                    }
                }

                if (i0 == indices.size() - 1)
                {
                    distance = allDist.last();
                }
                else
                {
                    // linearly interpolate to determine distance
                    scalar alpha = (mCurrent - mThreshold)/allMass[indices[i0]];
                    distance =
                        allDist[i0] + alpha*(allDist[i0+1] - allDist[i0]);
                }
            }
        }
        else
        {
            distance = allDist.first();
        }
    }

    Pstream::scatter(distance);

    return distance;
}


// ************************************************************************* //
