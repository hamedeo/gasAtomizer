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

inline Foam::scalar Foam::Fe::rho(scalar p, scalar T) const
{
    //Info << "Call for density f(p,T)" << endl;
    return rho_.f(p, T);
}


inline Foam::scalar Foam::Fe::pv(scalar p, scalar T) const
{
    //Info << "Call for vapor pressure f(p,T)" << endl;
    return pv_.f(p, T);
}


inline Foam::scalar Foam::Fe::hl(scalar p, scalar T) const
{
    //Info << "Call for heat of vaporisation f(p,T)" << endl;
    return hl_.f(p, T);
}


inline Foam::scalar Foam::Fe::Cp(scalar p, scalar T) const
{
    //Info << "Call for heat capacity f(p,T)" << endl;
    return Cp_.f(p, T);
}


inline Foam::scalar Foam::Fe::h(scalar p, scalar T) const
{
    //Info << "Call for liquid enthalpy f(p,T)" << endl;
    return h_.f(p, T);
}


inline Foam::scalar Foam::Fe::Cpg(scalar p, scalar T) const
{
    //Info << "Call for gas heat capacity f(p,T)" << endl;
    return Cpg_.f(p, T);
}


inline Foam::scalar Foam::Fe::B(scalar p, scalar T) const
{
    //Info << "Call for virial coefficient f(p,T)" << endl;
    return B_.f(p, T);
}


inline Foam::scalar Foam::Fe::mu(scalar p, scalar T) const
{
    //Info << "Call for viscosity f(p,T)" << endl;
    return mu_.f(p, T);
}


inline Foam::scalar Foam::Fe::mug(scalar p, scalar T) const
{
    //Info << "Call for gas viscosity f(p,T)" << endl;
    return mug_.f(p, T);
}


inline Foam::scalar Foam::Fe::kappa(scalar p, scalar T) const
{
    //Info << "Call for thermal conductivity f(p,T)" << endl;
    return kappa_.f(p, T);
}


inline Foam::scalar Foam::Fe::kappag(scalar p, scalar T) const
{
    //Info << "Call for vapor thermal conductivity f(p,T)" << endl;
    return kappag_.f(p, T);
}


inline Foam::scalar Foam::Fe::sigma(scalar p, scalar T) const
{
    //Info << "Call for surface tension f(p,T)" << endl;
    return sigma_.f(p, T);
}


inline Foam::scalar Foam::Fe::D(scalar p, scalar T) const
{
    //Info << "Call for vapor diffusivity f(p,T)" << endl;
    return D_.f(p, T);
}


inline Foam::scalar Foam::Fe::D(scalar p, scalar T, scalar Wb) const
{
    //Info << "Call for vapor diffusivity with binary pair f(p,T)" << endl;
    return D_.f(p, T, Wb);
}


// ************************************************************************* //
