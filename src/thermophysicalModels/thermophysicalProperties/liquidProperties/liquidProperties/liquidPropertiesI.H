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

inline Foam::scalar Foam::liquidProperties::limit(const scalar T) const
{
    return T;
}


inline Foam::scalar Foam::liquidProperties::Y() const
{
    return 1;
}


inline Foam::scalar Foam::liquidProperties::Tc() const
{
    return Tc_;
}


inline Foam::scalar Foam::liquidProperties::Pc() const
{
    return Pc_;
}


inline Foam::scalar Foam::liquidProperties::Vc() const
{
    return Vc_;
}


inline Foam::scalar Foam::liquidProperties::Zc() const
{
    return Zc_;
}


inline Foam::scalar Foam::liquidProperties::Tt() const
{
    return Tt_;
}


inline Foam::scalar Foam::liquidProperties::Pt() const
{
    return Pt_;
}


inline Foam::scalar Foam::liquidProperties::Tb() const
{
    return Tb_;
}


inline Foam::scalar Foam::liquidProperties::Tm() const
{
    return Tm_;
}

inline Foam::scalar Foam::liquidProperties::psi(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::CpMCv(scalar p, scalar T) const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::Ha(scalar p, scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::liquidProperties::Hs(scalar p, scalar T) const
{
    return h(p, T);
}


inline Foam::scalar Foam::liquidProperties::Hc() const
{
    return 0;
}


inline Foam::scalar Foam::liquidProperties::alphah(scalar p, scalar T) const
{
    return kappa(p, T)/Cp(p, T);
}


template<class Func>
inline void Foam::liquidProperties::readIfPresent
(
    Func& f,
    const word& name,
    const dictionary& dict
)
{
    if (dict.found(name))
    {
        f = Func(dict.subDict(name));
    }
}


template<class Liquid>
inline void Foam::liquidProperties::readIfPresent
(
    Liquid& l,
    const dictionary& dict
)
{
    l.liquidProperties::readIfPresent(dict);
    readIfPresent(l.rho_, "rho", dict);
    readIfPresent(l.pv_, "pv", dict);
    readIfPresent(l.hl_, "hl", dict);
    readIfPresent(l.Cp_, "Cp", dict);
    readIfPresent(l.h_, "h", dict);
    readIfPresent(l.mu_, "mu", dict);
    readIfPresent(l.kappa_, "K", dict);
    readIfPresent(l.sigma_, "sigma", dict);
}


template<class Liquid>
inline void Foam::liquidProperties::writeData
(
    const Liquid& l,
    Ostream& os
) const
{
    l.liquidProperties::writeData(os); os << nl;
    l.rho_.writeData(os); os << nl;
    l.pv_.writeData(os); os << nl;
    l.hl_.writeData(os); os << nl;
    l.Cp_.writeData(os); os << nl;
    l.h_.writeData(os); os << nl;
    l.mu_.writeData(os); os << nl;
    l.kappa_.writeData(os); os << nl;
    l.sigma_.writeData(os); os << endl;
}


// ************************************************************************* //
