/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "KelvinHelmholtz.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KelvinHelmholtz<CloudType>::KelvinHelmholtz
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName, true),
    b0_(0.61),
    b1_(1.73),//40.0
    cTau_(1.0),
    cRT_(0.1),
    msLimit_(0.03),
    weberLimit_(6.0)
{
	// Pout<< " Test KH Constructor 1 is used " << endl;
    if (!this->defaultCoeffs(true))
    {
        this->coeffDict().readEntry("B0", b0_);
        this->coeffDict().readEntry("B1", b1_);
        this->coeffDict().readEntry("Ctau", cTau_);
        this->coeffDict().readEntry("CRT", cRT_);
        this->coeffDict().readEntry("msLimit", msLimit_);
        this->coeffDict().readEntry("WeberLimit", weberLimit_);
    }
}


template<class CloudType>
Foam::KelvinHelmholtz<CloudType>::KelvinHelmholtz(const KelvinHelmholtz<CloudType>& bum)
:
    BreakupModel<CloudType>(bum),
    b0_(bum.b0_),
    b1_(bum.b1_),
    cTau_(bum.cTau_),
    cRT_(bum.cRT_),
    msLimit_(bum.msLimit_),
    weberLimit_(bum.weberLimit_)
{        // Pout<< " Test KH Constructor 2 is used " << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::KelvinHelmholtz<CloudType>::~KelvinHelmholtz()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::KelvinHelmholtz<CloudType>::update
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
)
{
    bool addParcel = false;

    const scalar averageParcelMass = this->owner().averageParcelMass();

    scalar r = 0.5*d;
    scalar d3 = pow3(d);
    scalar d03 = pow3(d0);

    scalar rhopi6 = rho*constant::mathematical::pi/6.0;
    scalar mass = nParticle*d3*rhopi6;
    scalar mass0 = nParticle*d03*rhopi6;

    // Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;

    scalar weGas = rhoc*sqr(Urmag)*r/sigma;
    scalar weLiquid = rho*sqr(Urmag)*r/sigma;

    // Note: Reitz is using radius instead of diameter for Re-number
    //scalar reLiquid = rho*Urmag*d/mu; 	//Use diameter for the Reynolds number of the liquid
    scalar reLiquid = rho*Urmag*r/mu;
    scalar ohnesorge = sqrt(weLiquid)/(reLiquid + VSMALL);
    scalar taylor = ohnesorge*sqrt(weGas);

    // frequency of the fastest growing KH-wave
    scalar omegaKH =
        (0.34 + 0.38*pow(weGas, 1.5))
       /((1.0 + ohnesorge)*(1.0 + 1.4*pow(taylor, 0.6)))
       *sqrt(sigma/(rho*pow3(r)));

    // corresponding KH wave-length.
    scalar lambdaKH =
        9.02
       *r
       *(1.0 + 0.45*sqrt(ohnesorge))
       *(1.0 + 0.4*pow(taylor, 0.7))
       /pow(1.0 + 0.87*pow(weGas, 1.67), 0.6);//0.865

    // characteristic Kelvin-Helmholtz breakup time
    scalar tauKH = 3.726*b1_*r/(omegaKH*lambdaKH);

    // stable KH diameter
    scalar dc = 2.0*b0_*lambdaKH;

    Pout << __FILE__ << ": " << __LINE__ << "Urmag: " << Urmag << endl;
    Pout << __FILE__ << ": " << __LINE__ << "rhoc: " <<  rhoc << endl;
    Pout << __FILE__ << ": " << __LINE__ << "sigma: " <<  sigma << endl;
    Pout << __FILE__ << ": " << __LINE__ << "dc: " <<  dc << endl;
    Pout << __FILE__ << ": " << __LINE__ << "d: " <<  d << endl;
    Pout << __FILE__ << ": " << __LINE__ << "weGas: " <<  weGas << endl;
    Pout << __FILE__ << ": " << __LINE__ << "weberLimit_: " <<  weberLimit_ << endl;

    if (dc < d)
    {
        // no breakup below Weber = 12
        if (weGas > weberLimit_)
        {
            scalar fraction = dt/tauKH;

            // reduce the diameter according to the rate-equation
            d = (fraction*dc + d)/(1.0 + fraction);
            // d = fraction*dc + d*(1.0 - fraction);

            // addParcel = true;

            // scalar ms0 = rho*pow3(dc)*mathematicalConstant::pi/6.0;
            scalar ms0 = mass0*(1.0 - pow3(d/d0));
            ms += ms0;

            if (ms/averageParcelMass > msLimit_)
            {
                // Correct evaluation of the number of child droplets and the
                // diameter of parcel droplets after breakup
                // Solution of cubic equation for the diameter of the parent
                // drops after breakup, see Eq. 18 in
                // Patterson & Reitz, SAE 980131
                bool br3 = true;
                scalar ae3 = 1.0;
                scalar be3 = -dc;
                scalar ce3 = 0.0;
                scalar de3 = d*d*(dc - d);
                scalar qe3 =
                    pow3(be3/(3.0*ae3)) - be3*ce3/(6.0*ae3*ae3) + de3/(2.0*ae3);
                scalar pe3 = (3.0*ae3*ce3 - be3*be3)/(9.0*ae3*ae3);
                scalar D3 = qe3*qe3 + pe3*pe3*pe3;

                if (D3 < 0) br3 = false;

                if (br3)
                {
                    D3 = sqrt(D3);
                    scalar ue3 = cbrt(-qe3 + D3);
                    scalar ve3 = cbrt(-qe3 - D3);
                    scalar dParenDrops = ue3 + ve3 - be3/3.;
                    scalar mc = nParticle*(pow3(d) - pow3(dParenDrops));
                    scalar nChildDrops = mc/pow3(dc);

                    if (nChildDrops >= nParticle)
                    {
                        addParcel = true;
                        d = dParenDrops;
                        ms = 0.0;
                        dChild = dc;
                        massChild = mc*rhopi6;

                        // reduce the parent mass by reducing nParticle
                        mass -= massChild;
                    }
                }
            }
        }
    }
 /*   else if (KHindex < 0.5)
    {
        // Case of larger drops after breakup (Reitz, Atomization & Spray
        // Technology 3 (1987) 309-337, p.322) pIndKH() should be introduced
        Pout << __FILE__ << ": " << __LINE__ << " KHindex is set to 1.0:" << KHindex << endl;
        scalar lengthScale =
            min(lambdaKH, constant::mathematical::twoPi*Urmag/omegaKH);
        scalar diameterLargerDrop = cbrt(1.5*d*d*lengthScale);
        d = diameterLargerDrop;
        ms = 0.0;
        KHindex = 1.0;
    }
*/
    // correct the number of parcels in parent
    scalar massDrop = pow3(d)*rhopi6;
    nParticle = mass/massDrop;

    return addParcel;
}


// ************************************************************************* //
