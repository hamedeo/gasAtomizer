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

#include "SolidifyingSprayParcel.H"
#include "BreakupModel.H"
#include "CompositionModel.H"
#include "AtomizationModel.H"
#include "mathematicalConstants.H"  // added
// #include "PhaseChangeModel.H"  // added


using namespace Foam::constant::mathematical;  // added
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>  // added
const Foam::label Foam::SolidifyingSprayParcel<ParcelType>::GAS(0);

template<class ParcelType>  // added
const Foam::label Foam::SolidifyingSprayParcel<ParcelType>::LIQ(1);

template<class ParcelType>  // added
const Foam::label Foam::SolidifyingSprayParcel<ParcelType>::SLD(2);

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //
/*
template<class ParcelType>
Foam::scalar Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
    scalarField& Y
) const
{
    scalar mass1 = mass0 - sum(dMass);

    // only update the mass fractions if the new particle mass is finite
    if (mass1 > ROOTVSMALL)
    {
        forAll(Y, i)
        {
            Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
        }
    }

    return mass1;
}
*/

template<class ParcelType>
bool Foam::SolidifyingSprayParcel<ParcelType>::Solidification
(
    const scalar T,
    const scalar Tm0 // Add Tm0 as a function parameter
    // scalarField& this->YLiquid_,
    // scalarField& this->YSolid_
) 
{
    scalarField& YMix = this->Y();

    // Introduce a flag to track whether the condition has been entered
    static bool conditionEntered = false;

    // Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    if (T <= Tm0 && !conditionEntered) 
    {
        // Pout << __FILE__ << ": " << __LINE__ << ": " << "Switch is activated." << endl;
        //forAll(YLIQ, i) {
            this->YSolid_[0] = 1;
            this->YLiquid_[0] = 0;
            YMix[1] = 0;
            YMix[2] = 1;
            // composition.YMixture0() = Ymix;
        Pout << __FILE__ << ": " << __LINE__ << ": " << "Switch is applied. YMix: " << YMix << endl;
        // Pout << __FILE__ << ": " << __LINE__ << ": " << "YLiquid_ : " << YLiquid_ << endl;
        // Pout << __FILE__ << ": " << __LINE__ << ": " << "YSolid_ : " << YSolid_ << endl;
        // Pout << __FILE__ << " : " << __LINE__ << " => YMixture0_: " << YMixture0_ << endl;
        conditionEntered = true;
        //}
        return true;
    }
    return false;
}


template<class ParcelType>  // added
template<class TrackCloudType>
Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::CpEff    // also 4 HsEff LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::HsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().L(idS, YSolid_, p, T);
}

/*
template<class ParcelType>  // added
Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, ROOTVSMALL);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = massSolid/massNew;

    scalar Ytotal = sum(YMix);

    YMix[GAS] /= Ytotal;
    YMix[LIQ] /= Ytotal;
    YMix[SLD] /= Ytotal;

    return massNew;
}

*/

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const auto& composition = cloud.composition();
    // const auto& liquids = composition.liquids();
    // added <<
    // Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    //typedef typename TrackCloudType::reactingCloudType reactingCloudType;
    //const CompositionModel<reactingCloudType>& composition =
    //    cloud.composition();


    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // const scalar np0 = this->nParticle_;
    // const scalar d0 = this->d_;
    // const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    // const scalar mass0 = this->mass();
    // const scalar rho0 = this->rho_;


    const scalar pc = td.pc();

    const scalarField& YMix = this->Y_;
    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();
    // const label YMx() = composition.YMixture0();
    // YMixture0() = this->composition.YMixture0();

    YGas() = composition.Y0(idG);     // added
    YLiquid() = composition.Y0(idL);  // added
    YSolid() = composition.Y0(idS);   // added

    const auto& liquids = composition.liquids();
    // const auto& solids = composition.solids();

    // Pout << __FILE__ << " : " << __LINE__ << " => this YSolid_ " << this->YSolid_[0] << endl;
    // Pout << __FILE__ << " : " << __LINE__ << " => this YLiquid_ " << this->YLiquid_[0] << endl;


    // added >>
    Pout << __FILE__ << " : " << __LINE__ << " => YLiquid_ " << YLiquid_ << endl;
    // Pout << __FILE__ << " : " << __LINE__ << " => YLiquid_[0] " << YLiquid_[0] << endl;
    Pout << __FILE__ << " : " << __LINE__ << " => YSolid_ " << YSolid_ << endl;
    // Pout << __FILE__ << " : " << __LINE__ << " => YSolid_[0] " << YSolid_[0] << endl;
    Pout << __FILE__ << " : " << __LINE__ << " => YMixture0_: " << composition.YMixture0() << ", this->Y(): " << this->Y() << endl;

    // Check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        // Liquid core parcels should not experience coupled forces
        cloud.forces().setCalcCoupled(false);
    }
    // Get old mixture composition
    scalarField X0(liquids.X(this->YLiquid()));

    const scalar Tm0 = 0; // liquids.Tm(X0); // 1811; // added

    // Check if we have critical or boiling conditions
    scalar TMax = liquids.Tc(X0);
    // const scalar T0 = this->T(); // Redeclaration
    const scalar pc0 = td.pc();
    if (liquids.pv(pc0, T0, X0) >= pc0*0.999)
    {    
        // Set TMax to boiling temperature
        TMax = liquids.pvInvert(pc0, X0);
    }

    // Set the maximum temperature limit
    cloud.constProps().setTMax(TMax);

    // Store the parcel properties
    // this->Cp() = liquids.Cp(pc0, T0, X0);
    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);
    // this->sigma_ = this->cloud.composition().sigma(idC, YGas_, p, T);
    sigma_ = liquids.sigma(pc0, T0, X0);    //- Liquid surface tension
    const scalar rho0 = liquids.rho(pc0, T0, X0);    // Redeclaration
    this->rho() = rho0;
    const scalar mass0 = this->mass();    // Redeclaration
    mu_ = liquids.mu(pc0, T0, X0);  //- Liquid dynamic viscosity

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // ParcelType::calc(cloud, td, dt); ///////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;

    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);

    // Sources
    // ~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Phase change
    // ~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), Zero);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), Zero);

    // Calc mass and enthalpy transfer due to phase change
    ParcelType::calcPhaseChange
    (
        cloud,
        td,
        dt,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        rho0,
        0,
        1.0,
        YLiquid_,
        scalarField(),
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scalarField dMass(dMassPC);
    scalar mass1 = ParcelType::updateMassFraction(mass0, dMass, YLiquid_);

    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);
    // this->Cp() = liquids.Cp(pc0, T0, X0);

    if
    (
        cloud.constProps().volUpdateType()
    == constantProperties::volumeUpdateType::mUndefined
    )
    {
        // Update particle density or diameter
        if (cloud.constProps().constantVolume())
        {
            this->rho_ = mass1/this->volume();
        }
        else
        {
            this->d_ = cbrt(mass1/this->rho_*6/pi);
        }
    }
    else
    {
        switch (cloud.constProps().volUpdateType())
        {
            case constantProperties::volumeUpdateType::mConstRho :
            {
                this->d_ = cbrt(mass1/this->rho_*6/pi);
                break;
            }
            case constantProperties::volumeUpdateType::mConstVol :
            {
                this->rho_ = mass1/this->volume();
                break;
            }
            case constantProperties::volumeUpdateType::mUpdateRhoAndVol :
            {
                scalar deltaVol =
                    ParcelType::updatedDeltaVolume
                    (
                        cloud,
                        dMass,
                        td.pc(),
                        T0
                    );
                this->rho_ = mass1/(this->volume() - deltaVol);
                this->d_ = cbrt(mass1/this->rho_*6/pi);
                break;
            }
        }
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(YLiquid_, i)
            {
                scalar dmi = dm*YLiquid_[i];
                label gid = composition.localToCarrierId(LIQ, i);
                // scalar hs = composition.carrier().Hs(gid, td.pc(), T0);

                cloud.rhoTrans(gid)[this->cell()] += dmi;
                // cloud.hsTrans()[this->cell()] += dmi*hs;
            }
            forAll(YSolid_, i)
            {
                label gid = composition.localToCarrierId(SLD, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[SLD]*YSolid_[i];
            }

            cloud.UTrans()[this->cell()] += dm*U0;

            cloud.hsTrans()[this->cell()] +=
                dm*HsEff(cloud, td, pc, T0, idG, idL, idS);

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    ParcelType::correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(rhos, U0, td.Uc(), this->d(), mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );

    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);
    // this->CpEff_ = composition.Cp(0, Y_, td.pc(), T0);
    // this->Cp() = liquids.Cp(pc0, T0, X0);
    
    if (YMix[idL] != 0) // added condition
    {
        (void)Solidification(this->T_, Tm0); //, YLiquid_, YSolid_);
    }

    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (cloud.solution().coupled())
    {
        // Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl; 
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        /*forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];
            label gid = composition.localToCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }
        forAll(YSolid_, i)
        {
            scalar dm = np0*dMassSolid[i];
            label gid = composition.localToCarrierId(SLD, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }
        forAll(dMass, i)
        {
            scalar dm = np0*dMass[i];
            scalar hs = composition.carrier().Hs(i, pc, T0);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
            scalar dm = np0*dMass[i];
            label gid = composition.localToCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, td.pc(), T0); //should I ???

            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
            cloud.hsTrans()[this->cell()] += dm*hs;
        }*/
        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeff()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {    
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;

        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////

    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);
    // if ()
    // {
    if (td.keepParticle && YMix[1] == 1)
    {
        // Reduce the stripped parcel mass due to evaporation
        // assuming the number of particles remains unchanged
        this->ms() -= this->ms()*(mass0 - this->mass())/mass0;

        // Update Cp, sigma, density and diameter due to change in temperature
        // and/or composition
        scalar T1 = this->T();
        scalarField X1(liquids.X(this->YLiquid()));

        this->Cp() = liquids.Cp(td.pc(), T1, X1);

        sigma_ = liquids.sigma(td.pc(), T1, X1);

        scalar rho1 = liquids.rho(td.pc(), T1, X1);
        this->rho() = rho1;

        mu_ = liquids.mu(td.pc(), T1, X1);

        scalar d1 = this->d()*cbrt(rho0/rho1);
        this->d() = d1;

        if (liquidCore() > 0.5)
        {
            calcAtomization(cloud, td, dt);

            // Preserve the total mass/volume by increasing the number of
            // particles in parcels due to breakup
            scalar d2 = this->d();
            this->nParticle() *= pow3(d1/d2);
        }
        else
        {
            // Pout << __FILE__ << ": " << __LINE__ << " is reached" << endl;
            // if (YLiquid_[0] == 1) calcBreakup(cloud, td, dt);
            calcBreakup(cloud, td, dt);
        }
    }
    // }

    // Restore coupled forces
    cloud.forces().setCalcCoupled(true);

}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::calcAtomization
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    const auto& atomization = cloud.atomization();

    if (!atomization.active())
    {
        return;
    }

    const auto& composition = cloud.composition();
    const auto& liquids = composition.liquids();

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();
    scalar R = RR/Wc;
    scalar Tav = atomization.Taverage(this->T(), td.Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = td.pc()/(R*Tav);

    scalar soi = cloud.injectors().timeStart();
    scalar currentTime = cloud.db().time().value();
    const vector pos(this->position());
    const vector& injectionPos = this->position0();

    // Disregard the continuous phase when calculating the relative velocity
    // (in line with the deactivated coupled assumption)
    scalar Urel = mag(this->U());

    scalar t0 = max(0.0, currentTime - this->age() - soi);
    scalar t1 = min(t0 + dt, cloud.injectors().timeEnd() - soi);

    // This should be the vol flow rate from when the parcel was injected
    scalar volFlowRate = cloud.injectors().volumeToInject(t0, t1)/dt;

    scalar chi = 0.0;
    if (atomization.calcChi())
    {
        chi = this->chi(cloud, td, liquids.X(this->Y()));
    }

    atomization.update
    (
        dt,
        this->d(),
        this->liquidCore(),
        this->tc(),
        this->rho(),
        mu_,
        sigma_,
        volFlowRate,
        rhoAv,
        Urel,
        pos,
        injectionPos,
        cloud.pAmbient(),
        chi,
        cloud.rndGen()
    );
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::calcBreakup
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
    // const scalarField& YLiquidEff,
)
{
    Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    auto& breakup = cloud.breakup();

    if (!breakup.active())
    {
        return;
    }

    if (breakup.solveOscillationEq())
    {
        solveTABEq(cloud, td, dt);
    }

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = td.rhoc()*RR*td.Tc()/td.pc();
    scalar R = RR/Wc;
    scalar Tav = cloud.atomization().Taverage(this->T(), td.Tc());

    // Calculate average gas density based on average temperature
    scalar rhoAv = td.pc()/(R*Tav);
    scalar muAv = td.muc();
    vector Urel = this->U() - td.Uc();
    Pout << __FILE__ << ": " << __LINE__ << "this->U(): " <<  this->U() << endl;
    Pout << __FILE__ << ": " << __LINE__ << "td.Uc(): " <<  td.Uc() << endl;
    Pout << __FILE__ << ": " << __LINE__ << "Urel: " <<  Urel << endl;
    scalar Urmag = mag(Urel);
    scalar Re = this->Re(rhoAv, this->U(), td.Uc(), this->d(), muAv);

    const typename TrackCloudType::parcelType& p =
        static_cast<const typename TrackCloudType::parcelType&>(*this);
    typename TrackCloudType::parcelType::trackingData& ttd =
        static_cast<typename TrackCloudType::parcelType::trackingData&>(td);
    const scalar mass = p.mass();
    const typename TrackCloudType::forceType& forces = cloud.forces();
    const forceSuSp Fcp = forces.calcCoupled(p, ttd, dt, mass, Re, muAv);
    const forceSuSp Fncp = forces.calcNonCoupled(p, ttd, dt, mass, Re, muAv);
    this->tMom() = mass/(Fcp.Sp() + Fncp.Sp() + ROOTVSMALL);

    const vector g = cloud.g().value();

    scalar parcelMassChild = 0.0;
    scalar dChild = 0.0;
    if
    (
        breakup.update
        (
            dt,
            g,
            this->d(),
            this->tc(),
            this->ms(),
            this->nParticle(),
            this->KHindex(),
            this->y(),
            this->yDot(),
            this->d0(),
            this->rho(),
            mu_,
            sigma_,
            this->U(),
            rhoAv,
            muAv,
            Urel,
            Urmag,
            this->tMom(),
            dChild,
            parcelMassChild
        )
    )
    {
        scalar Re = rhoAv*Urmag*dChild/muAv;

        // Add child parcel as copy of parent
        SolidifyingSprayParcel<ParcelType>* child = new SolidifyingSprayParcel<ParcelType>(*this);
        child->origId() = this->getNewParticleID();
        child->d() = dChild;
        child->d0() = dChild;
        const scalar massChild = child->mass();
        child->mass0() = massChild;
        child->nParticle() = parcelMassChild/massChild;

        const forceSuSp Fcp =
            forces.calcCoupled(*child, ttd, dt, massChild, Re, muAv);
        const forceSuSp Fncp =
            forces.calcNonCoupled(*child, ttd, dt, massChild, Re, muAv);

        child->age() = 0.0;
        child->liquidCore() = 0.0;
        child->KHindex() = 1.0;
        child->y() = cloud.breakup().y0();
        child->yDot() = cloud.breakup().yDot0();
        child->tc() = 0.0;
        child->ms() = -GREAT;
        child->injector() = this->injector();
        child->tMom() = massChild/(Fcp.Sp() + Fncp.Sp());
        child->user() = 0.0;
        child->calcDispersion(cloud, td, dt);

        cloud.addParticle(child);
    }
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::SolidifyingSprayParcel<ParcelType>::chi
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalarField& X
) const
{
    // Modifications to take account of the flash boiling on primary break-up

    const auto& composition = cloud.composition();
    const auto& liquids = composition.liquids();

    scalar chi = 0.0;
    scalar T0 = this->T();
    scalar p0 = td.pc();
    scalar pAmb = cloud.pAmbient();

    scalar pv = liquids.pv(p0, T0, X);

    forAll(liquids, i)
    {
        if (pv >= 0.999*pAmb)
        {
            // Liquid is boiling - calc boiling temperature

            const liquidProperties& liq = liquids.properties()[i];
            scalar TBoil = liq.pvInvert(p0);

            scalar hl = liq.hl(pAmb, TBoil);
            scalar iTp = liq.h(pAmb, T0) - pAmb/liq.rho(pAmb, T0);
            scalar iTb = liq.h(pAmb, TBoil) - pAmb/liq.rho(pAmb, TBoil);

            chi += X[i]*(iTp - iTb)/hl;
        }
    }

    chi = min(1.0, max(chi, 0.0));

    return chi;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::SolidifyingSprayParcel<ParcelType>::solveTABEq
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Pout << __FILE__ << ": " << __LINE__ << ": " <<  __FUNCTION__<< " is reached" << endl;
    const scalar& TABCmu = cloud.breakup().TABCmu();
    const scalar& TABtwoWeCrit = cloud.breakup().TABtwoWeCrit();
    const scalar& TABComega = cloud.breakup().TABComega();

    scalar r = 0.5*this->d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    // Inverse of characteristic viscous damping time
    scalar rtd = 0.5*TABCmu*mu_/(this->rho()*r2);

    // Oscillation frequency (squared)
    scalar omega2 = TABComega*sigma_/(this->rho()*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar We =
            this->We(td.rhoc(), this->U(), td.Uc(), r, sigma_)/TABtwoWeCrit;

        // Initial values for y and yDot
        scalar y0 = this->y() - We;
        scalar yDot0 = this->yDot() + y0*rtd;

        // Update distortion parameters
        scalar c = cos(omega*dt);
        scalar s = sin(omega*dt);
        scalar e = exp(-rtd*dt);

        this->y() = We + e*(y0*c + (yDot0/omega)*s);
        this->yDot() = (We - this->y())*rtd + e*(yDot0*c - omega*y0*s);
    }
    else
    {
        // Reset distortion parameters
        this->y() = 0;
        this->yDot() = 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel(const SolidifyingSprayParcel<ParcelType>& p)
:
    ParcelType(p),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_),
    YGas_(p.YGas_), //added
    YLiquid_(p.YLiquid_), //added
    YSolid_(p.YSolid_) //added
{}


template<class ParcelType>
Foam::SolidifyingSprayParcel<ParcelType>::SolidifyingSprayParcel
(
    const SolidifyingSprayParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    d0_(p.d0_),
    position0_(p.position0_),
    sigma_(p.sigma_),
    mu_(p.mu_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_),
    YGas_(p.YGas_), //added
    YLiquid_(p.YLiquid_), //added
    YSolid_(p.YSolid_) //added
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SolidifyingSprayParcelIO.C"


// ************************************************************************* //
