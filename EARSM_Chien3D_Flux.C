/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
	\\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
	 \\/     M anipulation  |
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
#include "ChienKEpsilonEARSM_Flux.H"
#include "wallDist.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
	namespace RASModels
	{

		// * * * * * * * * * * * * Static Data Members * * * * * * * * * * //

		// * * * * * * * * * * Private Member Functions * * * * * * * * * //
		template <class BasicMomentumTransportModel>
		void ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::correctNut()
		{
			this->nut_ == Cmu_ *k_ *tau_;
			this->nut_.correctBoundaryConditions();
			fvConstraints::New(this->mesh_).constrain(this->nut_);
		}
		template <class BasicMomentumTransportModel>
		void ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::correctDt()
		{
			if (fluxType_ == "DalyHarlow")
			{
				Dt_ == -Cc_ *(tensor::I & R_) * tau_;
			}
			else if (fluxType_ == "SugaAbe")
			{
				Dt_ == -Cc_ *(R_ & R_) / k_ *tau_;
			}
			else if (fluxType_ == "WWJ")
			{
				Dt_ == -(1.0 - Cc_) * (Bij_ & (aij_ + 2.0 / 3 * I) * k_) * tau_; // TODO:check
			}
			else
			{
				WarningInFunction
					<< "fluxType are DalyHarlow, SugaAbe, WWJ"
					<< endl;
			}

			Info << "Dt max/min/ave: "
				 << max(mag(Dt_)).value() << ","
				 << min(mag(Dt_)).value() << ","
				 << average(mag(Dt_)).value() << endl;

			Dt_.correctBoundaryConditions();
			fvConstraints::New(this->mesh_).constrain(Dt_);
		}
		template <class BasicMomentumTransportModel>
		tmp<volScalarField> ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::fMu() const
		{
			const volScalarField &y = wallDist(this->mesh_).y();
			volScalarField Rey = pow(k_, 0.5) * (y / this->nu());
			volScalarField yStar = pow(Rey, 2) * 0.003 + pow(Rey, 0.5) * 2.4;
			return scalar(1) - exp(-0.0115 * yStar);
		}
		template <class BasicMomentumTransportModel>
		tmp<volScalarField> ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::f2() const
		{
			return scalar(1) - 0.22 * exp(-sqr(sqr(k_) / (this->nu() * epsilonTilda_ * 6))); // TODO:check
		}
		// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * //
		template <class BasicMomentumTransportModel>
		ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::ChienKEpsilonEARSM_Flux(
			const alphaField &alpha,
			const rhoField &rho,
			const volVectorField &U,
			const surfaceScalarField &alphaRhoPhi,
			const surfaceScalarField &phi,
			const transportModel &transport,
			const word &type) : eddyViscosity<RASModel<BasicMomentumTransportModel>>(type,
																					 alpha,
																					 rho,
																					 U,
																					 alphaRhoPhi,
																					 phi,
																					 transport),
								fluxType_(
									this->coeffDict_.lookup("fluxType")),
								Cmu_(
									IOobject(
										"Cmu",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::NO_WRITE),
									this->mesh_,
									dimensionedScalar("Cmu", dimless, 0.09)),
								Cs_(
									dimensioned<scalar>::lookupOrAddToDict(
										"Cs",
										this->coeffDict_,
										1.35)),
								Cw_(
									dimensioned<scalar>::lookupOrAddToDict(
										"Cw",
										this->coeffDict_,
										1.8)),
								Cc_(
									dimensioned<scalar>::lookupOrAddToDict(
										"Cc",
										this->coeffDict_,
										0.22)),
								Cc1_(
									dimensioned<scalar>::lookupOrAddToDict(
										"Cc1",
										this->coeffDict_,
										1.8)),
								C1_(
									dimensioned<scalar>::lookupOrAddToDict(
										"C1",
										this->coeffDict_,
										1.35)),
								C2_(
									dimensioned<scalar>::lookupOrAddToDict(
										"C2",
										this->coeffDict_,
										1.8)),
								sigmaEps_(
									dimensioned<scalar>::lookupOrAddToDict(
										"sigmaEps",
										this->coeffDict_,
										1.3)),
								sigmak_(
									dimensioned<scalar>::lookupOrAddToDict(
										"sigmak",
										this->coeffDict_,
										1.0)),
								k_(
									IOobject(
										"k",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::MUST_READ,
										IOobject::AUTO_WRITE),
									this->mesh_),
								epsilonTilda_(
									IOobject(
										"epsilon",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::MUST_READ,
										IOobject::AUTO_WRITE),
									this->mesh_),
								R_(
									IOobject(
										"R",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::MUST_READ,
										IOobject::AUTO_WRITE),
									this->mesh_),
								Nij_(
									IOobject(
										"Nij",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::NO_WRITE),
									R_),
								aij_(
									IOobject(
										"aij",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::NO_WRITE),
									this->mesh_,
									dimensionedSymmTensor("aij", dimless, symmTensor::zero)),
								Bij_(
									IOobject(
										"Bij",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::NO_WRITE),
									this->mesh_,
									dimensionedTensor("Bij", dimless, tensor::zero) // TODO:check
									),
								Dt_(
									IOobject(
										"Dt",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::AUTO_WRITE),
									this->mesh_,
									dimensionedTensor("Dt", dimViscosity, tensor::zero)),
								tau_(
									IOobject(
										"tau",
										this->runTime_.timeName(),
										this->mesh_,
										IOobject::NO_READ,
										IOobject::AUTO_WRITE),
									k_ / epsilonTilda_)
		{
			bound(k_, this->kMin_);
			bound(epsilonTilda_, this->epsilonMin_);

			if (type == typeName)
			{
				this->printCoeffs(type);
			}
		}
		// * * * * * * * * * * * * * Member Functions * * * * * * * * * * //
		template <class BasicMomentumTransportModel>
		tmp<volSymmTensorField> ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::sigma() const
		{
			return tmp<volSymmTensorField>(
				new volSymmTensorField(
					IOobject(
						"R",
						this->runTime_.timeName(),
						this->mesh_,
						IOobject::NO_READ,
						IOobject::NO_WRITE),
					R_,
					R_.boundaryField().types()));
		}
		template <class BasicMomentumTransportModel>
		tmp<volSymmTensorField> ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::devSigma() const // devReff()
		{

			return tmp<volSymmTensorField>(
				new volSymmTensorField(
					IOobject(
						"devRhoReff",
						this->runTime_.timeName(),
						this->mesh_,
						IOobject::NO_READ,
						IOobject::NO_WRITE),
					-this->nuEff() * dev(twoSymm(fvc::grad(this->U_)))));
		}
		template <class BasicMomentumTransportModel>
		tmp<fvVectorMatrix> ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::divDevSigma(volVectorField &U) const
		{
			Info << "Using divDevSigma from ChienKEpsilonEARSM_Flux" << endl;
			return (
				fvc::laplacian(
					this->alpha_ * this->rho_ * this->nut(),
					U,
					"laplacian(nuEff,U)") +
				fvc::div(this->alpha_ * this->rho_ * Nij_) - fvc::div(this->alpha_ * this->rho_ * this->nu() * dev2(T(fvc::grad(U)))) - fvm::laplacian(this->alpha_ * this->rho_ * this->nuEff(), U));
		}

		template <class BasicMomentumTransportModel>
		bool ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::read()
		{
			if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
			{
				dimensionedScalar Cs_;
				dimensionedScalar Cw_;
				dimensionedScalar Cc_;
				dimensionedScalar Cc1_;
				dimensionedScalar C1_;
				dimensionedScalar C2_;
				dimensionedScalar sigmaEps_;
				dimensionedScalar sigmak_;
				Cs_.readIfPresent(this->coeffDict());
				Cw_.readIfPresent(this->coeffDict());
				Cc_.readIfPresent(this->coeffDict());
				Cc1_.readIfPresent(this->coeffDict());
				C1_.readIfPresent(this->coeffDict());
				C2_.readIfPresent(this->coeffDict());
				sigmaEps_.readIfPresent(this->coeffDict());
				sigmak_.readIfPresent(this->coeffDict());
				return true;
			}
			else
			{
				return false;
			}
		}
		template <class BasicMomentumTransportModel>
		void ChienKEpsilonEARSM_Flux<BasicMomentumTransportModel>::correct()
		{
			// RASModel::correct();
			if (!this->turbulence_)
			{
				return;
			}
			const alphaField &alpha = this->alpha_;
			const rhoField &rho = this->rho_;
			const volVectorField &U = this->U_;
			const surfaceScalarField &phi = this->phi_;
			const volScalarField &nu = this->nu();
			volTensorField gradU = dev(fvc::grad(U)); // fvc::grad(U) //TODO:check

			tau_ = max(k_ / epsilonTilda_, 6.0 * sqrt(nu / epsilonTilda_));

			//**** calcEARSM ***//
			volSymmTensorField S = tau_ * symm(gradU);
			volTensorField W = tau_ * skew(gradU);

			volTensorField WW = W & W;
			volTensorField SS = S & S;
			volTensorField SSS = SS & S;
			volTensorField SW = S & W;
			volTensorField WS = W & S;
			volTensorField SWW = SW & W;
			volTensorField WWS = WW & S;
			volTensorField WSWW = (WS & W) & W;
			volTensorField WWSW = (WW & S) & W;

			volScalarField IIs = tr(SS);
			volScalarField IIw = tr(WW);
			volScalarField IIIs = tr(SSS);

			volScalarField P1 = (1.8 * (sqr(1.8) / 27.0 + (9.0 / 20.0) * IIs - (2.0 / 3.0) * IIw));
			volScalarField P2 = sqr(P1) - pow((sqr(1.8) / 9.0 + (9.0 / 10.0) * IIs + (2.0 / 3.0) * IIw), 3);

			// N is initialized first and then recalculated later.
			// volScalarField N = 1.8+9/4*pow(2*0.09*IIs,0.5);
			volScalarField N(
				IOobject(
					"N",
					this->runTime_.timeName(),
					this->mesh_,
					IOobject::NO_READ,
					IOobject::NO_WRITE),
				this->mesh_,
				dimensionedScalar("N", dimless, 0),
				zeroGradientFvPatchScalarField::typeName);

			forAll(P2, cellI)
			{
				if (P2[cellI] < 0.0)
				{
					N[cellI] = 1.8 / 3.0 + 2.0 * pow(pow(P1[cellI], 2.0) - P2[cellI], 1.0 / 6.0) *
											   cos(1.0 / 3.0 * acos(P1[cellI] / (pow(pow(P1[cellI], 2.0) - P2[cellI], 0.5))));
				}
				else
				{
					N[cellI] = 1.8 / 3.0 + pow(P1[cellI] + pow(P2[cellI], 0.5), 1.0 / 3.0) + sign(P1[cellI] - pow(P2[cellI], 0.5)) * pow(mag(P1[cellI] - pow(P2[cellI], 0.5)), 1.0 / 3.0);
				}
			}
			N.correctBoundaryConditions();

			volScalarField Q = 5.0 / 6 * (sqr(N) - 2 * IIw) * (2 * sqr(N) - IIw);
			volScalarField beta_1 = -(N * (2 * sqr(N) - 7 * IIw) / Q);
			volScalarField beta_3 = -12 * (tr(SWW) / (Q * N));
			volScalarField beta_4 = -2 * ((sqr(N) - 2 * IIw) / Q);
			volScalarField beta_6 = -6 * (N / Q);
			volScalarField beta_9 = 6.0 / Q;

			volScalarField y = wallDist(this->mesh_).y();
			volScalarField Re_t = pow(k_, 0.5) * (y / nu);
			volScalarField yStar = 2.4 * pow(Re_t, 0.5) + 0.003 * pow(Re_t, 2.0);
			volScalarField f1 = 1.0 - exp(-yStar / 26.0);

			volScalarField beta_2_lowre = (3.0 * 1.8 - 4.0) / (max(IIs, 5.74)) * (1 - pow(f1, 2));
			volScalarField beta_4_lowre = pow(f1, 2) * beta_4 - 1.8 / (2 * max(IIs, 5.74)) * (1 - pow(f1, 2));

			aij_ = symm(
				f1 * beta_1 * S + beta_2_lowre * (SS - (1.0 / 3.0) * IIs * I) + pow(f1, 2) * beta_3 * (WW - (1.0 / 3.0) * IIw * I) + beta_4_lowre * (SW - WS) + f1 * beta_6 * (SWW + WWS - 2.0 / 3 * tr(SWW) * I) + pow(f1, 2) * beta_9 * (WSWW - WWSW));

			Cmu_ = -f1 * (beta_1 + IIw * beta_6) / 2.0;

			Nij_ = k_ * (aij_ - 2.0 * Cmu_ * S); // Non-linear part of Reynolds Stresses
			R_ = Nij_ + (2.0 / 3) * k_ * I;

			if (true)
			{
				Info << "aij_ TRACE  "
					 << " max:" << max(tr(aij_)).value()
					 << " min:" << min(tr(aij_)).value()
					 << " avg:" << average(tr(aij_)).value()
					 << endl;

				Info << "Nij_ TRACE  "
					 << " max:" << max(tr(Nij_)).value()
					 << " min:" << min(tr(Nij_)).value()
					 << " avg:" << average(tr(Nij_)).value()
					 << endl;

				Info << "R_ TRACE  "
					 << " max:" << max(tr(R_)).value()
					 << " min:" << min(tr(R_)).value()
					 << " avg:" << average(tr(R_)).value()
					 << endl;
			}

			//***** Calc KEpslion ****//
			volScalarField G(this->GName(), -R_ && T(gradU));

			// calculate yPlus
			volScalarField yPlus(
				IOobject(
					"yPlus",
					this->runTime_.timeName(),
					this->mesh_,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE),
				this->mesh_,
				dimensionedScalar("yPlus", dimless, 0));

			const fvPatchList &patches = this->mesh_.boundary();
			forAll(patches, patchi)
			{
				const fvPatch &patch = patches[patchi];

				if (isA<wallFvPatch>(patch))
				{
					yPlus.boundaryFieldRef()[patchi] = y.boundaryField()[patchi] * sqrt(nu.boundaryField()[patchi] * mag(U.boundaryField()[patchi].snGrad())) / nu.boundaryField()[patchi];
				}
			}

			const volScalarField E(-2.0 * nu * epsilonTilda_ / pow(y, 2) * exp(-0.5 * yStar) * exp(-0.04 * yPlus)); // TODO:CHECK
			const volScalarField D(2.0 * nu * k_ / pow(y, 2) * exp(-0.04 * yPlus));									// TODO:CHECK

			// Dissipation rate equation
			tmp<fvScalarMatrix> epsEqn(
				fvm::ddt(epsilonTilda_) + fvm::div(phi, epsilonTilda_) - fvm::laplacian(DepsilonEff(), epsilonTilda_) ==
				C1_ * G * epsilonTilda_ / k_ - fvm::Sp(C2_ * f2() * epsilonTilda_ / k_, epsilonTilda_) + E // TODO:check
			);
			epsEqn.ref().relax();
			solve(epsEqn);
			bound(epsilonTilda_, this->epsilonMin_);

			// Turbulent kinetic energy equation
			tmp<fvScalarMatrix> kEqn(
				fvm::ddt(k_) + fvm::div(phi, k_) - fvm::laplacian(DkEff(), k_) ==
				G - fvm::Sp((epsilonTilda_ + D) / k_, k_) // TODO:check
			);
			kEqn.ref().relax();
			solve(kEqn);
			bound(k_, this->kMin_);

			// Re-calculate turbulent viscosity
			correctNut();

			//****flux Calculation ****//
			volScalarField Q1 = sqr(Cs_) * IIs + sqr(Cw_) * IIw;
			volScalarField Q2 = 2.0 / 3 * pow(Cs_, 3) * IIIs + 2 * Cs_ * sqr(Cw_) * tr(SWW);

			volScalarField Gc = 0.5 * (2 * Cc1_ - 2.19 + G / epsilonTilda_);
			volTensorField CSCW = Cs_ * S + Cw_ * W;
			// Bij_ = (((sqr(Gc) - Q1/2.0)*I - Gc*CSCW*tau_+ (CSCW & CSCW)*sqr(tau_))/(Gc*(sqr(Gc) - Q1/2.0) + Q2/2.0);
			// Info << "test" << endl;
			Bij_ = ((sqr(Gc) - Q1 / 2.0) * I - Gc * CSCW + (CSCW & CSCW)) / (Gc * (sqr(Gc) - Q1 / 2.0) + Q2 / 2.0); // TODO:

			// Re-calculate turbulent transport diffusivity
			correctDt();
		}

		// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	} // End namespace    RASModels_
	//} // End namespace incompressible
} // End namespace Foam
// *************************************************************** //
