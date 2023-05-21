/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "bslKOmegaEARSM_Flux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// Anmerkung: I == delta_ij = Einheitsmatrix == Kronecker-Delta

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/* = spezielle Initialisierungs-Fkt., die aufgerufen wird um ein neues Objekt der Klasse zu erstellen; 
	kann auch zur Initialisierung von Attributen eines Objkts genutzt werden
*/

//static const dimensionedScalar smallOmega_ ("smallOmega",dimensionSet (0,0,-1,0,0,0,0), 1e-9); //Skalar mit Einheit 1/s
//static const dimensionedScalar smallDimless_ ("smallDimless",dimless, 1e-9); //dimensionsloses Skalar
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::correctNut()
{
    //this->nut_ = this->k_/(this->omega_+smallOmega_);
    this->nut_ = this->k_/(this->omega_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}
template<class BasicMomentumTransportModel>
void bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::correctDt()
{
      if (fluxType_ == "DalyHarlow")
      { 
      	Dt_ == -Cc_* (tensor::I & R_) *tau_;
      }
      else if ( fluxType_ == "SugaAbe")
      {
      	Dt_ == -Cc_*(R_ & R_)/this->k_*tau_;
      }
      else if (fluxType_ == "WWJ")
      {
     	Dt_ == -(1.0 - Cc_)*(Bij_ & (aij_+ 2.0/3*I)*this->k_)*tau_;//TODO:check
      }
      else
      {
      	 WarningInFunction 
      	 <<"fluxType are DalyHarlow, SugaAbe, WWJ"
      	 << endl;	
      }
      
      Info << "Dt max/min/ave: " 
      	   << max(mag(Dt_)).value() << ","
      	   << min(mag(Dt_)).value() << ","
      	   << average(mag(Dt_)).value() <<endl;
      
      Dt_.correctBoundaryConditions();
      fvConstraints::New(this->mesh_).constrain(Dt_);
}


//Vorlage, die alle Turbulenzmodelle umfasst
template<class BasicMomentumTransportModel> 
bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::bslKOmegaEARSM_Flux
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& type
)
:
    Foam::kOmegaSST // die folgenden Templates sind vom kOmegaSST-Template unmittelbar abhängig
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
	// erstelle die folgenden Variablen/Methoden aus den obigen Komponenten
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport
    ),
    fluxType_
    (
	  this->coeffDict_.lookup("fluxType")
    ),
    Cmu_
    (
	IOobject
	(
		"Cmu",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
        this->mesh_,
        dimensionedScalar("Cmu",dimless,0.09)
   ),
   Cs_
   (
	dimensioned<scalar>::lookupOrAddToDict
	(
		"Cs",
		this->coeffDict_,
		1.35
	)
   ),
   Cw_
   (
	dimensioned<scalar>::lookupOrAddToDict
	(
		"Cw",
		this->coeffDict_,
		1.8
	)
   ),    
   Cc_
   ( 
	dimensioned<scalar>::lookupOrAddToDict
	(
		"Cc",
		this->coeffDict_,
		0.22
	)
   ),
   Cc1_
   (
	dimensioned<scalar>::lookupOrAddToDict
	(
		"Cc1",
		this->coeffDict_,
		1.8
	)
   ),    
    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.4403
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            this->coeffDict_,
            1.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.8
        )
    ),

    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            this->coeffDict_,
            1.245
        )
    ),
    R_
    (
	IOobject
	(
		"R",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
        this->mesh_
    ),
    Nij_
    (
	IOobject
	(
		"Nij",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
        R_
    ),
    aij_
    (
	IOobject
	(
		"aij",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
        this->mesh_,
        dimensionedSymmTensor ("aij",dimless, symmTensor::zero)
     ),
     Bij_
     (
	IOobject
	(
		"Bij",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
        this->mesh_,
        dimensionedTensor("Bij",dimless,tensor::zero) 
    ),
    Dt_
    (
	IOobject
	(
		"Dt",
		this->runTime_.timeName(),
		this->mesh_,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
        this->mesh_,
        dimensionedTensor("Dt",dimViscosity,tensor::zero)
    ),
    tau_
    (
    	//6.*sqrt( this->nu() / (this->betaStar_*this->k_*this->omega_))
       max((1./(this->betaStar_*this->omega_)),6.*sqrt( this->nu() / (this->betaStar_*this->k_*this->omega_)) ) 
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
        Info << "BSL-EARSM Turbulence Model." << endl;
    }
}

// Reynoldsspannungen in symmetrischem Tensor, oberes rechtes Dreieck besetzt, 6 statt 9 Einträge
template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::sigma() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*this->k_ - (this->nut_)*dev(twoSymm(fvc::grad(this->U_)))+Nij_/this->rho_,
            this->k_.boundaryField().types()
        )
    );
}

// Methode, die durch Template "CViscousStress" bereitgestellt wird
template<class BasicMomentumTransportModel>
tmp<volSymmTensorField> bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::devSigma() const
{
    Info << "Using devSigma from BSL-EARSM" << endl;
	
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
        //    this->alpha_*this->rho_*R_
       // - (this->alpha_*this->rho_*this->nu())
       //*dev(twoSymm(fvc::grad(this->U_)))

           dev(Nij_) - this->rho_*this->nuEff()*dev(twoSymm(fvc::grad(this->U_))) //TODO:check
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvVectorMatrix> bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::divDevSigma(volVectorField& U) const
{
    Info << "Using divDevSigma from EARSM BSL KOmega_SST" << endl;
      
    
    return
    (
           fvc::laplacian
            (
                this->alpha_*this->rho_*this->nut(),
                U,
                "laplacian(nuEff,U)"
            )
          + fvc::div(this->alpha_*this->rho_*Nij_)
          - fvc::div(this->alpha_*this->rho_*this->nu()*dev2(T(fvc::grad(U))))
          - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );
}

template<class BasicMomentumTransportModel>
void bslKOmegaEARSM_Flux<BasicMomentumTransportModel>::correct()
{

	if (!this->turbulence_)
    	{
        	return;
    	}

    	// Local references
    	const alphaField& alpha = this->alpha_;
    	const rhoField& rho = this->rho_;
    	const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    	const volVectorField& U = this->U_;
    	volScalarField& nut = this->nut_;
    	const volScalarField& nu = this->nu();
    
    	const Foam::fvConstraints& fvConstraints
    	(
        	Foam::fvConstraints::New(this->mesh_)
    	);

    	BasicMomentumTransportModel::correct();

	// skalares Feld mit Werten im Zellmittelpunkt
    	volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

     	volTensorField gradU = dev(fvc::grad(U));
	
    // Zeitskala mit Kolmogorov-Limiter (Gl. 5) 
	tau_ = max((1./(this->betaStar_*this->omega_)),6.*sqrt( nu /(this->betaStar_*this->k_*this->omega_)));
	// 

	volSymmTensorField S = tau_*symm(gradU);
	volTensorField W = tau_*skew(gradU);
	
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
        
    	volScalarField P1 =(1.8*(sqr(1.8)/27.0 + (9.0/20.0)*IIs - (2.0/3.0)*IIw));
	volScalarField P2 = sqr(P1) - pow((sqr(1.8)/9.0 + (9.0/10.0)*IIs + (2.0/3.0)*IIw), 3);

	volScalarField N (
    		IOobject
        	(
            		"N",
            		this->runTime_.timeName(),
            		this->mesh_,
            		IOobject::NO_READ,
            		IOobject::NO_WRITE 
        	),
        	this->mesh_,
        	dimensionedScalar ("N",dimless,0),
        	zeroGradientFvPatchScalarField::typeName
    		);
	
	forAll(P2, cellI)
	{
		if (P2[cellI] < 0.0)
		{
			N[cellI] = 1.8/3.0 + 2.0*pow(pow(P1[cellI],2.0)-P2[cellI],1.0/6.0)*
			cos(1.0/3.0*acos(P1[cellI]/(pow(pow(P1[cellI],2.0)-P2[cellI],0.5))));
		}
		else
		{
			N[cellI] = 1.8/3.0 + pow(P1[cellI]+pow(P2[cellI],0.5),1.0/3.0) 
			+sign(P1[cellI]-pow(P2[cellI],0.5))
			*pow(mag(P1[cellI]- pow(P2[cellI],0.5)), 1.0/3.0);
		}
	}
	N.correctBoundaryConditions();


//TODO: check
	volScalarField Q0 = (sqr(N)-2.*IIw)/A1_;
    	volScalarField Q = Q0*(2.*sqr(N)-IIw)/6.;
    	
    	volScalarField beta1 = -1.*N/Q0;
    	volScalarField beta2 = 0.*N;     //0 work around
    	volScalarField beta3 = -2.*tr(SWW)/(N*Q);
    	volScalarField beta4 = -1./Q0;
    	volScalarField beta6 = -1.*N/Q;
    	volScalarField beta9 = 1./Q;
    	
    	/*volScalarField beta1 = -1.*N/(Q0+smallDimless_);
    	volScalarField beta2 = 0.*N;     //0 work around
    	volScalarField beta3 = -2.*tr(SWW)/(N*Q+smallDimless_);
    	volScalarField beta4 = -1./(Q0+smallDimless_);
    	volScalarField beta6 = -1.*N/(Q+smallDimless_);
    	volScalarField beta9 = 1./(Q+smallDimless_);*/

    	/*volScalarField Q = 5.0/6*(sqr(N)-2*IIw)*(2*sqr(N)-IIw);
	volScalarField beta_1=-(N*(2*sqr(N)-7*IIw)/Q); 
	volScalarField beta_3=-12*(tr(SWW)/(Q*N));
	volScalarField beta_4=-2*((sqr(N)-2*IIw)/Q);  
	volScalarField beta_6=-6*(N/Q);
	volScalarField beta_9= 6.0/Q;*/


	// Objekt des Anisotropietensor wird definiert und durch (Gl. 2) initialisiert
    	aij_ = symm(
            beta1*S
            //+ beta2*(SS-(1./3.)*IIs_*I)
            + beta3*(WW-(1./3.)*IIw*I)
            + beta4*(SW-WS)
            + beta6*(SWW + WWS - (2./3.)*tr(SWW)*I - IIw*S)
            //+ beta9*(WSWW - WWSW + 0.5*IIw_*(SW - WS))
            );
    
        Nij_ =  this->k_ * aij_; // Non-linear part of Reynolds Stresses  
	R_ = Nij_ + (2.0/3) * this->k_ * I;
	
	
	if (true) {                     
	Info<< "aij_ TRACE  "
       <<" max:"<< max(tr(aij_)).value() 
       <<" min:"<< min(tr(aij_)).value()
       <<" avg:"<< average(tr(aij_)).value()
       << endl;
             
       Info<< "Nij_ TRACE  "
       <<" max:"<< max(tr(Nij_)).value() 
       <<" min:"<< min(tr(Nij_)).value()
       <<" avg:"<< average(tr(Nij_)).value()
       << endl;
    
       Info<< "R_ TRACE  "
       <<" max:"<< max(tr(R_)).value() 
       <<" min:"<< min(tr(R_)).value()
       <<" avg:"<< average(tr(R_)).value()
       << endl;
       }
	
      volScalarField G(this->GName(), -R_ && T(gradU));	

    // skalares Feld mit Werten im Zellmittelpunkt 
    // GName = Hilfs-Fkt. zur Wiedergabe des Namens des turbulenten G-Feldes
    //volScalarField::Internal G(this->GName(), -(tauij_&&gradU_));
    //tgradU.clear();

    // Update omega and G at the wall
    this->omega_.boundaryFieldRef().updateCoeffs();

    // alphaOmega2 = Sigma_(Omega 2) in Gl. 
    volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation (Gl. 13_2)
        tmp<fvScalarMatrix> omegaEqn
        (
		// zeitliche Ableitung = D omega / Dt 			--> Zeitableitungs-Term
            fvm::ddt(alpha, rho, this->omega_)
		// Divergenz						--> ????????????????
          + fvm::div(alphaRhoPhi, this->omega_)
		// Laplace, DomegaEff = effektive Diffusivität von omega --> Laplace-Term
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), this->omega_) //TODO: omega
         ==
		// = (gamma * omaga / k) * P_k 				--> P_k-Term
			// mit P_k = min(G, ...) (Gl. 14) 
            alpha()*rho()*gamma
           *min
            (
                G,
                this->c1_*this->betaStar_*this->k_()*this->omega_() //TODO: omega
            )*this->omega_()/this->k_()
		// = Quellterm = 2/3 * gamma * divU * omega		--> ????????????????
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, this->omega_)
		// = Quellterm = - beta * omega²			--> beta - omega²-Term
          - fvm::Sp(alpha()*rho()*beta*this->omega_(), this->omega_)
		// = Quellterm = - (F1 - 1) * CDkOmega / omega * omega	--> CDkOmega-Term
          - fvm::SuSp //TODO:check
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/this->omega_(),
                this->omega_
            )
        );

        omegaEqn.ref().relax();
        fvConstraints.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(this->omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvConstraints.constrain(this->omega_);
        bound(this->omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation (Gl. 13_1)
    tmp<fvScalarMatrix> kEqn
    (
	// zeitliche Ableitung = Dk / Dt				--> Zeitableitungs-Term
        fvm::ddt(alpha, rho, this->k_)
	// Divergenz							--> ??????????????
      + fvm::div(alphaRhoPhi, this->k_)
	// Laplace, DkEff = effektive Diffusivität von k 		--> Laplace-Term
      - fvm::laplacian(alpha*rho*this->DkEff(F1), this->k_) //TODO:check
     ==
	// = P_k = min(G, ...) 						--> P_k-Term
        min
	(
		alpha()*rho()*G, 
		(this->c1_*this->betaStar_)*alpha()*rho()*this->k_()*this->omega_() //TODO:check
	)
	// = Quellterm = - 2/3 * divU * k				--> ???????????????
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, this->k_)
	// = Quellterm = - betaStar * omega * k				--> beta*-k-Omega-Term
      - fvm::Sp(alpha()*rho()*this->betaStar_*this->omega_(),this->k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(this->k_);
    bound(this->k_, this->kMin_);
    correctNut();
    
    //****flux Calculation ****//
	volScalarField Q1 = sqr(Cs_)*IIs + sqr(Cw_)*IIw;
	volScalarField Q2 = 2.0/3*pow(Cs_,3)*IIIs + 2*Cs_*sqr(Cw_)*tr(SWW);
	
	volScalarField Gc = 0.5*(2*Cc1_ - 2.19 + G/(0.09*this->k_*this->omega_) ); //TODO:check
	volTensorField CSCW = Cs_*S + Cw_*W;
	//Bij_ = (((sqr(Gc) - Q1/2.0)*I - Gc*CSCW*tau_+ (CSCW & CSCW)*sqr(tau_))/(Gc*(sqr(Gc) - Q1/2.0) + Q2/2.0);
	//Info << "test" << endl;
	Bij_ = ((sqr(Gc) - Q1/2.0)*I - Gc*CSCW + (CSCW & CSCW))/(Gc*(sqr(Gc) - Q1/2.0) + Q2/2.0);//TODO:
	
  // Re-calculate turbulent transport diffusivity
	correctDt();
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

//
