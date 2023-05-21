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

Class
    Foam::RASModels::bslKOmega_Flux

Group
    grpRASTurbulence

See also
    Foam::kOmegaSST

SourceFiles
    bslKOmega_Flux.C

\*---------------------------------------------------------------------------*/

#ifndef bslKOmega_Flux_H
#define bslKOmega_Flux_H

#include "kOmegaSSTBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class bslKOmega_Flux Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class bslKOmega_Flux
:
    public Foam::kOmegaSST
    <
        eddyViscosity<RASModel<BasicMomentumTransportModel>>,
        BasicMomentumTransportModel
    >
{

protected:
	
	word  fluxType_;
        volScalarField Cmu_;
        dimensionedScalar Cs_;
        dimensionedScalar Cw_;
        dimensionedScalar Cc_;
        dimensionedScalar Cc1_;
        
        dimensionedScalar alphaOmega_;
	dimensionedScalar beta_;
        dimensionedScalar gamma_;	
        dimensionedScalar alphaK1_;
        dimensionedScalar alphaOmega2_;
        dimensionedScalar gamma1_;
        dimensionedScalar gamma2_;
        dimensionedScalar Prt_;
        dimensionedScalar kappa_;
        dimensionedScalar C1_;
        dimensionedScalar A1_;

     	volSymmTensorField R_;
        volSymmTensorField Nij_;
        volSymmTensorField aij_;
        
        volTensorField Bij_;        
        //turbulence diffusion coefficient
        volTensorField Dt_;
            
        //turbulent timescale
        volScalarField tau_;
            
        virtual void correctNut();
        virtual void correctDt();
        
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("bslKOmega_Flux");


    // Constructors

        //- Construct from components
        bslKOmega_Flux
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~bslKOmega_Flux()
    {}
    
       //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                alphaK1_*this->nut_ + this->nu()
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return volScalarField::New
            (
                "DomegaEff",
                alphaOmega_*this->nut_ + this->nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return this->k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return this->omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return volScalarField::New
            (
                "epsilon",
                this->betaStar_*this->k_*this->omega_,
                this->omega_.boundaryField().types()
            );
        }
    
    
    //- Return the source term for the momentum equation
    virtual tmp<volSymmTensorField> sigma() const;

    //- Return the source term for the momentum equation
    virtual tmp<volSymmTensorField> devSigma() const;

    //- Return the source term for the momentum equation
    virtual tmp<fvVectorMatrix> divDevSigma(volVectorField& U) const;

    //- Solve the turbulence equations and correct the turbulence viscosity
    virtual void correct();
    
    //- Return the turbulent diffusion 
    tmp<volTensorField> Dt() const {	
    		return Dt_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "bslKOmega_Flux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
