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
    Foam::RASModels::bslKOmegaEARSM_Flux

Group
    grpRASTurbulence

See also
    Foam::kOmegaSST

SourceFiles
    bslKOmegaEARSM_Flux.C

\*---------------------------------------------------------------------------*/

#ifndef bslKOmegaEARSM_Flux_H
#define bslKOmegaEARSM_Flux_H

#include "kOmegaSSTBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class bslKOmegaEARSM_Flux Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class bslKOmegaEARSM_Flux
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
public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("bslKOmegaEARSM_Flux");


    // Constructors

        //- Construct from components
        bslKOmegaEARSM_Flux
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
    virtual ~bslKOmegaEARSM_Flux()
    {}
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
    #include "bslKOmegaEARSM_Flux.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
