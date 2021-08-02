/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 Alex Rattner
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
    Foam::thermalPhaseChangeModel::constantPhaseChange

Description
    Dummy model for no phase change. Just set none as the model type.

SourceFiles
    constantPhaseChange.C

\*--------------------------------------------------------------------*/

#ifndef constantPhaseChange_H
#define constantPhaseChange_H

#include "thermalPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{

/*--------------------------------------------------------------------*\
                         Class constantPhaseChange
\*--------------------------------------------------------------------*/

class constantPhaseChange
:
    public thermalPhaseChangeModel
{

// Private Member Functions

	    // in kg/(m2 s)
    	dimensionedScalar mEvapFlux_;

    //- This phase change model needs access to the mesh:
    const fvMesh& mesh_;

    //- Phase change heat
    volScalarField Q_pc_;

        volScalarField calcGradAlphal() const;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

    //- construct from components
    constantPhaseChange
    (
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
    );



    //- Destructor
    ~constantPhaseChange()
    {}


    // Member Functions

        //- Correct the thermalPhaseChange model
        void correct()
        {

        }
        
        //- Return phase change heat
        tmp<volScalarField> Q_pc() const
        {
            return Q_pc_;
        }


        //- Read the transportProperties dictionary and update
        bool read(const dictionary& thermalPhaseChangeProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
