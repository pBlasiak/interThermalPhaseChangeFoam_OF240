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

\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "constantPhaseChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace thermalPhaseChangeModels
{
    defineTypeNameAndDebug(constantPhaseChange, 0);
    addToRunTimeSelectionTable
    (
        thermalPhaseChangeModel, 
        constantPhaseChange, 
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeModels::constantPhaseChange::constantPhaseChange
(
        const word& name,
        const dictionary& thermalPhaseChangeProperties,
        const twoPhaseThermalMixture& twoPhaseProperties,
        const volScalarField& T,
        const volScalarField& alpha1
)
:
    thermalPhaseChangeModel
    (
        name, 
        thermalPhaseChangeProperties, 
        twoPhaseProperties, 
        T, 
        alpha1
    ),
    mEvapFlux_(thermalPhaseChangeProperties.lookup("evapMassFlux")),
    mesh_(T.mesh()),
    Q_pc_
    (
        IOobject
        (
            "PhaseChangeHeat",
            T_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar( "dummy", dimensionSet(1,-1,-3,0,0,0,0), 0 )
    )
{
	Q_pc_ = mEvapFlux_*h_lv_*calcGradAlphal();
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam ::volScalarField Foam::thermalPhaseChangeModels::constantPhaseChange::calcGradAlphal() const
{
	volScalarField limitedAlpha1 = min(max(alpha1_, scalar(0)), scalar(1));
	return mag(fvc::grad(limitedAlpha1));
}


bool Foam::thermalPhaseChangeModels::constantPhaseChange::
read(const dictionary& thermalPhaseChangeProperties)
{
    thermalPhaseChangeModel::read(thermalPhaseChangeProperties);
        thermalPhaseChangeProperties.lookup("evapMassFlux") >> mEvapFlux_;
		Q_pc_ = mEvapFlux_*h_lv_*calcGradAlphal();

    return true;
}


// ************************************************************************* //
