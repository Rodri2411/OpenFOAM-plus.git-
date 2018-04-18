/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::assembleDecomp

Description
    Domain decomposition which conserves mapped patches locally in
    processors. This is necessary in cases where the coupled option
    for energy is used in order to make internal all faces
    of the mapped patch.

    This method uses a run-time selectable method for the designated main
    region. This decomposition will write the decomposition list for each
    of the remaining region which need to use 'manual' decomposition.

    \verbatim
        numberOfSubdomains  4;
        method              assembleDecomp;
        coeffs
        {
            numberOfSubdomains  4;
            method              scotch;
            dataFile            "cellDecompose.dat";
        }
    \endverbatim

SourceFiles
    assembleDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef assembleDecomp_H
#define assembleDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class assembleDecomp Declaration
\*---------------------------------------------------------------------------*/

class assembleDecomp
:
    public decompositionMethod
{
    // Private data

        //- Method dictionary
        dictionary methodDict_;

        //- Method
        autoPtr<decompositionMethod> method_;

        //- Data file for regions
        fileName dataFile_;



    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const assembleDecomp&) = delete;
        assembleDecomp(const assembleDecomp&) = delete;


public:

      //- Runtime type information
        TypeName("assembleDecomp");

    // Constructors

        //- Construct given the decomposition dictionary
        assembleDecomp(const dictionary& decompDict);

    //- Destructor
    virtual ~assembleDecomp()
    {}


    // Member Functions

        //- This decompose does not care about proc boundaries - is all
        //  up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number.
        //  Uses the mesh connectivity (if needed).
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights
        );

        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres,
            const scalarField& cellWeights
        )
        {
            NotImplemented;
            return labelList(0);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
