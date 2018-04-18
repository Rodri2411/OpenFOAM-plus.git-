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
    Foam::fvMatrixAssemble

Description


SourceFiles
    fvMatrixAssemble.C


\*---------------------------------------------------------------------------*/

#ifndef fvMatrixAssemble_H
#define fvMatrixAssemble_H

#include "fvMatrix.H"
#include "fvMatricesFwd.H"
#include "fvSolution.H"
#include "lduPrimitiveMeshAssemble.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mappedPatchBase;
class cyclicLduInterface;

/*---------------------------------------------------------------------------*\
                           Class fvMatrixAssemble Declaration
\*---------------------------------------------------------------------------*/

class fvMatrixAssemble
:
    public lduMatrix
{
private:

    // Private Data

        //- Reference to lduPrimitiveMeshAssemble
        const lduPrimitiveMeshAssemble& primitiveMesh_;

        //- Pointer list to fvMatrix
        PtrList<fvMatrix<scalar> > matrices_;

        //- Pointer list to fields
        UPtrList<volScalarField> psis_;

        //- Number of fvMatrices
        label nMatrix_;

        //- Name of the single field solved.
        word psiName_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Source term
        Field<scalar> source_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for internal cells
        FieldField<Field, scalar> internalCoeffs_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //  for boundary cells
        FieldField<Field, scalar> boundaryCoeffs_;

        //- List of interfaces for the single field
        lduInterfaceFieldPtrsList interfaces_;


        //- Face areas of the assembled mesh
        vectorField* faceAreasPtr_;

        //- Cell volumes the assembled mesh
        scalarField* cellVolumesPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMatrixAssemble(const fvMatrixAssemble&);

        //- Disallow default bitwise assignment
        void operator=(const fvMatrixAssemble&);

        //- Update matrix coefficients on the added faces
        //void updateCoeffs(const label matrixI, const label patchI);

        //- Add boundary contribution to diagonal
        void addBoundaryDiag(scalarField& diag) const;

        //- Add boundary sources
        void addBoundarySource(Field<scalar>& source, const bool couples) const;

        //- Transfer lower, upper and source fields and clean individual
        /// components
        void transferFieldsAndClean();

        //- Update assembled source, internalCoeffs and boundaryCoeffs
        //  once all the fvMatrix have added appended
        void update();

        //- Update to global patchID neigbourds for cyclics, etc polyPatches
        // This is needed for the interface init/updates
        //void updateNnbrPatchID(bool reset);


        //- Calculate face areas and cell volumes for the assembled matrix
        // to be used by the GAMG solver
        void calFaceAreasCellVolumes();


public:

    //- Runtime type information
    TypeName("fvMatrixAssemble");


    // Constructors


         //- Construc from lduPrimitiveMeshAssemble, dimension and psiName
        fvMatrixAssemble
        (
            const lduPrimitiveMeshAssemble&,
            const dimensionSet& ds,
            word psiName
        );


    //- Destructor
    virtual ~fvMatrixAssemble();


    // Member functions

        // Editing

            //- Adding fvMatrix
            void addFvMatrix(const fvMatrix<scalar>&);

            //- Clear
            void clear();


        // Solution

            //- Solve
            SolverPerformance<scalar> solve(const dictionary&);


        // Access

            //- Return the LDU mesh from which the addressing is obtained
            const lduMesh& mesh() const
            {
                return primitiveMesh_;
            }


            //- Return faceAreas
            const vectorField& Sf() const
            {
                return *faceAreasPtr_;
            }


            //- Return cell Volumes
            const scalarField& V() const
            {
                return *cellVolumesPtr_;
            }

            //- Return non-const source
            Field<scalar>& source()
            {
                return source_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
