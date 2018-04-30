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

\*---------------------------------------------------------------------------*/

#include "fvMatrixAssemble.H"
#include "dictionary.H"

#include "GAMGSolver.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMatrixAssemble, 0);
}

//* * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
/*
void Foam::fvMatrixAssemble::updateCoeffs
(
    const label matrixI,
    const label patchI
)
{
    const volScalarField::Boundary& tbf =
        matrices_[matrixI].mesh().thisDb().lookupObject
        <
            volScalarField
        >("T").boundaryField();

    typedef typename
        compressible::turbulentTemperatureRadCoupledMixedFvPatchScalarField
            mappedTemp;

    const polyPatch& pp = psis_[matrixI].mesh().boundaryMesh()[patchI];


    if (isA<mappedTemp>(tbf[patchI]))
    {
        const mappedTemp& tempFvPatch = refCast<const mappedTemp>(tbf[patchI]);

        const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

        tmp<scalarField> vf = tempFvPatch.alphaDeltaVf();

        // Correct source in cell next to owner mappPolyPatch
        // required by jump deltaH
        tmp<scalarField> deltaH = tempFvPatch.deltaH();

        const scalarField sourceCorrection
        (
            tempFvPatch.alphaSfDelta()
           *(
               deltaH()*vf()
             + tempFvPatch.deltaQflux()/tempFvPatch.beta()
            )
        );

        tmp<scalarField> alphaSfDelta = tempFvPatch.alphaSfDelta();

        forAll (tbf[patchI], faceI)
        {
            if (primitiveMesh_.faceBoundMap()[matrixI][patchI].size() > 0)
            {
                label globalFaceI =
                    primitiveMesh_.faceBoundMap()[matrixI][patchI][faceI];

                if (globalFaceI != -1)
                {
                    const scalar corr(vf()[faceI]*alphaSfDelta()[faceI]);

                    if (mpp.owner())
                    {
                        const labelUList& l = lduAddr().lowerAddr();
                        diag()[l[globalFaceI]] += corr;
                        lower()[globalFaceI] -= corr;
                    }
                    else
                    {
                        upper()[globalFaceI] -= corr;
                        const labelUList& u = lduAddr().upperAddr();
                        diag()[u[globalFaceI]] += corr;

                    }
                    // negSumDiag();
                    // I need to negSum the new upper value into the diag
                    // Emulating the laplacian operator

                    //const labelUList& l = primitiveMesh_.lduAddr().lowerAddr();
                    //const labelUList& u = primitiveMesh_.lduAddr().upperAddr();
                    //diag()[l[globalFaceI]] += corr;
                    //diag()[u[globalFaceI]] += corr;

                }
                else
                {
                    FatalErrorInFunction
                        << "Can't find globalFaceI"
                        << exit(FatalError);
                }
            }
        }



        const labelUList& fc =
            psis_[matrixI].boundaryField()[patchI].patch().faceCells();

        forAll(fc, i)
        {
            label localCelli = fc[i];
            label globalCelli = primitiveMesh_.cellOffsets()[matrixI] + localCelli;
            source_[globalCelli] += sourceCorrection[i];
        }

    }
    else
    {
       FatalErrorInFunction
            << "Patch is not a "
            << "compressible::turbulentTemperatureRadCoupledMixedFvPatchScalarField."
            << "It is : "  << tbf[patchI].type()
            << exit(FatalError);
    }
}
*/

void Foam::fvMatrixAssemble::addBoundaryDiag
(
    scalarField& diag
) const
{
    forAll(internalCoeffs_, patchi)
    {
        // Using any of the fvMatrix to addToInternal
        matrices_[0].addToInternalField
        (
            primitiveMesh_.lduAddr().patchAddr(patchi),
            internalCoeffs_[patchi],
            diag
        );
    }
}


void Foam::fvMatrixAssemble::addBoundarySource
(
    Field<scalar>& source,
    const bool couples
) const
{
    forAll(psis_, i)
    {
        const GeometricField<scalar, fvPatchField, volMesh>& psi = psis_[i];

        forAll(psi.boundaryField(), patchi)
        {
            const fvPatchField<scalar>& ptf = psi.boundaryField()[patchi];

            label globalPatchi = primitiveMesh_.patchMap()[i][patchi];

            // In global field mappedWall is not patch
            if (globalPatchi != -1)
            {
                const Field<scalar>& pbc = boundaryCoeffs_[globalPatchi];

                if (!ptf.coupled())
                {
                    matrices_[0].addToInternalField
                    (
                        primitiveMesh_.lduAddr().patchAddr(globalPatchi),
                        pbc,
                        source
                    );
                }
                else if (couples)
                {
                    const tmp<Field<scalar>> tpnf = ptf.patchNeighbourField();
                    const Field<scalar>& pnf = tpnf();

                    const labelUList& addr =
                        primitiveMesh_.lduAddr().patchAddr(globalPatchi);

                    forAll(addr, facei)
                    {
                        source[addr[facei]] +=
                            cmptMultiply(pbc[facei], pnf[facei]);
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMatrixAssemble::fvMatrixAssemble
(
    const lduPrimitiveMeshAssemble& mesh,
    const dimensionSet& ds,
    word psiName
)
:
    lduMatrix(mesh),
    primitiveMesh_(mesh),
    matrices_(0),
    psis_(0),
    nMatrix_(0),
    psiName_(psiName),
    dimensions_(ds),
    source_(primitiveMesh_.lduAddr().size(), 0.0),
    internalCoeffs_(primitiveMesh_.patchAddr().size()),
    boundaryCoeffs_(internalCoeffs_.size()),
    interfaces_(),
    faceAreasPtr_(nullptr),
    cellVolumesPtr_(nullptr)
{
    lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), 0.0);
    upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), 0.0);
    diag().setSize(primitiveMesh_.lduAddr().size(), 0.0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMatrixAssemble::transferFieldsAndClean()
{

    //Reset local mesh
    lower().setSize(primitiveMesh_.lduAddr().upperAddr().size(), 0.0);
    upper().setSize(primitiveMesh_.lduAddr().upperAddr().size(), 0.0);
    diag().setSize(primitiveMesh_.lduAddr().size(), 0.0);
    source_.setSize(primitiveMesh_.lduAddr().size(), 0.0);

    const labelListList& procFaceMap = primitiveMesh_.faceMap();
    const labelList& cellMap = primitiveMesh_.cellOffsets();

    // Move append contents into intermediate list
    for (label i=0; i < nMatrix_; i++)
    {
        scalarField& lowerSub = matrices_[i].lower();
        scalarField& upperSub = matrices_[i].upper();
        scalarField& diagSub = matrices_[i].diag();
        scalarField& sourceSub = matrices_[i].source();

        forAll (lowerSub, facei)
        {
            lower()[procFaceMap[i][facei]] = lowerSub[facei];
            upper()[procFaceMap[i][facei]] = upperSub[facei];
        }

        forAll (diagSub, celli)
        {
            const label globalCelli = cellMap[i] + celli;
            diag()[globalCelli] = diagSub[celli];
            source_[globalCelli] = sourceSub[celli];
        }

        lowerSub.clear();
        upperSub.clear();
        diagSub.clear();
        sourceSub.clear();
    }
}

/*
void Foam::fvMatrixAssemble::updateNnbrPatchID(bool reset)
{
    for (label i=0; i < nMatrix_; i++)
    {
        polyBoundaryMesh& patches =
            const_cast<polyBoundaryMesh&>(psis_[i].mesh().boundaryMesh());

        forAll (patches, patchI)
        {
            if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                cyclicPolyPatch& pp = refCast<cyclicPolyPatch>(patches[patchI]);

                if (!reset)
                {
                    label nbrId = pp.neighbPatchID();

                    pp.neighbPatchID() =
                        primitiveMesh_.patchMap()[i][nbrId];
                }
                else
                {
                    pp.neighbPatchID() = -1;
                }
            }
            else if (isA<cyclicAMIPolyPatch>(patches[patchI]))
            {
                cyclicAMIPolyPatch& pp =
                    refCast<cyclicAMIPolyPatch>(patches[patchI]);

                if (!reset)
                {
                    label nbrId = pp.neighbPatchID();
                    pp.neighbPatchID() =
                        primitiveMesh_.patchMap()[i][nbrId];
                }
                else
                {
                    pp.neighbPatchID() = -1;
                }
            }
        }
    }
}
*/

void Foam::fvMatrixAssemble::calFaceAreasCellVolumes()
{
    const labelListList& procFaceMap = primitiveMesh_.faceMap();

    // Create the storage
    faceAreasPtr_ = new vectorField(lduAddr().upperAddr().size());
    vectorField& faceAreas = *faceAreasPtr_;

    cellVolumesPtr_ = new scalarField(lduAddr().size());
    scalarField& cellVolumes = *cellVolumesPtr_;

    for (label i=0; i < nMatrix_; i++)
    {
        scalarField& lowerSub = matrices_[i].lower();
        const vectorField& areas = psis_[i].mesh().Sf();

        forAll (lowerSub, facei)
        {
            faceAreas[procFaceMap[i][facei]] = areas[facei];
        }

        const polyBoundaryMesh& patches = psis_[i].mesh().boundaryMesh();

        // Fill faceAreas for new faces
        forAll (patches, patchI)
        {
            forAll (patches[patchI], faceI)
            {
                if (primitiveMesh_.faceBoundMap()[i][patchI].size() > 0)
                {
                    const mappedPatchBase& mpp =
                        refCast<const mappedPatchBase>(patches[patchI]);

                    const label globalFaceI =
                        primitiveMesh_.faceBoundMap()[i][patchI][faceI];

                    if (globalFaceI != -1 && mpp.owner())
                    {
                        faceAreas[globalFaceI] =
                            psis_[i].mesh().boundary()[patchI].Sf()[faceI];
                    }
                }
            }
        }

        // Fill cellVolumes
        const scalarField& V = psis_[i].mesh().V();
        const label cellOffset = primitiveMesh_.cellOffsets()[i];

        forAll (psis_[i], localCellI)
        {
            cellVolumes[cellOffset + localCellI] = V[localCellI];
        }
    }
}


void Foam::fvMatrixAssemble::update()
{

    if (interfaces_.empty())
    {
        interfaces_.setSize(internalCoeffs_.size());
        for (label i=0; i < nMatrix_; i++)
        {
            const typename volScalarField::Boundary& bpsi =
                psis_[i].boundaryField();

            forAll (bpsi, patchI)
            {
                if (matrices_[i].mesh().interfaces().set(patchI))
                {
                    if (isA<lduInterfaceField>(bpsi[patchI]))
                    {
                        label globalPatchID = primitiveMesh_.patchMap()[i][patchI];
                        if (isA<cyclicLduInterfaceField>(bpsi[patchI]))
                        {
                            interfaces_.set
                            (
                                globalPatchID,
                                new cyclicFvPatchField<scalar>
                                (
                                    refCast<const fvPatch>
                                    (
                                        primitiveMesh_.interfaces()[globalPatchID]
                                    ),
                                    bpsi[patchI].internalField()
                                )
                            );
                        }
                        else
                        {
                            interfaces_.set
                            (
                                globalPatchID,
                                &refCast<const lduInterfaceField>(bpsi[patchI])
                            );
                        }
                    }
                }
            }
        }
    }

    // Transfer fields from matrices to local assembly
    // lower, upper, diag, source
    transferFieldsAndClean();

    for (label i=0; i < nMatrix_; i++)
    {
        forAll (psis_[i].mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = psis_[i].mesh().boundaryMesh()[patchI];
            label globalPatchId = primitiveMesh_.patchMap()[i][patchI];

            if (globalPatchId != -1)
            {
                if (matrices_[i].internalCoeffs().set(patchI))
                {
                    internalCoeffs_.set
                    (
                        globalPatchId,
                        matrices_[i].internalCoeffs()[patchI]
                    );

                    boundaryCoeffs_.set
                    (
                        globalPatchId,
                        matrices_[i].boundaryCoeffs()[patchI]
                    );

                }
            }
            else
            {
                if (isA<mappedPatchBase>(pp))
                {
                    volScalarField::Boundary& tbf =
                        matrices_[i].mesh().thisDb().lookupObjectRef
                        <
                            volScalarField
                        >("T").boundaryFieldRef();

                    tbf[patchI].manipulateMatrix
                    (
                        *this,
                        primitiveMesh_.faceBoundMap()[i][patchI],
                        primitiveMesh_.cellOffsets()[i]
                    );

                    //updateCoeffs(i, patchI);
                }
                else
                {
                    FatalErrorInFunction
                        << "Local patch  : " << pp.index() << nl
                        << " it is not a mappedBased patch ID" << nl
                        << " in mesh : " << i
                        << exit(FatalError);
                }
            }
        }
    }
}


void Foam::fvMatrixAssemble::addFvMatrix(const fvMatrix<scalar>& matrix)
{
    matrices_.append(matrix);
    psis_.append(&const_cast<volScalarField&>(matrix.psi()));
    nMatrix_++;
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMatrixAssemble::solve
(
    const dictionary& solverControls
)
{
    // Modify nbrPatchID for cyclic,etc patches
    //updateNnbrPatchID(false);

    word solverType(solverControls.subDict(psiName_).lookup("solver"));

    if (solverType == GAMGSolver::typeName)
    {
        if (!faceAreasPtr_)
        {
            calFaceAreasCellVolumes();
        }
    }

    // Update lower/upper/diag/source to this assembled matrix
    // and delete individual matrices
    update();

    scalarField saveDiag(diag());

    addBoundaryDiag(diag());

    addBoundarySource(source_, false);

    scalarField psi(lduAddr().size(), 0.0);

    forAll (psis_, i)
    {
        label cellOffset = primitiveMesh_.cellOffsets()[i];

        forAll (psis_[i], localCellI)
        {
            psi[cellOffset + localCellI] = psis_[i][localCellI];
        }
    }

    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psiName_,
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces_,
        solverControls.subDict(psiName_)
    )->solve(psi, source_);

    forAll (psis_, i)
    {
        label cellOffset = primitiveMesh_.cellOffsets()[i];

        forAll (psis_[i], localCellI)
        {
            psis_[i][localCellI] = psi[localCellI + cellOffset];
        }
    }

    if (solverPerformance::debug)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    // Modify nbrPatchID for cyclic,etc patches
    //updateNnbrPatchID(true);

    // This will evaluate all BC's. The mapped BC should not alter the matrix.
    forAll (psis_, i)
    {
        psis_[i].correctBoundaryConditions();
        psis_[i].mesh().setSolverPerformance(psiName_, solverPerf);
    }

    clear();

    return solverPerf;
}


void Foam::fvMatrixAssemble::clear()
{
    matrices_.clear();
    psis_.clear();
    lower().clear();
    upper().clear();
    diag().clear();
    nMatrix_ = 0;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMatrixAssemble::~fvMatrixAssemble()
{}


// ************************************************************************* //
