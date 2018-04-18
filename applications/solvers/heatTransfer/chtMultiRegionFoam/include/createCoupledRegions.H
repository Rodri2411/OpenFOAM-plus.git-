fvSolution solutionDict(runTime);

bool coupled = false;

solutionDict.lookup("coupled") >> coupled;

autoPtr<lduPrimitiveMeshAssemble> meshPtr;
autoPtr<fvMatrixAssemble> fvMatrixAssemblePtr;

if (coupled)
{
    UPtrList<fvMesh> meshes(fluidRegions.size() + solidRegions.size());

    label iMesh = 0;

    forAll(fluidRegions, i)
    {
        meshes.set(iMesh, &fluidRegions[i]);
        iMesh++;
    }

    forAll(solidRegions, i)
    {
        meshes.set(iMesh, &solidRegions[i]);
        iMesh++;
    }

    Info<< "Creating coupled mesh.."<< nl << endl;

    IOobject io
    (
        "assembleLdu",
        runTime.timeName(),
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    );

    meshPtr.reset
    (
        new lduPrimitiveMeshAssemble(meshes, io)
    );

    Info<< "Creating coupled matrix.."<< nl << endl;

    fvMatrixAssemblePtr.reset
    (
        new fvMatrixAssemble
        (
            meshPtr(),
            dimEnergy/dimTime,
            "he"
        )
    );

    Info<< "Created coupled matrix.."<< nl << endl;
}
