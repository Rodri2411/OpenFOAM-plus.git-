bool listOptions = false ;

if
(
    args.optionFound("list-switches")
)
{
    debug::listSwitches(args.optionFound("include-unset-switches"));
    listOptions = true;
}

if
(
    args.optionFound("list-registered-switches")
)
{
    debug::listRegisteredSwitches(args.optionFound("include-unset-switches"));
    listOptions = true;
}

#ifdef fvPatchField_H
if (args.optionFound("list-scalar-bcs"))
{
    Info<< "scalarBCs"
        << fvPatchField<scalar>::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}

if (args.optionFound("list-vector-bcs"))
{
    Info<< "vectorBCs"
        << fvPatchField<vector>::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef functionObject_H
if (args.optionFound("list-function-objects"))
{
    Info<< "functionObjects"
        << functionObject::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef fvOption_H
if (args.optionFound("list-fvoptions"))
{
    Info<< "fvOptions"
        << fv::option::dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

#ifdef turbulentTransportModel_H
if (args.optionFound("list-turbulence-models"))
{
    Info<< "Turbulence models"
        << incompressible::turbulenceModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "RAS models"
        << incompressible::RASModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "LES models"
        << incompressible::LESModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#elif defined(turbulentFluidThermoModel_H)
if (args.optionFound("list-turbulence-models"))
{
    Info<< "Turbulence models"
        << compressible::turbulenceModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "RAS models"
        << compressible::RASModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;

    Info<< "LES models"
        << compressible::LESModel::
           dictionaryConstructorTablePtr_->sortedToc()
        << endl;
    listOptions = true;
}
#endif

if (listOptions)
{
    exit(0);
}
