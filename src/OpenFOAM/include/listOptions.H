argList::addBoolOption
(
    "list-switches",
    "List switches declared in libraries but not set in etc/controlDict"
);
argList::addBoolOption
(
    "list-registered-switches",
    "List switches registered for run-time modification"
);
argList::addBoolOption
(
    "list-unset-switches",
    "List switches declared in libraries but not set in etc/controlDict"
);

#ifdef fvPatchField_H
argList::addBoolOption
(
    "list-scalar-bcs",
    "List scalar field boundary conditions (fvPatchField<scalar>)"
);
argList::addBoolOption
(
    "list-vector-bcs",
    "List vector field boundary conditions (fvPatchField<vector>)"
);
#endif

#ifdef functionObject_H
argList::addBoolOption
(
    "list-function-objects",
    "List function objects"
);
#endif

#ifdef fvOption_H
argList::addBoolOption
(
    "list-fvoptions",
    "List fvOptions"
);
#endif

#if defined(turbulentTransportModel_H) || defined(turbulentFluidThermoModel_H)
argList::addBoolOption
(
    "list-turbulence-models",
    "List turbulenceModels"
);
#endif
