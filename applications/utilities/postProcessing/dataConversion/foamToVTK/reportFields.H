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
    Foam::foamToVtkReportFields

Description
    Collection of simple static methods for reporting field names
    by category, which is used by foamToVTK.

\*---------------------------------------------------------------------------*/

#ifndef foamToVtkReportFields_H
#define foamToVtkReportFields_H

#include "Ostream.H"
#include "IOobjectNames.H"
#include "wordList.H"
#include "UPtrList.H"

#include "areaFields.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class foamToVtkReportFields Declaration
\*---------------------------------------------------------------------------*/

struct foamToVtkReportFields
{
    template<class GeoField>
    static void print
    (
        const char* msg,
        Ostream& os,
        const UPtrList<const GeoField>& flds
    )
    {
        if (flds.size())
        {
            os  << msg;
            for (const auto& fld : flds)
            {
                os  << ' ' << fld.name();
            }
            os  << endl;
        }
    }


    static void print
    (
        const char* msg,
        Ostream& os,
        const wordList& fieldNames
    )
    {
        if (fieldNames.size())
        {
            os  << msg;
            for (const word& fieldName : fieldNames)
            {
                os  << ' ' << fieldName;
            }
            os  << endl;
        }
    }


    template<class FieldType>
    static void print
    (
        const char* msg,
        Ostream& os,
        const IOobjectNames& objNames
    )
    {
        print(msg, os, objNames.sortedNames<FieldType>());
    }


    //- Supported volume field types
    static void vols(Ostream& os, const IOobjectNames& objNames)
    {
        print<volScalarField>
        (
            "    volScalar    :", os, objNames
        );
        print<volVectorField>
        (
            "    volVector    :", os, objNames
        );
        print<volSphericalTensorField>
        (
            "    volSphTensor :", os, objNames
        );
        print<volSymmTensorField>
        (
            "    volSymTensor :", os, objNames
        );
        print<volTensorField>
        (
            "    volTensor    :", os, objNames
        );
    }


    //- Supported dimensioned field types
    static void dims(Ostream& os, const IOobjectNames& objNames)
    {
        print<volScalarField::Internal>
        (
            "    volScalar:Internal    :", os, objNames
        );
        print<volVectorField::Internal>
        (
            "    volVector:Internal    :", os, objNames
        );
        print<volSphericalTensorField::Internal>
        (
            "    volSphTensor:Internal :", os, objNames
        );
        print<volSymmTensorField::Internal>
        (
            "    volSymTensor:Internal :", os, objNames
        );
        print<volTensorField::Internal>
        (
            "    volTensor:Internal    :", os, objNames
        );
    }


    //- Supported point field types
    static void points(Ostream& os, const IOobjectNames& objNames)
    {
    }


    //- Supported area field types
    static void areas(Ostream& os, const IOobjectNames& objNames)
    {
        print<areaScalarField>
        (
            "    areaScalar    :", os, objNames
        );
        print<areaVectorField>
        (
            "    areaVector    :", os, objNames
        );
        print<areaSphericalTensorField>
        (
            "    areaSphTensor :", os, objNames
        );
        print<areaSymmTensorField>
        (
            "    areaSymTensor :", os, objNames
        );
        print<areaTensorField>
        (
            "    areaTensor    :", os, objNames
        );
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
