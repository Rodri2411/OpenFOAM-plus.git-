/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::sampledSurfaces::writeSurface
(
    sampledSurfaceProxy& proxy,
    const word& fieldName,
    const Field<Type>& values
)
{
    if (proxy.staged())
    {
        // Trigger any changes
        formatter_->updateMesh(proxy.surface().name());
        proxy.handled();
    }

    Info<<"sampled " << fieldName << " on "
        << proxy.name() << " = " << proxy << nl;


    fileName sampleFile;
    if (Pstream::parRun())
    {
        // Collect values from all processors
        List<Field<Type>> gatheredValues(Pstream::nProcs());
        gatheredValues[Pstream::myProcNo()] = values;
        Pstream::gatherList(gatheredValues);

        if (Pstream::master())
        {
            // Combine values into single field
            Field<Type> allValues
            (
                ListListOps::combine<Field<Type>>
                (
                    gatheredValues,
                    accessOp<Field<Type>>()
                )
            );

            // Renumber (point data) to correspond to merged points
            proxy.renumberPointData(allValues);

            // Skip surface without faces (eg, a failed cut-plane)
            if (proxy.hasGeometry())
            {
                sampleFile = formatter_->write
                (
                    proxy.surface().name(),
                    proxy.geometry(),
                    fieldName,
                    allValues,
                    proxy.surface().interpolate()
                );
            }
        }

        Pstream::scatter(sampleFile);
    }
    else
    {
        // Skip surface without faces (eg, a failed cut-plane)
        if (proxy.hasGeometry())
        {
            sampleFile = formatter_->write
            (
                proxy.surface().name(),
                proxy.geometry(),
                fieldName,
                values,
                proxy.surface().interpolate()
            );
        }
    }

    if (sampleFile.size())
    {
        dictionary propsDict;
        propsDict.add("file", sampleFile);
        setProperty(fieldName, propsDict);
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    // Formatter timeValue/timeName already updated in caller

    // The sampler/interpolator for this field
    autoPtr<interpolation<Type>> interp;

    const word& fieldName = vField.name();

    const auto& surfs = surfaces();

    forAll(surfs, surfi)
    {
        const sampledSurface& s = surfs[surfi];
        sampledSurfaceProxy& proxy = proxies_[surfi];

        Field<Type> values;

        if (s.interpolate())
        {
            if (!interp)
            {
                interp = interpolation<Type>::New(sampleNodeScheme_, vField);
            }

            values = s.interpolate(*interp);
        }
        else
        {
            if (!interp)
            {
                interp = interpolation<Type>::New(sampleFaceScheme_, vField);
            }

            values = s.sample(*interp);
        }

        writeSurface<Type>(proxy, fieldName, values);
    }
}


template<class Type>
void Foam::sampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    // Formatter timeValue/timeName already updated in caller

    const word& fieldName = sField.name();
    const fileName outputDir = outputPath_/sField.time().timeName();

    const auto& surfs = surfaces();

    forAll(surfs, surfi)
    {
        const sampledSurface& s = surfs[surfi];
        sampledSurfaceProxy& proxy = proxies_[surfi];

        Field<Type> values(s.sample(sField));
        writeSurface<Type>(proxy, fieldName,  values);
    }
}


template<class GeoField>
void Foam::sampledSurfaces::sampleAndWrite(const IOobjectList& objects)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames(GeoField::typeName, fieldSelection_);
    }
    else
    {
        fieldNames = mesh_.thisDb().sortedNames<GeoField>(fieldSelection_);

        writeOriginalIds();
    }

    for (const word& fieldName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "sampleAndWrite: " << fieldName << endl;
        }

        if (loadFromFiles_)
        {
            const GeoField fld
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    mesh_,
                    IOobject::MUST_READ
                ),
                mesh_
            );

            sampleAndWrite(fld);
        }
        else
        {
            sampleAndWrite
            (
                mesh_.thisDb().lookupObject<GeoField>(fieldName)
            );
        }
    }
}


// ************************************************************************* //
