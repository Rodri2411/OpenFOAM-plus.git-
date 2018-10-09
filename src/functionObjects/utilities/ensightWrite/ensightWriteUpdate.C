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

#include "ensightWrite.H"
#include "dictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::searchableSurfaces>
Foam::functionObjects::ensightWrite::loadSearchableSurfaces
(
    const dictionary& dict,
    const word& dictName
) const
{
    const auto finder = dict.csearch(dictName, keyType::LITERAL);

    if (!finder.isDict())
    {
        return nullptr;
    }

    auto surfsPtr = autoPtr<searchableSurfaces>::New
    (
        IOobject
        (
            "abc",              // dummy name
            time_.constant(),   // instance
            "triSurface",       // local
            time_,              // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        finder.dict(),
        true // singleRegionName
    );


    // Check/warn for non-enclosed

    auto& surfs = *surfsPtr;

    const label nsurfs = surfs.size();

    labelList oldToNew(nsurfs);

    label ngood = 0, badpos = nsurfs-1;

    forAll(surfs, surfi)
    {
        const searchableSurface& s = surfs[surfi];
        if (s.hasVolumeType())
        {
            oldToNew[surfi] = ngood;
            ++ngood;
        }
        else
        {
            // Bad. move to the end
            oldToNew[surfi] = badpos;
            --badpos;

            WarningInFunction
                << nl << "The surface '" << s.name() << "' of type '"
                << s.type() << "' appears to be unclosed"
                << nl << endl;
        }
    }

    if (ngood < nsurfs)
    {
        surfs.reorder(oldToNew);
        surfs.resize(ngood);
    }

    if (surfs.empty())
    {
        return nullptr;
    }

    return surfsPtr;
}


bool Foam::functionObjects::ensightWrite::updateSubset
(
    fvMeshSubset& subsetter
) const
{
    if (selectZones_.empty() && !selectEnclosed_.valid() && bounds_.empty())
    {
        return false;
    }

    const fvMesh& baseMesh = subsetter.baseMesh();

    const auto& cellCentres = baseMesh.C();

    const label len = baseMesh.nCells();

    bitSet cellsToSelect;

    // Mark all cells with the enclosed volumes
    if (selectEnclosed_.valid())
    {
        cellsToSelect.resize(len);

        List<volumeType> volTypes;
        for (const searchableSurface& s : *selectEnclosed_)
        {
            if (s.hasVolumeType())
            {
                s.getVolumeType(cellCentres, volTypes);

                for (label celli=0; celli < len; ++celli)
                {
                    if (volTypes[celli] == volumeType::INSIDE)
                    {
                        cellsToSelect.set(celli);
                    }
                }
            }
        }
    }


    // The intersection with cellZones
    if (selectZones_.size())
    {
        bitSet zoned(len);

        for (const label zonei : baseMesh.cellZones().indices(selectZones_))
        {
            const labelList& cellIds = baseMesh.cellZones()[zonei];

            zoned.set(cellIds);
        }

        if (cellsToSelect.empty())
        {
            // No enclosing volumes => use the whole mesh
            cellsToSelect = std::move(zoned);
        }
        else
        {
            cellsToSelect &= zoned;
            zoned.clear();
        }
    }


    // Clip with bounding box
    if (!bounds_.empty())
    {
        if (cellsToSelect.empty())
        {
            cellsToSelect.resize(len);

            for (label celli=0; celli < len; ++celli)
            {
                const point& cc = cellCentres[celli];

                if (bounds_.contains(cc))
                {
                    cellsToSelect.set(celli);
                }
            }
        }
        else
        {
            for (const label celli : cellsToSelect)
            {
                const point& cc = cellCentres[celli];

                if (!bounds_.contains(cc))
                {
                    cellsToSelect.unset(celli);
                }
            }
        }
    }

    subsetter.setCellSubset(cellsToSelect);

    return true;
}


#if 0
Foam::labelList Foam::functionObjects::ensightWrite::getSelectedPatches
(
    const polyBoundaryMesh& patches
) const
{
    DynamicList<label> patchIDs(patches.size());

    for (const polyPatch& pp : patches)
    {
        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (isType<processorPolyPatch>(pp))
        {
            break; // No processor patches
        }

        if
        (
            selectPatches_.size()
          ? selectPatches_.match(pp.name())
          : true
        )
        {
            patchIDs.append(pp.index());
        }
    }

    return patchIDs.shrink();
}
#endif


bool Foam::functionObjects::ensightWrite::update()
{
    if (meshState_ == polyMesh::UNCHANGED)
    {
        return false;
    }

    // This is heavy-handed, but with a bounding-box limited sub-mesh,
    // we don't readily know if the updates affect the subsetted mesh.

    // if (meshSubset_.hasSubMesh())
    // {
    //     ensMesh_.clear();
    //     meshSubset_.clear();
    // }
    // else if (ensMesh_.valid())
    // {
    //     ensMesh_->expire();
    // }

    meshSubset_.clear();

    updateSubset(meshSubset_);

    meshState_ = polyMesh::UNCHANGED;

    if (!ensMesh_.valid())
    {
        ensMesh_.reset(new ensightMesh(meshSubset_.mesh(), writeOpts_));
    }
    else if (ensMesh_().needsUpdate())
    {
        ensMesh_().correct();
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::readSelection(const dictionary& dict)
{
    // Ensure consistency
    ensMesh_.clear();

    meshSubset_.clear();
    meshState_ = polyMesh::TOPO_CHANGE;

    selectZones_.clear();
    if (!dict.readIfPresent("zones", selectZones_) && dict.found("zone"))
    {
        selectZones_.resize(1);
        dict.readEntry("zone", selectZones_.first());
    }

    selectFields_.clear();
    dict.readEntry("fields", selectFields_);
    selectFields_.uniq();

    selectEnclosed_.clear();
    selectEnclosed_ = loadSearchableSurfaces(dict, "geometry");

    bounds_.clear();
    dict.readIfPresent("bounds", bounds_);

    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh&)
{
    meshState_ = polyMesh::TOPO_CHANGE;
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh&)
{
    // Only move to worse states
    if (meshState_ == polyMesh::UNCHANGED)
    {
        meshState_ = polyMesh::POINTS_MOVED;
    }
}


// ************************************************************************* //
