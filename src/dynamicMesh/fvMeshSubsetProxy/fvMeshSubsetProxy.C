/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "fvMeshSubsetProxy.H"
#include "cellSet.H"
#include "cellZone.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshSubsetProxy::fvMeshSubsetProxy(fvMesh& baseMesh)
:
    fvMeshSubsetProxy(baseMesh, boundBox::invertedBox)
{}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const boundBox& bounds
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(-1),
    type_(NONE),
    name_(),
    zones_(),
    bounds_(bounds)
{
    if (!bounds.empty())
    {
        type_ |= 0x100;   // Use bounding box
    }

    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    subsetType type,
    const word& name,
    const boundBox& bounds,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(name.empty() ? NONE : type),
    name_(name),
    zones_(),
    bounds_(bounds)
{
    if (type_ == ZONES)
    {
        // Ensure wordRes is populated for ZONES
        zones_.resize(1);
        zones_.first() = name_;
    }

    if (!bounds.empty())
    {
        type_ |= 0x100;   // Use bounding box
    }

    if (useSubMesh())
    {
        correct();
    }
}


Foam::fvMeshSubsetProxy::fvMeshSubsetProxy
(
    fvMesh& baseMesh,
    const wordRes& zones,
    const word& name,
    const boundBox& bounds,
    label exposedPatchId
)
:
    baseMesh_(baseMesh),
    subsetter_(baseMesh),
    exposedPatchId_(exposedPatchId),
    type_(ZONES),
    name_(name),
    zones_(zones),
    bounds_(bounds)
{
    if (!bounds.empty())
    {
        type_ |= 0x100;   // Use bounding box
    }

    if (useSubMesh())
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshSubsetProxy::bounds(const boundBox& bounds)
{
    if (!bounds.empty())
    {
        bounds_ = bounds;
        type_ |= 0x100;   // Use bounding box
    }
}


void Foam::fvMeshSubsetProxy::correct(bool verbose)
{
    if (type_ == NONE)
    {
        subsetter_.clear();
        return;
    }

    const label nCells = baseMesh_.nCells();

    bitSet selectedCells;

    if ((type_ & SET) != 0)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellSet " << name_ << endl;
        }

        cellSet cset(baseMesh_, name_);

        selectedCells.resize(nCells);
        for (const label idx : cset)
        {
            selectedCells.set(idx);
        }
    }
    else if ((type_ & ZONE) != 0)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZone " << name_ << endl;
        }

        selectedCells.resize(nCells);
        selectedCells.set(baseMesh_.cellZones()[name_]);
    }
    else if ((type_ & ZONES) != 0)
    {
        if (verbose)
        {
            Info<< "Subsetting mesh based on cellZones" << endl;
        }

        selectedCells = baseMesh_.cellZones().selection(zones_);
    }

    if (((type_ & 0x100) != 0) && !bounds_.empty())
    {
        // Accept/reject based on cell centres
        const auto& cc = baseMesh_.C();

        if (verbose)
        {
            Info<< "Subsetting mesh based on bounding box" << endl;
        }

        if ((type_ & ~0x100) == NONE)
        {
            // Select all
            selectedCells.resize(nCells, true);

            for (label celli=0; celli < nCells; ++celli)
            {
                if (!bounds_.contains(cc[celli]))
                {
                    selectedCells.unset(celli);
                }
            }
        }
        else
        {
            // Sub-selection
            for (const label celli : selectedCells)
            {
                if (!bounds_.contains(cc[celli]))
                {
                    selectedCells.unset(celli);
                }
            }
        }
    }

    subsetter_.setCellSubset(selectedCells, exposedPatchId_);
}


Foam::polyMesh::readUpdateState Foam::fvMeshSubsetProxy::readUpdate()
{
    const polyMesh::readUpdateState meshState = baseMesh_.readUpdate();

    if
    (
        meshState == polyMesh::TOPO_CHANGE
     || meshState == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        correct(true);
    }

    return meshState;
}


// ************************************************************************* //
