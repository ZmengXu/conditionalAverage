/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "conditionalAverage.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conditionalAverage, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        conditionalAverage,
        dictionary
    );
}

bool Foam::conditionalAverage::verbose_ = true;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conditionalAverage::conditionalAverage
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    loadFromFiles_(false),
    outputPath_(fileName::null),
    dict_(dict),
    conditionalFieldNames_(dict.lookup("conditionalFields")),
    conditionalFields_(conditionalFieldNames_.size()),
	weightedAveragedFieldName_(dict.lookup("weightedAveragedField")),
    nBins_(readLabel(dict.lookup("nBins"))),
	maxF_(dict.lookupOrDefault<scalar>("maxValue", 0.0)),
	minF_(dict.lookupOrDefault<scalar>("minValue", 0.0)),
	conditionalFieldOutputs_(conditionalFieldNames_.size()),
	fieldSelection_(dict.lookup("averagedFields")),
    writeFormat_(dict.lookup("writeFormat"))
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/"postProcessing"/name;
    }
    else
    {
        outputPath_ = mesh_.time().path()/"postProcessing"/name;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

	clearFieldGroups();

    if (Pstream::master() && debug)
    {
		Pout<< "conditionalFields are:" << conditionalFieldNames_ << endl;
        Pout<< "conditionalAverage fields:" << fieldSelection_ << endl;
    }
}


Foam::conditionalAverage::conditionalAverage
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObject(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    dict_(dict),
    conditionalFieldNames_(dict.lookup("conditionalFields")),
    conditionalFields_(conditionalFieldNames_.size()),
	weightedAveragedFieldName_(dict.lookup("weightedAveragedField")),
    nBins_(readLabel(dict.lookup("nBins"))),
	maxF_(dict.lookupOrDefault<scalar>("maxValue", 0.0)),
	minF_(dict.lookupOrDefault<scalar>("minValue", 0.0)),
	conditionalFieldOutputs_(conditionalFieldNames_.size()),
	fieldSelection_(dict.lookup("averagedFields")),
    writeFormat_(dict.lookup("writeFormat"))
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/"postProcessing"/name;
    }
    else
    {
        outputPath_ = mesh_.time().path()/"postProcessing"/name;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

	clearFieldGroups();

    if (Pstream::master() && debug)
    {
		Pout<< "conditionalFields are:" << conditionalFieldNames_ << endl;
        Pout<< "conditionalAverage fields:" << fieldSelection_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conditionalAverage::~conditionalAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conditionalAverage::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::conditionalAverage::execute()
{
    return true;
}


bool Foam::conditionalAverage::write()
{
	const label nFields = classifyFields();

	if (Pstream::master())
	{
		if (debug)
		{
			Pout<< "timeName = " << mesh_.time().timeName() << nl
				<< "scalarFields    " << scalarFields_ << nl
				<< "vectorFields    " << vectorFields_ << nl
				<< "sphTensorFields " << sphericalTensorFields_ << nl
				<< "symTensorFields " << symmTensorFields_ <<nl
				<< "tensorFields    " << tensorFields_ <<nl;
		}

		if (nFields)
		{
			if (debug)
			{
				Pout<< "Creating directory "
					<< outputPath_/mesh_.time().timeName()
						<< nl << endl;
			}

			mkDir(outputPath_/mesh_.time().timeName());
		}
	}

	if (nFields)
	{
		sampleAndWrite(scalarFields_);
		sampleAndWrite(vectorFields_);
		sampleAndWrite(sphericalTensorFields_);
		sampleAndWrite(symmTensorFields_);
		sampleAndWrite(tensorFields_);
	}

    return true;
}


bool Foam::conditionalAverage::read(const dictionary& dict)
{    
    return true;
}

// ************************************************************************* //
