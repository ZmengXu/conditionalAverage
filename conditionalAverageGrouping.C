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
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::conditionalAverage::clearFieldGroups()
{
    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();
}


Foam::label Foam::conditionalAverage::appendFieldGroup
(
    const word& fieldName,
    const word& fieldType
)
{
    if (fieldType == volScalarField::typeName)
    {
        scalarFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volVectorField::typeName)
    {
        vectorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volSphericalTensorField::typeName)
    {
        sphericalTensorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volSymmTensorField::typeName)
    {
        symmTensorFields_.append(fieldName);
        return 1;
    }
    else if (fieldType == volTensorField::typeName)
    {
        tensorFields_.append(fieldName);
        return 1;
    }

    return 0;
}


Foam::label Foam::conditionalAverage::classifyFields()
{
    label nFields = 0;
    clearFieldGroups();

    if (loadFromFiles_)
    {
        // Check files for a particular time
        IOobjectList objects(fvMeshFunctionObject::mesh_, fvMeshFunctionObject::mesh_.time().timeName());
        wordList allFields = objects.sortedNames();

        forAll(fieldSelection_, i)
        {
            labelList indices = findStrings(fieldSelection_[i], allFields);

            if (indices.size())
            {
                forAll(indices, fieldi)
                {
                    const word& fieldName = allFields[indices[fieldi]];

                    nFields += appendFieldGroup
                    (
                        fieldName,
                        objects.find(fieldName)()->headerClassName()
                    );
                }
            }
            else
            {
                WarningInFunction
                    << "Cannot find field file matching "
                    << fieldSelection_[i] << endl;
            }
        }
    }
    else
    {
        // Check currently available fields
        wordList allFields = fvMeshFunctionObject::mesh_.sortedNames();
        labelList indices = findStrings(fieldSelection_, allFields);

        forAll(fieldSelection_, i)
        {
            labelList indices = findStrings(fieldSelection_[i], allFields);

            if (indices.size())
            {
                forAll(indices, fieldi)
                {
                    const word& fieldName = allFields[indices[fieldi]];

                    nFields += appendFieldGroup
                    (
                        fieldName,
                        fvMeshFunctionObject::mesh_.find(fieldName)()->type()
                    );
                }
            }
            else
            {
                WarningInFunction
                    << "Cannot find registered field matching "
                    << fieldSelection_[i] << endl;
            }
        }
    }

    return nFields;
}


// ************************************************************************* //
