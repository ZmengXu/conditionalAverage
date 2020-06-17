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
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::conditionalAverage::writeSampleFile
(
    const coordSet& masterSampleSet,
    const PtrList<Field<Type>>& masterFields,
	const wordList nameList,
    const fileName& timeDir,
    const writer<Type>& formatter
)
{
    wordList valueSetNames(masterFields.size());
    List<const Field<Type>*> valueSets(masterFields.size());

    forAll(masterFields, fieldi)
    {
        valueSetNames[fieldi] = nameList[fieldi];
        valueSets[fieldi] = &masterFields[fieldi];
    }

    fileName fName
    (
        timeDir/formatter.getFileName(masterSampleSet, valueSetNames)
    );

    OFstream ofs(fName);
    if (ofs.opened())
    {
        formatter.write
        (
            masterSampleSet,
            valueSetNames,
            valueSets,
            ofs
        );
    }
    else
    {
        WarningInFunction
            << "File " << ofs.name() << " could not be opened. "
            << "No data will be written" << endl;
    }
}


template<class T>
void Foam::conditionalAverage::combineSampledValues
(
    PtrList<Field<T>>& averagedFields,
    PtrList<PtrList<Field<T>>>& masterFields
)
{			
	forAll(conditionalFields_, conditionalFieldi)
	{
		word conditionalFieldName(conditionalFieldNames_[conditionalFieldi]);

		if
		(
			(conditionalFieldName == "x") ||
			(conditionalFieldName == "y") ||
			(conditionalFieldName == "z") ||
			(conditionalFieldName == "direction")
		)
		{
			volScalarField tempconditionalField
			(
				IOobject
				(
					conditionalFieldName,
					mesh_.time().timeName(),
					mesh_,
					IOobject::NO_READ,
					IOobject::NO_WRITE,
					false
				),
				mesh_,
				dimensionedScalar("direction", dimLength, 0)
			);
			
			const volVectorField centres = mesh_.C();
			
			if(conditionalFieldName == "x") tempconditionalField = centres.component(vector::X);
			else if(conditionalFieldName == "y") tempconditionalField = centres.component(vector::Y);
			else if(conditionalFieldName == "z") tempconditionalField = centres.component(vector::Z);
			else
			{
				vector directionVector(dict_.lookup("directionVector"));
				directionVector = directionVector/(mag(directionVector)+SMALL);
				Info << "The normalized directionVector is:" << directionVector << endl;

				tempconditionalField = (directionVector & centres);	
			}
			
			conditionalFields_.set
			(
				conditionalFieldi,
				new volScalarField
				(
					tempconditionalField
				)
			);	
		}
		else
		{
			if (loadFromFiles_)
			{
				conditionalFields_.set
				(
					conditionalFieldi,
					new volScalarField
					(
						IOobject
						(
							conditionalFieldName,
							mesh_.time().timeName(),
							mesh_,
							IOobject::MUST_READ,
							IOobject::NO_WRITE,
							false
						),
						mesh_
					)
				);
			}
			else
			{
				if(mesh_.foundObject<volScalarField>(conditionalFieldName))
				{
					conditionalFields_.set
					(
						conditionalFieldi,
						new volScalarField
						(
							mesh_.lookupObject<volScalarField>(conditionalFieldName)
						)
					);
				}
				else
				{
					FatalErrorInFunction
						<< "Cannot find volScalarField "
						<< conditionalFieldName
						<< " for a conditionalField." << nl
						<< "It is not in objectRegistry or not a volScalarField."
						<< abort(FatalError);
				}
			}
		}

		const volScalarField& conditionalField = conditionalFields_[conditionalFieldi];

		scalarList totalCounts_(nBins_,scalar(0.0));// Only for of6, need scalar(0)
		scalarList localCellCounts(nBins_,scalar(0.0));// of4 and 7 can use, scalarList totalCounts_(nBins_,0);
		List<Field<T>> localAveragedFields(averagedFields.size());
		List<Field<T>> averagedFieldsOutput_(averagedFields.size());

		masterFields.set
		(
			conditionalFieldi,
			new PtrList<Field<T>>(averagedFields.size())
		);

		conditionalFieldOutputs_.set
		(
			conditionalFieldi,
			new scalarList(nBins_)
		);
		forAll(localAveragedFields, averagedFieldi)
		{
			localAveragedFields[averagedFieldi] = Field<T>(nBins_);
			averagedFieldsOutput_[averagedFieldi] = Field<T>(nBins_);
			localAveragedFields[averagedFieldi] = 0.0*localAveragedFields[averagedFieldi];
			averagedFieldsOutput_[averagedFieldi] = 0.0*averagedFieldsOutput_[averagedFieldi];
		}
		const scalar start = gMin(conditionalField);
		const scalar end = gMax(conditionalField);
		const scalar offset = (end - start)/nBins_;

		forAll( conditionalField, trackCelli )
		{
			if((conditionalField[trackCelli]>=start)&&(conditionalField[trackCelli]<=end))
			{
				//get the in iBin_th localSamplingCells, assuming it is averaged distributed
				label iBin = 0;
				if(offset!=0) iBin = label((conditionalField[trackCelli] - start)/offset);
				
				localCellCounts[iBin] ++;
				
				forAll( averagedFields, averagedFieldi )
				{
					// New averaged value, avoied exceed value
					localAveragedFields[averagedFieldi][iBin] =
							localAveragedFields[averagedFieldi][iBin]*(localCellCounts[iBin]-1)/localCellCounts[iBin]
							+ averagedFields[averagedFieldi][trackCelli]/localCellCounts[iBin];
				}
			}
		}
	
		//- Get the value from all of the processors
		forAll( conditionalFieldOutputs_[conditionalFieldi], iBin )
		{
			conditionalFieldOutputs_[conditionalFieldi][iBin] = start + iBin*offset;
			totalCounts_[iBin] = localCellCounts[iBin];
			reduce( totalCounts_[iBin], sumOp<scalar>()); 
			//- if there is no data in this region, did not change averagedFieldsOutput_, use 0 instead
			if(totalCounts_[iBin] > 0)
			{
				forAll( averagedFields, averagedFieldi )
				{
					T totalAveraged = (localCellCounts[iBin]/totalCounts_[iBin])*localAveragedFields[averagedFieldi][iBin];
					reduce( totalAveraged, sumOp<T>()); 
					averagedFieldsOutput_[averagedFieldi][iBin] = totalAveraged;
				}
			}
		}

		forAll( averagedFields, averagedFieldi )
		{
			//- In the condition of "conditionalFieldi", the averaged field "averagedFieldi"		
			//- masterFields[averagedFieldi].setname() = averagedFields[averagedFieldi].name();
			//- masterFields[conditionalFieldi] = ;
			masterFields[conditionalFieldi].set
			(
				averagedFieldi,
				new Field<T>(averagedFieldsOutput_[averagedFieldi])
			);
		}
	}
}


template<class Type>
void Foam::conditionalAverage::sampleAndWrite
(
    fieldGroup<Type>& fields
)
{
    if (fields.size())
    {
        // Create or use existing writer
        // fields has two properties, DynamicList and formatter
        // here we check the formatter
        if (fields.formatter.empty())
        {
            fields = writeFormat_;
        }

        // Storage for interpolated values
        PtrList<Field<Type>> averagedFields(fields.size());

        forAll(fields, fieldi)
        {
            if (Pstream::master() && verbose_)
            {
                Pout<< "conditionalAverage::sampleAndWrite: "
                    << fields[fieldi] << endl;
            }

            if (loadFromFiles_)
            {
                GeometricField<Type, fvPatchField, volMesh> vf
                (
                    IOobject
                    (
                        fields[fieldi],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh_
                );

				averagedFields.set
				(
					fieldi,
					new Field<Type>(vf)
				);
            }
            else
            {
				averagedFields.set
				(
					fieldi,
					new Field<Type>
					(
						mesh_.lookupObject
						<GeometricField<Type, fvPatchField, volMesh>>
						(fields[fieldi])
					)
				);
            }
        }

        // Combine sampled fields from processors.
        // Note: only master results are valid
        PtrList<PtrList<Field<Type>>> masterFields(conditionalFieldNames_.size());
		combineSampledValues(averagedFields, masterFields);

        if (Pstream::master())
        {
			forAll(conditionalFields_, conditionalFieldi)
			{
				coordSet masterCoordSets
				(
					conditionalFieldNames_[conditionalFieldi]+"_conditionalwith",
					"distance",
					List<point>(nBins_),
					conditionalFieldOutputs_[conditionalFieldi]
				);

				wordList nameList(masterFields[conditionalFieldi].size());

				forAll(nameList, fieldi)
				{
					nameList[fieldi] = fields[fieldi];
				}
		
				writeSampleFile
				(
					masterCoordSets,
					masterFields[conditionalFieldi],
					nameList,
					outputPath_/mesh_.time().timeName(),
					fields.formatter()
				);
			}

        }
    }
}


// ************************************************************************* //
