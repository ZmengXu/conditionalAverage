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
Foam::tmp<Foam::Field<Type>>
Foam::conditionalAverage::filterField
(
    const Field<Type>& field
) const
{
    if (isNull(cellIDs()))
    {
        return field;
    }
    else
    {
        return tmp<Field<Type>>(new Field<Type>(field, cellIDs()));
    }
}

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
					fvMeshFunctionObject::mesh_.time().timeName(),
					fvMeshFunctionObject::mesh_,
					IOobject::NO_READ,
					IOobject::NO_WRITE,
					false
				),
				fvMeshFunctionObject::mesh_,
				dimensionedScalar("direction", dimLength, 0)
			);
			
			const volVectorField centres = fvMeshFunctionObject::mesh_.C();
			
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
				filterField
				(
					tempconditionalField
				)
			);	
		}
		else
		{
			if (loadFromFiles_)
			{
                volScalarField tempconditionalField
                (
                    IOobject
                    (
                        conditionalFieldName,
                        fvMeshFunctionObject::mesh_.time().timeName(),
                        fvMeshFunctionObject::mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    fvMeshFunctionObject::mesh_
                );
				conditionalFields_.set
				(
					conditionalFieldi,
                    filterField
                    (
                        tempconditionalField
                    )
				);
			}
			else
			{
				if(fvMeshFunctionObject::mesh_.foundObject<volScalarField>(conditionalFieldName))
				{
					conditionalFields_.set
					(
						conditionalFieldi,
						filterField
						(
							fvMeshFunctionObject::mesh_.lookupObject<volScalarField>(conditionalFieldName)
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

		const scalarField& conditionalField = conditionalFields_[conditionalFieldi];
		
		scalarField weightedAveragedField_ = fvMeshFunctionObject::mesh_.V();
		if(weightedAveragedFieldName_ == "meshV")
		{
			forAll(fvMeshFunctionObject::mesh_.C(),celli)
			{
				weightedAveragedField_[celli] =  fvMeshFunctionObject::mesh_.V()[celli];
			}
		}
		else if (weightedAveragedFieldName_ == "rhoMeshV")
		{
			const volScalarField rho_ = fvMeshFunctionObject::mesh_.lookupObject<volScalarField>("rho");
			forAll(fvMeshFunctionObject::mesh_.C(),celli)
			{
				weightedAveragedField_[celli] =  fvMeshFunctionObject::mesh_.V()[celli]*rho_[celli];
			}
		}
		else if (weightedAveragedFieldName_ == "none")
		{
			forAll(fvMeshFunctionObject::mesh_.C(),celli)
			{
				weightedAveragedField_[celli] =  1.;
			}			
		}

		//scalarList totalCounts_(nBins_,scalar(0.0));// Only for of6, need scalar(0)
		//scalarList localCellCounts(nBins_,scalar(0.0));// of4 and 7 can use, scalarList totalCounts_(nBins_,0);
		scalarList localWeightedAveragedFields_(nBins_,scalar(0.0));
		scalarList totalWeightedAveragedFields_(nBins_,scalar(0.0));
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
			new scalarList(0)//Init as size 0, append when there are some data in this bin
		);
		
		forAll(localAveragedFields, averagedFieldi)
		{
			localAveragedFields[averagedFieldi] = Field<T>(nBins_);
			averagedFieldsOutput_[averagedFieldi] = Field<T>(0);//Init as size 0, append when there are some data in this bin
			localAveragedFields[averagedFieldi] = 0.0*localAveragedFields[averagedFieldi];
			averagedFieldsOutput_[averagedFieldi] = 0.0*averagedFieldsOutput_[averagedFieldi];
		}
		scalar start = gMin(conditionalField);
		scalar end = gMax(conditionalField);
        //- This if will be skiped if minF_ and maxF_ are not given
        if( maxF_ > minF_ )
        {
            start = max( minF_, start);
            end   = min( maxF_, end  );
        }
		const scalar offset = (end - start)/(nBins_ - 1);

		forAll( conditionalField, trackCelli )
		{
			if((conditionalField[trackCelli]>=start)&&(conditionalField[trackCelli]<=end))
			{
				//get the in iBin_th localSamplingCells, assuming it is averaged distributed
				label iBin = 0;
				if(offset!=0) iBin = label((conditionalField[trackCelli] - start)/offset);
				
				//localCellCounts[iBin] ++;
				localWeightedAveragedFields_[iBin] += weightedAveragedField_[trackCelli];
				forAll( averagedFields, averagedFieldi )
				{
					// New averaged value, avoid exceed value
				//	localAveragedFields[averagedFieldi][iBin] =
				//			localAveragedFields[averagedFieldi][iBin]*(localCellCounts[iBin]-1)/localCellCounts[iBin]
				//			+ weightedAveragedField_[trackCelli]*averagedFields[averagedFieldi][trackCelli]/localCellCounts[iBin];
					localAveragedFields[averagedFieldi][iBin] += weightedAveragedField_[trackCelli]*averagedFields[averagedFieldi][trackCelli];
				}
			}
		}
		
		// for weightedAveragedField
		forAll( averagedFields, averagedFieldi )
		{
			// New averaged value, avoid exceed value
			forAll(localAveragedFields[averagedFieldi], iBin)
			{
				if(localWeightedAveragedFields_[iBin] != 0)
				{
					localAveragedFields[averagedFieldi][iBin] /= localWeightedAveragedFields_[iBin];
				}
			}
		}
	
		//- Get the value from all of the processors
		//forAll( conditionalFieldOutputs_[conditionalFieldi], iBin )
        for(label iBin = 0; iBin < nBins_; iBin++ )
		{
			//totalCounts_[iBin] = localCellCounts[iBin];
			totalWeightedAveragedFields_[iBin] = localWeightedAveragedFields_[iBin];
			reduce( totalWeightedAveragedFields_[iBin], sumOp<scalar>()); 
			//- if there is no data in this region, did not append the value for conditionalFieldOutputs_ and averagedFieldsOutput_
			if(totalWeightedAveragedFields_[iBin] > 0)
			{
                conditionalFieldOutputs_[conditionalFieldi].append(start + iBin*offset);
				forAll( averagedFields, averagedFieldi )
				{
					T totalAveraged = (localWeightedAveragedFields_[iBin]/totalWeightedAveragedFields_[iBin])*localAveragedFields[averagedFieldi][iBin];
					reduce( totalAveraged, sumOp<T>()); 
					averagedFieldsOutput_[averagedFieldi].append(totalAveraged);
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
                        fvMeshFunctionObject::mesh_.time().timeName(),
                        fvMeshFunctionObject::mesh_,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    fvMeshFunctionObject::mesh_
                );

				averagedFields.set
				(
					fieldi,
					filterField(vf)// use region to filter part of the volScalarField
				);
            }
            else
            {
				averagedFields.set
				(
					fieldi,
					filterField// use region to filter part of the volScalarField
					(
						fvMeshFunctionObject::mesh_.lookupObject
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
					List<point>(conditionalFieldOutputs_[conditionalFieldi].size()),
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
					outputPath_/fvMeshFunctionObject::mesh_.time().timeName(),
					fields.formatter()
				);
			}

        }
    }
}


// ************************************************************************* //
