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

Class
    Foam::conditionalAverage

Description
    Set of sets to sample.
    Call conditionalAverage.write() to sample&write files.

SourceFiles
    conditionalAverage.C

\*---------------------------------------------------------------------------*/

#ifndef conditionalAverage_H
#define conditionalAverage_H

#include "functionObject.H"
#include "volFieldsFwd.H"
#include "coordSet.H"
#include "writer.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class fvMesh;

/*---------------------------------------------------------------------------*\
                         Class conditionalAverage Declaration
\*---------------------------------------------------------------------------*/

class conditionalAverage
:
    public functionObject
{
    
    // Private classes

        //- Class used for grouping field types
        template<class Type>
        class fieldGroup
        :
            public DynamicList<word>
        {
        public:

            //- The set formatter
            autoPtr<writer<Type>> formatter;

            //- Construct null
            fieldGroup()
            :
                DynamicList<word>(0),
                formatter(nullptr) // change for of1806, from formatter(NULL)
            {}

            //- Construct for a particular format
            fieldGroup(const word& writeFormat)
            :
                DynamicList<word>(0),
                formatter(writer<Type>::New(writeFormat))
            {}

            //- Reset format and field list
            void clear()
            {
                DynamicList<word>::clear();
                formatter.clear();
            }

            //- Assign a new formatter
            void operator=(const word& writeFormat)
            {
                formatter = writer<Type>::New(writeFormat);
            }

        };

    // Static data members

        //- Output verbosity
        static bool verbose_;


    // Private data

        //- Const reference to fvMesh
        const fvMesh& mesh_;

        //- Load fields from files (not from objectRegistry)
        bool loadFromFiles_;

        //- Output path
        fileName outputPath_;
        
        dictionary dict_;


    // Private data

        // Read from dictonary

			//- conditional field List
			//- must give, and must be a volScalarField not vector
			wordReList conditionalFieldNames_;

			PtrList<volScalarField> conditionalFields_;
			
			//- Number of Bins
			label nBins_;
			
			//- The conditional field is cutted into nBins_ piece uniformly
			//- From min to max
			//- The conditionalFieldOutputs_ record the list of these pieces
			PtrList<scalarList> conditionalFieldOutputs_;
			
            //- Names of fields to sample
            wordReList fieldSelection_;

            //- Output format to use
            word writeFormat_;


        // Categorized scalar/vector/tensor fields

            fieldGroup<scalar> scalarFields_;
            fieldGroup<vector> vectorFields_;
            fieldGroup<sphericalTensor> sphericalTensorFields_;
            fieldGroup<symmTensor> symmTensorFields_;
            fieldGroup<tensor> tensorFields_;


    // Private Member Functions

        //- Clear old field groups
        void clearFieldGroups();

        //- Append fieldName to the appropriate group
        label appendFieldGroup(const word& fieldName, const word& fieldType);

        //- Classify field types, returns the number of fields
        label classifyFields();

        //- Combine values from all processors.
        //  Valid result only on master processor.
        template<class T>
        void combineSampledValues
        (
            PtrList<Field<T>>& sampledFields,
            PtrList<PtrList<Field<T>>>& masterFields
        );

        template<class Type>
        void writeSampleFile
        (
            const coordSet& masterSampleSet,
            const PtrList<Field<Type>>& masterFields,
			const wordList nameList,
            const fileName& timeDir,
            const writer<Type>& formatter
        );

        template<class Type>
        void sampleAndWrite(fieldGroup<Type>& fields);


        //- Disallow default bitwise copy construct and assignment
        conditionalAverage(const conditionalAverage&);
        void operator=(const conditionalAverage&);


public:

    //- Runtime type information
    TypeName("conditionalAverage");


    // Constructors

        //- Construct from Time and dictionary
        conditionalAverage
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        conditionalAverage
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~conditionalAverage();


    // Member Functions

        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the conditionalAverage
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "conditionalAverageTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
