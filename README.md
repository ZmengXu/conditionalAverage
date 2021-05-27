# conditionalAverage
OpenFOAM functionObject for postProcessing, getting the conditional averaged fields with a given scalarField

It has been tested in OF4, OF6, OF7, OFv1712 and OFv1806
```
functions
{
	conditionalAverageTest
	{
		type conditionalAverage;
		libs ("libconditionalAverage.so");
        regionType all;//;cellZone;
        //name ignitionCells;//give a specific region to postProcess
		writeFormat       raw;//ensight,gnuplot,jplot,raw,vtk,xmgr
		conditionalFields 	(Z);// Must be volScalarField, can also be 'x','y','z','direction'.
		weightedAveragedField   rhoMeshV;//rhoMeshV;none;meshV
		//directionVector	(1 1 1);//If conditionalFields contains 'direction', a directionVector will be needed
		nBins			20;//The number of sampled data, the conditionalFields are averaged devided into nBins, [min:(max-min)/nBins:max].
		//maxValue               1;//The first column of output file is within [minValue, maxValue],divided into nBins number of equal parts 
        //minValue               0;//If the maxValue and the minValue is not presented, the max and min of the conditionalFields will be used instead
		averagedFields		( U C7H16 OH meanChemistryHRR );
	}
}
```
