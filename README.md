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
		writeFormat       raw;//ensight,gnuplot,jplot,raw,vtk,xmgr
		conditionalFields 	(Z);// Must be volScalarField, can also be 'x','y','z','direction'.
		weightedAveragedField   rhoMeshV;//rhoMeshV;none;meshV
		//directionVector	(1 1 1);//If conditionalFields contains 'direction', a directionVector will be needed
		nBins			20;//The number of sampled data, the conditionalFields are averaged devided into nBins, [min:(max-min)/nBins:max].
		maxF               1;
                minF               0;
		averagedFields		( U C7H16 OH meanChemistryHRR );
	}
}
```
