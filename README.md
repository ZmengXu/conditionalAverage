# conditionalAverage
OpenFOAM functionObject for postProcessing, getting the conditional averaged fields with a given scalarField

```
functions
{
	conditionalAverageTest
	{
		type conditionalAverage;
		libs ("libconditionalAverage.so");
		writeFormat       raw;//ensight,gnuplot,jplot,raw,vtk,xmgr
		conditionalFields 	(Z);// Must be volScalarField, can also be 'x','y','z','direction'.
		//directionVector	(1 1 1);//If conditionalFields contains 'direction', a directionVector will be needed
		nBins			20;//The number of sampled data, the conditionalFields are averaged devided into nBins, [min:(max-min)/nBins:max].
		averagedFields		( U C7H16 OH meanChemistryHRR );
	}
}
```
