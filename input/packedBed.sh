rm -f packedBed.xml
touch packedBed.xml
echo '<?xml version="1.0" encoding="utf-8"?>' >> packedBed.xml
echo '<models>' >> packedBed.xml
echo '	<heterogeneous>true</heterogeneous>' >> packedBed.xml
echo '	<pseudoHomogeneous>true</pseudoHomogeneous>' >> packedBed.xml
echo '</models>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<reaction>O-xylene-to-phthalic-complex</reaction>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<mole>' >> packedBed.xml
echo '	<O2>0.2081428571</O2>' >> packedBed.xml
echo '	<XYLENE>0.00924</XYLENE>' >> packedBed.xml
echo '	<N2>0.7826171429</N2>' >> packedBed.xml
echo '</mole>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<operatingConditions>' >> packedBed.xml
echo '	<pressure>100000</pressure>' >> packedBed.xml
echo '	<temperauture>' >> packedBed.xml
echo '		<coolant>625</coolant>' >> packedBed.xml
echo '		<feed>625</feed>' >> packedBed.xml
echo '	</temperauture>' >> packedBed.xml
echo '	<specificMassFlowRate>1.3</specificMassFlowRate>' >> packedBed.xml
echo '	<externalExchange>true</externalExchange>' >> packedBed.xml
echo '</operatingConditions>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<reactorDimensions>' >> packedBed.xml
echo '	<tubeDiameter>0.0254</tubeDiameter>' >> packedBed.xml
echo '	<particleDiameter>0.004</particleDiameter>' >> packedBed.xml
echo '	<length>3</length>' >> packedBed.xml
echo '</reactorDimensions>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<gasProperties>' >> packedBed.xml
echo '	<specificHeat>992</specificHeat>' >> packedBed.xml
echo '	<viscosity>4.3e-05</viscosity>' >> packedBed.xml
echo '	<diffusivity>7e-06</diffusivity>' >> packedBed.xml
echo '	<conductivity>0.058</conductivity>' >> packedBed.xml
echo '</gasProperties>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<catalystProperties>' >> packedBed.xml
echo '	<conductivity>1.0</conductivity>' >> packedBed.xml
echo '	<specificHeat>925</specificHeat>' >> packedBed.xml
echo '	<density>2100</density>' >> packedBed.xml
echo '	<voidFraction>0.68</voidFraction>' >> packedBed.xml
echo '	<tortuosity>8</tortuosity>' >> packedBed.xml
echo '</catalystProperties>' >> packedBed.xml
echo ' ' >> packedBed.xml
echo '<solver>OpenSMOKE</solver>' >> packedBed.xml
