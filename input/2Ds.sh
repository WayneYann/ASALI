rm -f 2Ds.xml
touch 2Ds.xml
echo '<?xml version="1.0" encoding="utf-8"?>' >> 2Ds.xml
echo '<models>' >> 2Ds.xml
echo '	<honeyComb>false</honeyComb>' >> 2Ds.xml
echo '	<microPackedbed>true</microPackedbed>' >> 2Ds.xml
echo '	<packedBed>false</packedBed>' >> 2Ds.xml
echo '</models>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<reaction>O-xylene-to-phthalic</reaction>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<mole>' >> 2Ds.xml
echo '	<O2>0.207</O2>' >> 2Ds.xml
echo '	<XYLENE>0.0147</XYLENE>' >> 2Ds.xml
echo '	<N2>0.7783</N2>' >> 2Ds.xml
echo '	<PHTHALIC>0</PHTHALIC>' >> 2Ds.xml
echo '	<CO>0.</CO>' >> 2Ds.xml
echo '	<CO2>0.</CO2>' >> 2Ds.xml
echo '	<H2>0.</H2>' >> 2Ds.xml
echo '	<H2O>0.</H2O>' >> 2Ds.xml
echo '</mole>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<operatingConditions>' >> 2Ds.xml
echo '	<pressure>100000</pressure>' >> 2Ds.xml
echo '	<temperauture>' >> 2Ds.xml
echo '		<coolant>625.15</coolant>' >> 2Ds.xml
echo '		<feed>625.15</feed>' >> 2Ds.xml
echo '	</temperauture>' >> 2Ds.xml
echo '	<specificMassFlowRate>0.867904008</specificMassFlowRate>' >> 2Ds.xml
echo '	<inertLength>0.3</inertLength>' >> 2Ds.xml
echo '</operatingConditions>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<packedBed>' >> 2Ds.xml
echo '	<tubeDiameter>0.0254</tubeDiameter>' >> 2Ds.xml
echo '	<particleDiameter>0.004</particleDiameter>' >> 2Ds.xml
echo '	<length>3</length>' >> 2Ds.xml
echo '</packedBed>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<honeyComb>' >> 2Ds.xml
echo '	<type>washcoated</type>' >> 2Ds.xml
echo '	<tubeDiameter>0.0254</tubeDiameter>' >> 2Ds.xml
echo '	<length>36</length>' >> 2Ds.xml
echo '	<CPSI>400</CPSI>' >> 2Ds.xml
echo '	<washcoatThickness>0.000018</washcoatThickness>' >> 2Ds.xml
echo '	<wallThickness>7</wallThickness>' >> 2Ds.xml
echo '</honeyComb>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<microBed>' >> 2Ds.xml
echo '	<tubeDiameter>0.0380456397</tubeDiameter>' >> 2Ds.xml
echo '	<particleDiameter>0.0006</particleDiameter>' >> 2Ds.xml
echo '	<length>3.0</length>' >> 2Ds.xml
echo '	<CPSI>25</CPSI>' >> 2Ds.xml
echo '	<wallThickness>25</wallThickness>' >> 2Ds.xml
echo '</microBed>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<gasProperties>' >> 2Ds.xml
echo '	<specificHeat>992</specificHeat>' >> 2Ds.xml
echo '	<viscosity>4.3e-05</viscosity>' >> 2Ds.xml
echo '	<diffusivity>1.5e-05</diffusivity>' >> 2Ds.xml
echo '	<conductivity>0.058</conductivity>' >> 2Ds.xml
echo '</gasProperties>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<catalystProperties>' >> 2Ds.xml
echo '	<conductivity>0.22</conductivity>' >> 2Ds.xml
echo '	<specificHeat>925</specificHeat>' >> 2Ds.xml
echo '	<density>2100</density>' >> 2Ds.xml
echo '	<voidFraction>0.68</voidFraction>' >> 2Ds.xml
echo '	<tortuosity>8</tortuosity>' >> 2Ds.xml
echo '</catalystProperties>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<supportProperties>' >> 2Ds.xml
echo '	<conductivity>400</conductivity>' >> 2Ds.xml
echo '	<specificHeat>925</specificHeat>' >> 2Ds.xml
echo '	<density>2100</density>' >> 2Ds.xml
echo '</supportProperties>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<grid>' >> 2Ds.xml
echo '	<axial>20</axial>' >> 2Ds.xml
echo '	<radial>30</radial>' >> 2Ds.xml
echo '</grid>' >> 2Ds.xml
echo ' ' >> 2Ds.xml
echo '<solver>' >> 2Ds.xml
echo '	<ODE>BzzMath</ODE>' >> 2Ds.xml
echo '	<BVP>BzzMath</BVP>' >> 2Ds.xml
echo '</solver>' >> 2Ds.xml
