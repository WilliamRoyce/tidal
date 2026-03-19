(*=======================*)
(*  SphericalPolarChart  *)
(*=======================*)

Comment@"Set up spherical polar coordinates.";
Switch[$Dimensions,
	4,
	DefChart[SphericalPolarChart,M4,{0,1,2,3},{CT[],CR[],CTheta[],CPhi[]}];
	Format@CT[]^=Symbol@"\[ScriptT]";
	Format@CR[]^=Symbol@"\[ScriptR]";
	Format@CTheta[]^=Symbol@"\[Theta]";
	Format@CPhi[]^=Symbol@"\[Phi]";
	DisplayExpression@{CT[],CR[],CTheta[],CPhi[]};,
	3,
	DefChart[SphericalPolarChart,M3,{1,2,3},{CR[],CTheta[],CPhi[]}];
	DefConstantSymbol[ConstantCT,PrintAs->"\[ScriptT]"];
	CT[]:=ConstantCT;
	Format@CR[]^=Symbol@"\[ScriptR]";
	Format@CTheta[]^=Symbol@"\[Theta]";
	Format@CPhi[]^=Symbol@"\[Phi]";
	DisplayExpression@{CR[],CTheta[],CPhi[]};
];
