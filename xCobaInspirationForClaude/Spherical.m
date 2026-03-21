(*==============*)
(*  Spherical  *)
(*==============*)

<<xAct`xPlain`;
Title@"Spherically symmetric field equations";

$ThisDirectory=If[NotebookDirectory[]==$Failed,
	Directory[],
	NotebookDirectory[],
	NotebookDirectory[]];

$PDFPath=FileNameJoin@{$ThisDirectory,"Spherical.pdf"};
$CellsPath=FileNameJoin@{$ThisDirectory,"SphericalCells.mx"};
$Cells={};
$ColorTitle=RGBColor[0.1,0.2,0.6];
$ColorSection=RGBColor[0.15,0.3,0.65];
$ColorSubsection=RGBColor[0.2,0.4,0.7];

AppendCell[style_,content_]:=Module[{opts={}},
	Switch[style,
		"Title",opts={FontColor->$ColorTitle},
		"Section",opts={FontColor->$ColorSection},
		"Subsection",opts={FontColor->$ColorSubsection}
	];
	$Cells=Append[$Cells,Cell[content,style,Sequence@@opts]];
];
AppendMath[expr_MatrixForm]:=Module[{boxes},
	boxes=ToBoxes@Text@expr;
	boxes=boxes/.(3/4*$Failed)->9;
	boxes=boxes/."("->""/.")"->"";
	$Cells=Append[$Cells,Cell[BoxData@boxes,"Output"]];
];
AppendMath[expr_]:=Module[{boxes},
	boxes=ToBoxes@Text@expr;
	boxes=boxes/.(3/4*$Failed)->9;
	$Cells=Append[$Cells,Cell[BoxData@boxes,"Output"]];
];
AppendText[text_]:=AppendCell["Text",text];

(*=========================================*)
(*  Part I: Lorentzian Einstein-Klein-Gordon *)
(*=========================================*)

AppendCell["Title","Part I: Lorentzian Einstein-Klein-Gordon"];
AppendText@"This section serves as a minimal test of the xAct pipeline: we derive the Einstein-Klein-Gordon field equations, project them to spherical-polar components, and verify against the Schwarzschild vacuum solution. The physics here is standard; the purpose is to confirm that the variational, chart, and component machinery is functioning correctly before proceeding to the Euclidean sector.";

Comment@"Load xAct packages and suppress verbose output.";
Get@FileNameJoin@{$ThisDirectory,"Spherical","Initialization.m"};

Comment@"Define a four-dimensional manifold and metric.";
Get@FileNameJoin@{$ThisDirectory,"Spherical","DefFourGeometry.m"};

Comment@"Load helper functions for component calculations.";
Get@FileNameJoin@{$ThisDirectory,"Spherical","SpecialFunctions.m"};

Section@"Scalar field and Klein-Gordon Lagrangian";
AppendCell["Section","Scalar field and Klein-Gordon Lagrangian"];

Comment@"Define the scalar field on the manifold.";
DefTensor[Phi[],M4,PrintAs->"\[Phi]"];

Comment@"Define the gravitational coupling constant.";
DefConstantSymbol[GravitationalCoupling,PrintAs->"\[Kappa]"];

Comment@"The kinetic scalar is X = -(1/2) g^{ab} nabla_a phi nabla_b phi.";
AppendText@"The kinetic scalar X is defined as minus one half the squared gradient of the scalar field.";
KineticX=-CD[-a]@Phi[]*CD[a]@Phi[]/2;
Print["X = ", InputForm@KineticX];
AppendMath[KineticX];

Comment@"Write the covariant Einstein-Klein-Gordon Lagrangian density (without the sqrt(-g) factor).";
AppendText@"The Lagrangian density for Einstein gravity with a free massless scalar field:";
Lagrangian=RicciScalarCD[]/(2*GravitationalCoupling)+KineticX;
Print["L = ", InputForm@Lagrangian];
AppendMath[Lagrangian];

Section@"Covariant field equations via VarD";
AppendCell["Section","Covariant field equations via VarD"];

Comment@"Vary the Lagrangian with respect to the inverse metric to obtain the Einstein equations.";
AppendText@"The metric field equation is obtained by varying the Lagrangian with respect to the inverse metric:";
MetricEOM=VarD[G[a,b],CD][Lagrangian];
MetricEOM//=ToCanonical;
MetricEOM//=ContractMetric;
MetricEOM//=ScreenDollarIndices;
MetricEOM//=CollectTensors;

Print["Metric EOM = ", InputForm@MetricEOM];
AppendMath[MetricEOM==0];

Comment@"Vary the Lagrangian with respect to the scalar field to obtain the Klein-Gordon equation.";
AppendText@"The scalar field equation is obtained by varying the Lagrangian with respect to the scalar field:";
ScalarEOM=VarD[Phi[],CD][Lagrangian];
ScalarEOM//=ToCanonical;
ScalarEOM//=ContractMetric;
ScalarEOM//=ScreenDollarIndices;
ScalarEOM//=CollectTensors;
Print["Scalar EOM = ", InputForm@ScalarEOM];
AppendMath[ScalarEOM==0];


Section@"Chart and metric components";
AppendCell["Section","Chart and metric components"];

Comment@"Set up spherical polar coordinates on the four-dimensional manifold.";
$Dimensions=4;
Get@FileNameJoin@{$ThisDirectory,"Spherical","SphericalPolarChart.m"};
$Chart=SphericalPolarChart;
$ChristoffelPDChart=ChristoffelPDSphericalPolarChart;

Comment@"Define the unknown scalar functions.";
DefScalarFunction[MetricF,PrintAs->"\[ScriptF]"];
DefScalarFunction[MetricH,PrintAs->"\[ScriptH]"];
DefScalarFunction[ScalarProfile,PrintAs->"\[Phi]"];
DefConstantSymbol[SchwarzschildMass,PrintAs->"\[ScriptCapitalM]"];

Comment@"The line element takes the general static spherically symmetric form.";
AppendText@"Lorentzian line element with signature (-,+,+,+):";
MetricComponents={
	{-MetricF[CR[]],0,0,0},
	{0,MetricH[CR[]],0,0},
	{0,0,CR[]^2,0},
	{0,0,0,CR[]^2*Sin[CTheta[]]^2}
};
AppendMath[MatrixForm@MetricComponents];

Comment@"Initialize the metric components and compute all associated quantities.";
Block[{Print=Null},
	MetricInBasis[G,-$Chart,MetricComponents];
	MetricCompute[G,$Chart,All];
];

Comment@"Set the scalar field components: Phi depends only on the radial coordinate.";
Block[{Print=Null},
	ComponentValue[ComponentArray@Phi[],ScalarProfile[CR[]]];
];

Section@"Component field equations";
AppendCell["Section","Component field equations"];

Subsection@"Metric equation components";
AppendCell["Subsection","Metric equation components"];

MetricEOMComp=MetricEOM;
MetricEOMComp//=ToBasis[$Chart];
MetricEOMComp=MetricEOMComp/.{$ChristoffelPDChart->Zero};
MetricEOMComp//=ToBasis[$Chart];
MetricEOMComp=MetricEOMComp/.{$ChristoffelPDChart->Zero};
MetricEOMComp//=Splinter;
MetricEOMComp//=Map[Simplify,#,{2}]&;

CoordLabels={CT[],CR[],CTheta[],CPhi[]};
NDim=Length@MetricEOMComp;
Do[
	comp=MetricEOMComp[[i,j]]//Simplify;
	Print["E(",CoordLabels[[i]],",",CoordLabels[[j]],") = 0 : ", InputForm@comp];
	AppendMath[Row[{"(",CoordLabels[[i]],",",CoordLabels[[j]],"): "}]];
	AppendMath[comp==0];
,{i,NDim},{j,i,NDim}];

EinsteinTT=MetricEOMComp[[1,1]]//Simplify;
EinsteinRR=MetricEOMComp[[2,2]]//Simplify;
EinsteinThTh=MetricEOMComp[[3,3]]//Simplify;

Subsection@"Scalar equation components";
AppendCell["Subsection","Scalar equation components"];

ScalarEOMComp=ScalarEOM;
ScalarEOMComp//=ToBasis[$Chart];
ScalarEOMComp=ScalarEOMComp/.{$ChristoffelPDChart->Zero};
ScalarEOMComp//=ToBasis[$Chart];
ScalarEOMComp=ScalarEOMComp/.{$ChristoffelPDChart->Zero};
ScalarEOMComp//=Splinter;
ScalarEOMComp//=Simplify;

Print["KG = 0 : ", InputForm@ScalarEOMComp];
AppendMath[ScalarEOMComp==0];

Section@"Verification: Schwarzschild solution (vacuum)";
AppendCell["Section","Verification: Schwarzschild solution (vacuum)"];

AppendText@"We substitute the Schwarzschild metric functions and set the scalar field to zero, then verify that all field equations vanish identically.";
SchwarzschildRule={
	MetricF->(1-2*SchwarzschildMass/#&),
	MetricH->(1/(1-2*SchwarzschildMass/#)&),
	ScalarProfile->(0&)
};

CheckETT=EinsteinTT/.SchwarzschildRule//FullSimplify;
CheckERR=EinsteinRR/.SchwarzschildRule//FullSimplify;
CheckEThTh=EinsteinThTh/.SchwarzschildRule//FullSimplify;
CheckKG=ScalarEOMComp/.SchwarzschildRule//FullSimplify;

Print["E_tt(Schwarzschild) = ", CheckETT];
Print["E_rr(Schwarzschild) = ", CheckERR];
Print["E_thth(Schwarzschild) = ", CheckEThTh];
Print["KG(Schwarzschild) = ", CheckKG];

If[CheckETT===0&&CheckERR===0&&CheckEThTh===0&&CheckKG===0,
	Print["PASS: All field equations vanish for Schwarzschild with trivial scalar."];
	AppendText@"All field equations vanish identically for the Schwarzschild vacuum with trivial scalar. Verified.",
	Print["FAIL: Field equations do not all vanish."];
	AppendText@"WARNING: Not all field equations vanish for Schwarzschild."
];

(*=========================================*)
(*  Hand off to fresh kernel for Part II   *)
(*=========================================*)

Section@"Handing off to fresh kernel for Euclidean computation.";
Print["Saving cells to ", $CellsPath];
Export[$CellsPath,$Cells];

Print["Launching fresh kernel for Euclidean phase..."];
$EuclideanResult=RunProcess[{"wolfram","-script",
	FileNameJoin@{$ThisDirectory,"SphericalEuclidean.m"}},
	ProcessDirectory->$ThisDirectory];
Print[$EuclideanResult["StandardOutput"]];
If[$EuclideanResult["StandardError"]=!="",
	Print["STDERR: ",$EuclideanResult["StandardError"]]];
Print["Euclidean phase complete."];

Supercomment@"(This is the end of the Lorentzian script.)";
Quit[];
