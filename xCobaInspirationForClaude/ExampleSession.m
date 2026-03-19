(*========================*)
(*  ExampleSession.m      *)
(*========================*)

(* Standalone session: imports Spherical.mx and displays the equations. *)
(* No xAct required. *)

$ThisDirectory=If[NotebookDirectory[]==$Failed,
	Directory[],
	NotebookDirectory[],
	NotebookDirectory[]];

$MXPath=FileNameJoin@{$ThisDirectory,"Spherical.mx"};
$PDFPath=FileNameJoin@{$ThisDirectory,"ExampleSession.pdf"};

Print["Importing equations from ", $MXPath];
SphericalEquations=Import[$MXPath];
Print["Top-level keys: ", Keys@SphericalEquations];

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
AppendMath[expr_]:=Module[{boxes},
	boxes=ToBoxes@Text@expr;
	boxes=boxes/.(3/4*$Failed)->9;
	$Cells=Append[$Cells,Cell[BoxData@boxes,"Output"]];
];
AppendText[text_]:=AppendCell["Text",text];

AppendCell["Title","Spherical Equations (imported from Spherical.mx)"];
AppendCode[code_String]:=$Cells=Append[$Cells,Cell[code,"Text",FontFamily->"Courier",FontSize->11]];

AppendText@"This document lists the field equations exported by the Spherical computation pipeline. All equations are stored in a single Wolfram Language file, Spherical.mx, as a nested Association. No packages (such as xAct) are required to load or manipulate them.";
AppendCell["Subsection","Importing the equations"];
AppendText@"To import the equations into a fresh Wolfram Language session, evaluate:";
AppendCode["SphericalEquations = Import[\"Spherical.mx\"]"];
AppendText@"This returns an Association whose top-level keys name each equation set:";
AppendCode["Keys[SphericalEquations]"];
AppendText@"Each equation set is itself an Association keyed by coordinate component lists. For example, to access the (\[Tau],\[Tau]) component of the Euclidean metric equation of motion:";
AppendCode["SphericalEquations[\"EuclideanMetricEOM\"][{\"CTau\", \"CTau\"}]"];
AppendText@"Scalar quantities (such as the kinetic scalar or the scalar field equation) use an empty list {} as their key:";
AppendCode["SphericalEquations[\"EuclideanScalarEOM\"][{}]"];
AppendText@"The kinetic scalar X_E is stored directly (not as a sub-association):";
AppendCode["SphericalEquations[\"KineticXEScalar\"]"];
AppendText@"Each expression returned is an ordinary Wolfram Language expression involving the coordinate scalars CTau[], CRho[], CTheta[], CPhi[] and the functions ArealRadius, MetricA, ClockProfile, FuncG4, together with their derivatives. These can be manipulated with standard tools such as Simplify, Solve, DSolve, etc.";
AppendCell["Subsection","Equation listing"];

Do[
	fieldName=Keys[SphericalEquations][[f]];
	fieldVal=SphericalEquations[fieldName];
	AppendCell["Section",fieldName];
	Print["\n=== ", fieldName, " ==="];

	If[AssociationQ[fieldVal],
		Do[
			compKey=Keys[fieldVal][[c]];
			expr=fieldVal[compKey];
			accessStr=ToString[InputForm@SphericalEquations[[Key[fieldName],Key[compKey]]]];
			Print[accessStr, ":"];
			Print["  ", Short[expr,2]];
			AppendCode["SphericalEquations[\""<>fieldName<>"\"]["<>ToString[InputForm@compKey]<>"]"];
			AppendMath[expr==0];
		,{c,Length@fieldVal}];
	,
		accessStr=ToString[InputForm@SphericalEquations[[Key[fieldName]]]];
		Print[accessStr, ":"];
		Print["  ", Short[fieldVal,2]];
		AppendCode["SphericalEquations[\""<>fieldName<>"\"]"];
		AppendMath[fieldVal];
	];
,{f,Length@SphericalEquations}];

Print["\nExporting to ", $PDFPath];
UsingFrontEnd[
	nb=NotebookPut@Notebook[$Cells];
	NotebookPrint[nb,$PDFPath];
	NotebookClose[nb];
];
Print["Export complete."];

Quit[];
