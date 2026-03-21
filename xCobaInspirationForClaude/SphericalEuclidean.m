(*========================*)
(*  SphericalEuclidean   *)
(*========================*)

<<xAct`xPlain`;
Title@"Euclidean phase (fresh kernel)";

$ThisDirectory=If[NotebookDirectory[]==$Failed,
	Directory[],
	NotebookDirectory[],
	NotebookDirectory[]];

$PDFPath=FileNameJoin@{$ThisDirectory,"Spherical.pdf"};
$CellsPath=FileNameJoin@{$ThisDirectory,"SphericalCells.mx"};

Comment@"Load accumulated cells from Lorentzian phase.";
$Cells=Import[$CellsPath];
If[!ListQ[$Cells],$Cells={}];
Print["Loaded ", Length@$Cells, " cells from Lorentzian phase."];

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
(*  Part II: Euclidean scalar-gravity      *)
(*=========================================*)

AppendCell["Title","Part II: Euclidean Riemannian gravity (Mukohyama-Uzan)"];
AppendText@"This section contains the core scientific computation. We work in Euclidean (Riemannian) signature following Mukohyama-Uzan (1301.1361), deriving the field equations for the Horndeski G4 sector with a Lemaitre line element (cf. Eq. 4.2 of 2208.02943). The metric and scalar field equations are projected to spherical-polar components and exported as a standalone equation set for further analysis.";

Comment@"Load xAct packages in a fresh kernel.";
Off@ValidateSymbol::used;
<<xAct`xTensor`;
<<xAct`xTras`;
<<xAct`xCoba`;
On@ValidateSymbol::used;
Off@PrintAsCharacter::argx;
Off@CollectTensors::denominator;
$DefInfoQ=False;
Unprotect@AutomaticRules;
Options[AutomaticRules]={Verbose->False};
Protect@AutomaticRules;

Section@"Euclidean manifold and metric";
AppendCell["Section","Euclidean manifold and metric"];

Comment@"Define a four-dimensional manifold with a positive-definite (Euclidean) metric.";
AppendText@"We define a Riemannian manifold with signature (+,+,+,+).";
DefManifold[M4E,4,IndexRange[a,z]];
DefMetric[1,GE[-a,-b],CDE,
	{";","\*OverscriptBox[\(\[ScriptCapitalD]\),\(_\)]"},
	PrintAs->"\*OverscriptBox[\(\[ScriptG]\),\(\[ScriptCapitalE]\)]"];

Comment@"Define the clock field on the Euclidean manifold.";
DefTensor[ClockField[],M4E,PrintAs->"\[Phi]"];

Comment@"Define the scalar functions G4(X_E) and K(X_E) which parameterise the gravitational action.";
DefScalarFunction[FuncG4,PrintAs->"\!\(\*SubscriptBox[\(\[ScriptCapitalG]\), \(4\)]\)"];
DefScalarFunction[FuncK,PrintAs->"\[ScriptCapitalK]"];

Comment@"Define X_E as an abstract scalar for cleaner output of the final equations.";
DefTensor[KineticXEScalar[],M4E,PrintAs->"\!\(\*SubscriptBox[\(X\), \(\[ScriptCapitalE]\)]\)"];

Section@"Riemannian gravity Lagrangian";
AppendCell["Section","Riemannian gravity Lagrangian"];

AppendText@"The Lagrangian combines the Horndeski G4 sector (Eq. 3.6 of 1301.1361) with the quadratic-curvature invariants of the ELST theory (Eq. 4 of 2310.17266).";

Comment@"The kinetic scalar of the clock field is X_E = g_E^{ab} nabla_a phi nabla_b phi.";
AppendText@"The kinetic scalar of the clock field:";
KineticXE=CDE[-a]@ClockField[]*CDE[a]@ClockField[];
Print["X_E = ", InputForm@KineticXE];
AppendMath[KineticXE];

Comment@"Define the constant couplings for the quadratic-curvature sector (cf. Eq. 4 of 2310.17266).";
DefConstantSymbol[CouplingC1,PrintAs->"\!\(\*SubscriptBox[\(c\), \(1\)]\)"];
DefConstantSymbol[CouplingC2,PrintAs->"\!\(\*SubscriptBox[\(c\), \(2\)]\)"];
DefConstantSymbol[CouplingC3,PrintAs->"\!\(\*SubscriptBox[\(c\), \(3\)]\)"];

Comment@"Construct the Lagrangian. Each contraction uses distinct dummy indices.";
BoxPhiBoxPhi=CDE[-c]@CDE[c]@ClockField[]*CDE[-d]@CDE[d]@ClockField[];
NablaNablaPhi2=CDE[-c]@CDE[-d]@ClockField[]*CDE[c]@CDE[d]@ClockField[];

Comment@"Horndeski G4 sector (Eq. 3.6 of 1301.1361).";
EinsteinHilbertTerm=FuncG4[KineticXE]*RicciScalarCDE[];
KTerm=FuncK[KineticXE];
GalileonTerm=-2*Derivative[1][FuncG4][KineticXE]*(BoxPhiBoxPhi-NablaNablaPhi2);

Comment@"Quadratic-curvature invariants (L4 of Eq. 4 of 2310.17266).";
RicciScalarSquared=RicciScalarCDE[]^2;
RicciTensorSquared=RicciCDE[-a,-b]*RicciCDE[a,b];
KretschnerScalar=RiemannCDE[-a,-b,-c,-d]*RiemannCDE[a,b,c,d];

AppendText@"The Horndeski G4 sector consists of the Ricci scalar coupled to G4, the shift-symmetric potential K, and the Galileon combination involving second covariant derivatives of the clock field.";
AppendText@"The quadratic-curvature sector adds three invariants: the square of the Ricci scalar, the contraction of the Ricci tensor with itself, and the Kretschner scalar (the full contraction of the Riemann tensor), with constant couplings c1, c2, c3.";

EuclideanLagrangian=(EinsteinHilbertTerm+KTerm+GalileonTerm
	+CouplingC1*RicciScalarSquared
	+CouplingC2*RicciTensorSquared
	+CouplingC3*KretschnerScalar);
Print["L_E = ", InputForm@EuclideanLagrangian];
AppendMath[EuclideanLagrangian];

Section@"Euclidean field equations via VarD";
AppendCell["Section","Euclidean field equations via VarD"];

Comment@"Vary the Lagrangian with respect to the Euclidean inverse metric to obtain the gravitational field equation.";
AppendText@"Metric field equation from varying the gravitational action with respect to the inverse Euclidean metric:";
EuclideanMetricEOM=VarD[GE[a,b],CDE][EuclideanLagrangian];
EuclideanMetricEOM//=ToCanonical;
EuclideanMetricEOM//=ContractMetric;
EuclideanMetricEOM//=ScreenDollarIndices;
EuclideanMetricEOM//=CollectTensors;

Comment@"Replace the expanded kinetic scalar with the abstract X_E in all function arguments.";
CleanXE[expr_]:=expr/.{Scalar[CDE[-c_]@ClockField[]*CDE[c_]@ClockField[]]->KineticXEScalar[],
	CDE[-c_]@ClockField[]*CDE[c_]@ClockField[]->KineticXEScalar[]};
EuclideanMetricEOM//=CleanXE;
Print["Euclidean Metric EOM (raw) = ", InputForm@EuclideanMetricEOM];
AppendText@"Metric field equation (before introducing shorthand tensors):";
AppendMath[EuclideanMetricEOM==0];

Comment@"Vary the Lagrangian with respect to the clock field to obtain the scalar field equation.";
AppendText@"Scalar field equation from varying the gravitational action with respect to the clock field:";
Print["Computing scalar EOM via VarD..."];
EuclideanScalarEOM=VarD[ClockField[],CDE][EuclideanLagrangian];
EuclideanScalarEOM//=ScreenDollarIndices;
EuclideanScalarEOM//=ToCanonical;
EuclideanScalarEOM//=ContractMetric;
EuclideanScalarEOM//=ScreenDollarIndices;
EuclideanScalarEOM//=CollectTensors;
EuclideanScalarEOM//=CleanXE;
Print["Euclidean Scalar EOM = ", InputForm@EuclideanScalarEOM];
AppendText@"Scalar field equation (covariant):";
AppendMath[EuclideanScalarEOM==0];

Comment@"Define shorthand tensors for covariant derivatives of the clock field.";
DefTensor[DClockField[-a],M4E,PrintAs->"\[Del]\[Phi]"];
DefTensor[DDClockField[-a,-b],M4E,PrintAs->"\[Del]\[Del]\[Phi]"];
DefTensor[DDDClockField[-a,-b,-c],M4E,PrintAs->"\[Del]\[Del]\[Del]\[Phi]"];
DefTensor[DDDDClockField[-a,-b,-c,-d],M4E,PrintAs->"\[Del]\[Del]\[Del]\[Del]\[Phi]"];

ToDDDDClockField=MakeRule[{CDE[-a]@CDE[-b]@CDE[-c]@CDE[-d]@ClockField[],DDDDClockField[-a,-b,-c,-d]},
	MetricOn->All,ContractMetrics->True];
ToDDDClockField=MakeRule[{CDE[-a]@CDE[-b]@CDE[-c]@ClockField[],DDDClockField[-a,-b,-c]},
	MetricOn->All,ContractMetrics->True];
ToDDClockField=MakeRule[{CDE[-a]@CDE[-b]@ClockField[],DDClockField[-a,-b]},
	MetricOn->All,ContractMetrics->True];
ToDClockField=MakeRule[{CDE[-a]@ClockField[],DClockField[-a]},
	MetricOn->All,ContractMetrics->True];

ApplyShorthand[expr_]:=Module[{e=expr},
	e=e/.ToDDDDClockField;
	e=e/.ToDDDClockField;
	e=e/.ToDDClockField;
	e=e/.ToDClockField;
	e//=ToCanonical;
	e//=ContractMetric;
	e//=ScreenDollarIndices;
	e//=CollectTensors;
e];

Comment@"Substitute shorthand tensors into the metric field equation (rank 4 first, then 3, 2, 1).";
EuclideanMetricEOM//=ApplyShorthand;
Print["Euclidean Metric EOM (reduced) = ", InputForm@EuclideanMetricEOM];
AppendText@"Metric field equation with shorthand tensors:";
AppendMath[EuclideanMetricEOM==0];

Comment@"Substitute shorthand tensors into the scalar field equation.";
EuclideanScalarEOM//=ApplyShorthand;
Print["Euclidean Scalar EOM (reduced) = ", InputForm@EuclideanScalarEOM];
AppendText@"Scalar field equation with shorthand tensors:";
AppendMath[EuclideanScalarEOM==0];

(*=========================================*)
(*  Pre-compute shorthand tensor components *)
(*=========================================*)

Section@"Euclidean chart and pre-computed components";
AppendCell["Section","Euclidean chart and pre-computed components"];

Comment@"Set up Lemaitre coordinates on the Euclidean manifold.";
DefChart[EuclideanChart,M4E,{1,2,3,4},{CTau[],CRho[],CTheta[],CPhi[]},
	ChartColor->RGBColor[0,0,1]];
Format@CTau[]^=Symbol@"\[Tau]";
Format@CRho[]^=Symbol@"\[Rho]";
Format@CTheta[]^=Symbol@"\[Theta]";
Format@CPhi[]^=Symbol@"\[Phi]";

Comment@"Define the unknown metric functions for the Euclidean Lemaitre line element.";
DefScalarFunction[ArealRadius,PrintAs->"\[ScriptR]"];
DefScalarFunction[MetricA,PrintAs->"\[ScriptCapitalA]"];
DefScalarFunction[ClockProfile,PrintAs->"\[Phi]"];

Comment@"The Euclidean Lemaitre line element with signature (+,+,+,+), cf. Eq. (4.2) of 2208.02943.";
AppendText@"Euclidean Lemaitre line element with signature (+,+,+,+):";
EuclideanMetricComponents={
	{1,0,0,0},
	{0,1-MetricA[CTau[],CRho[]],0,0},
	{0,0,ArealRadius[CTau[],CRho[]]^2,0},
	{0,0,0,ArealRadius[CTau[],CRho[]]^2*Sin[CTheta[]]^2}
};
AppendMath[MatrixForm@EuclideanMetricComponents];

Comment@"Initialize the Euclidean metric components.";
Block[{Print=Null},
	MetricInBasis[GE,-EuclideanChart,EuclideanMetricComponents];
	MetricCompute[GE,EuclideanChart,All];
];

Comment@"Set the clock field components: depends on tau and rho.";
Block[{Print=Null},
	ComponentValue[ComponentArray@ClockField[],ClockProfile[CTau[],CRho[]]];
];

Comment@"Helper to compute components of an expression via the full projection pipeline.";
EuclideanSplinter[expr_]:=Module[{e=expr},
	e//=ToBasis[EuclideanChart];
	e=e/.{ChristoffelPDEuclideanChart->Zero};
	e//=ToBasis[EuclideanChart];
	e=e/.{ChristoffelPDEuclideanChart->Zero};
	e//=SeparateMetric[GE];
	e//=ToBasis[EuclideanChart];
	e//=ComponentArray;
	e//=ToValues;
	e//=TraceBasisDummy;
	e//=ToValues;
e];

Comment@"Helper to verify component assignment by comparing projection with reference.";
VerifyComponents[name_,tensor_,reference_]:=Module[{check,diff},
	check=tensor;
	check//=ToBasis[EuclideanChart];
	check//=ComponentArray;
	check//=ToValues;
	diff=Flatten[check-reference]//FullSimplify//DeleteDuplicates;
	If[diff==={0},
		Print["PASS: ", name],
		Print["FAIL: ", name, " (", Length@DeleteCases[diff,0], " nonzero entries)"]
	];
];

Comment@"Compute downstairs components of CDE[-a]@ClockField[] via the standard pipeline.";
DPhiDown=EuclideanSplinter[CDE[-a]@ClockField[]];

Comment@"Assign these as the component values of the shorthand DClockField vector.";
Block[{Print=Null},
	ComponentValue[ComponentArray@DClockField[-{a,EuclideanChart}],DPhiDown];
];

Comment@"Compute upstairs components by raising with the inverse metric.";
DPhiUp=EuclideanSplinter[DClockField[a]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DClockField[{a,EuclideanChart}],DPhiUp];
];

Comment@"Verify: project DClockField[-a] to components and check agreement.";
DPhiCheck=DClockField[-a];
DPhiCheck//=ToBasis[EuclideanChart];
DPhiCheck//=ComponentArray;
DPhiCheck//=ToValues;
Print["DClockField downstairs check = ", InputForm@DPhiCheck];

DPhiDiff=Flatten[DPhiCheck-DPhiDown]//FullSimplify//DeleteDuplicates;
If[DPhiDiff==={0},
	Print["PASS: DClockField downstairs components agree."],
	Print["FAIL: DClockField downstairs components disagree."];
	Print["  Nonzero entries: ", InputForm@DeleteCases[DPhiDiff,0]];
];

Subsection@"DDClockField (rank 2)";

Comment@"Compute DDClockField as the covariant derivative of the pre-computed DClockField vector.";
DDPhiDown=EuclideanSplinter[CDE[-a]@DClockField[-b]];
DDPhiDown=Map[Simplify,DDPhiDown,{2}];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDClockField[-{a,EuclideanChart},-{b,EuclideanChart}],DDPhiDown];
];

DDPhiUp=EuclideanSplinter[DDClockField[a,b]];
DDPhiUp=Map[Simplify,DDPhiUp,{2}];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDClockField[{a,EuclideanChart},{b,EuclideanChart}],DDPhiUp];
];

DDPhiMixedDown1=EuclideanSplinter[DDClockField[-a,b]];
DDPhiMixedDown1=Map[Simplify,DDPhiMixedDown1,{2}];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDClockField[-{a,EuclideanChart},{b,EuclideanChart}],DDPhiMixedDown1];
];

DDPhiMixedDown2=EuclideanSplinter[DDClockField[a,-b]];
DDPhiMixedDown2=Map[Simplify,DDPhiMixedDown2,{2}];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDClockField[{a,EuclideanChart},-{b,EuclideanChart}],DDPhiMixedDown2];
];

VerifyComponents["DDClockField[-a,-b]",DDClockField[-a,-b],DDPhiDown];
VerifyComponents["DDClockField[a,b]",DDClockField[a,b],DDPhiUp];
VerifyComponents["DDClockField[-a,b]",DDClockField[-a,b],DDPhiMixedDown1];
VerifyComponents["DDClockField[a,-b]",DDClockField[a,-b],DDPhiMixedDown2];

Subsection@"DDDClockField (rank 3)";

Comment@"Compute DDDClockField as the covariant derivative of the pre-computed DDClockField tensor.";
DDDPhiDown=EuclideanSplinter[CDE[-a]@DDClockField[-b,-c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[-{a,EuclideanChart},-{b,EuclideanChart},-{c,EuclideanChart}],DDDPhiDown];
];

DDDPhiUp=EuclideanSplinter[DDDClockField[a,b,c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[{a,EuclideanChart},{b,EuclideanChart},{c,EuclideanChart}],DDDPhiUp];
];

DDDPhiMixed100=EuclideanSplinter[DDDClockField[a,-b,-c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[{a,EuclideanChart},-{b,EuclideanChart},-{c,EuclideanChart}],DDDPhiMixed100];
];

DDDPhiMixed010=EuclideanSplinter[DDDClockField[-a,b,-c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[-{a,EuclideanChart},{b,EuclideanChart},-{c,EuclideanChart}],DDDPhiMixed010];
];

DDDPhiMixed001=EuclideanSplinter[DDDClockField[-a,-b,c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[-{a,EuclideanChart},-{b,EuclideanChart},{c,EuclideanChart}],DDDPhiMixed001];
];

DDDPhiMixed110=EuclideanSplinter[DDDClockField[a,b,-c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[{a,EuclideanChart},{b,EuclideanChart},-{c,EuclideanChart}],DDDPhiMixed110];
];

DDDPhiMixed101=EuclideanSplinter[DDDClockField[a,-b,c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[{a,EuclideanChart},-{b,EuclideanChart},{c,EuclideanChart}],DDDPhiMixed101];
];

DDDPhiMixed011=EuclideanSplinter[DDDClockField[-a,b,c]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDClockField[-{a,EuclideanChart},{b,EuclideanChart},{c,EuclideanChart}],DDDPhiMixed011];
];

VerifyComponents["DDDClockField[-a,-b,-c]",DDDClockField[-a,-b,-c],DDDPhiDown];
VerifyComponents["DDDClockField[a,b,c]",DDDClockField[a,b,c],DDDPhiUp];
VerifyComponents["DDDClockField[a,-b,-c]",DDDClockField[a,-b,-c],DDDPhiMixed100];
VerifyComponents["DDDClockField[-a,b,-c]",DDDClockField[-a,b,-c],DDDPhiMixed010];
VerifyComponents["DDDClockField[-a,-b,c]",DDDClockField[-a,-b,c],DDDPhiMixed001];
VerifyComponents["DDDClockField[a,b,-c]",DDDClockField[a,b,-c],DDDPhiMixed110];
VerifyComponents["DDDClockField[a,-b,c]",DDDClockField[a,-b,c],DDDPhiMixed101];
VerifyComponents["DDDClockField[-a,b,c]",DDDClockField[-a,b,c],DDDPhiMixed011];

Subsection@"DDDDClockField (rank 4, only needed placements)";

Comment@"Only two index placements of DDDDClockField appear in the scalar EOM: (1,0,1,0) and (1,1,0,0).";
Comment@"Compute all-downstairs first, then raise to the needed placements.";
DDDDPhiDown=EuclideanSplinter[CDE[-a]@DDDClockField[-b,-c,-d]];
Block[{Print=Null},
	ComponentValue[ComponentArray@DDDDClockField[-{a,EuclideanChart},-{b,EuclideanChart},-{c,EuclideanChart},-{d,EuclideanChart}],DDDDPhiDown];
];

DDDDPhiMixed1010=EuclideanSplinter[DDDDClockField[a,-b,c,-d]];
Block[{Print=Null},ComponentValue[ComponentArray@DDDDClockField[{a,EuclideanChart},-{b,EuclideanChart},{c,EuclideanChart},-{d,EuclideanChart}],DDDDPhiMixed1010];];

DDDDPhiMixed1100=EuclideanSplinter[DDDDClockField[a,b,-c,-d]];
Block[{Print=Null},ComponentValue[ComponentArray@DDDDClockField[{a,EuclideanChart},{b,EuclideanChart},-{c,EuclideanChart},-{d,EuclideanChart}],DDDDPhiMixed1100];];

VerifyComponents["DDDDClockField[-a,-b,-c,-d]",DDDDClockField[-a,-b,-c,-d],DDDDPhiDown];
VerifyComponents["DDDDClockField[a,-b,c,-d]",DDDDClockField[a,-b,c,-d],DDDDPhiMixed1010];
VerifyComponents["DDDDClockField[a,b,-c,-d]",DDDDClockField[a,b,-c,-d],DDDDPhiMixed1100];

(*=============================================*)
(*  Term-by-term component projection of EOM  *)
(*=============================================*)

(*Comment@"Early exit for fast iteration on PDF formatting.";
Print["Exporting to ", $PDFPath];
UsingFrontEnd[
	nb=NotebookPut@Notebook[$Cells];
	NotebookPrint[nb,$PDFPath];
	NotebookClose[nb];
];
Print["Export complete (fast mode)."];
DeleteFile[$CellsPath];
Supercomment@"(Early exit before component projection.)";
Quit[];*)

Section@"Component field equations (term-by-term)";
AppendCell["Section","Component field equations (term-by-term)"];

Comment@"Split the covariant metric EOM into additive terms and project each individually.";
EOMTerms=List@@EuclideanMetricEOM;
NTerms=Length@EOMTerms;
Print["Euclidean metric EOM has ", NTerms, " terms."];

EuclideanMetricEOMComp=ConstantArray[0,{4,4}];
Do[
	Print["Term ", k, "/", NTerms, ": ", Short[EOMTerms[[k]],1]];
	$TermStart=AbsoluteTime[];
	termComp=EuclideanSplinter[EOMTerms[[k]]];
	$TermEnd=AbsoluteTime[];
	Print["  -> done in ", NumberForm[$TermEnd-$TermStart,{4,1}], " s"];
	EuclideanMetricEOMComp+=termComp;
,{k,NTerms}];

Comment@"Simplify the summed component equations.";
Print["Simplifying final component equations..."];
EuclideanMetricEOMComp=Map[Simplify,EuclideanMetricEOMComp,{2}];
Print["Simplification complete."];

Comment@"Display the component field equations.";
AppendText@"Euclidean metric field equations (component form):";
CoordLabelsE={CTau[],CRho[],CTheta[],CPhi[]};
Do[
	comp=EuclideanMetricEOMComp[[i,j]];
	If[comp=!=0,
		Print["E_E(",CoordLabelsE[[i]],",",CoordLabelsE[[j]],") = ", InputForm@comp];
		AppendText@StringJoin["(",ToString@CoordLabelsE[[i]],",",ToString@CoordLabelsE[[j]],"):"];
		AppendMath[comp==0];
	];
,{i,4},{j,i,4}];

Subsection@"Scalar field equation components";
AppendCell["Subsection","Scalar field equation components"];

Comment@"Split the covariant scalar EOM into additive terms and project each individually.";
ScalarEOMTerms=List@@EuclideanScalarEOM;
NScalarTerms=Length@ScalarEOMTerms;
Print["Euclidean scalar EOM has ", NScalarTerms, " terms."];

EuclideanScalarEOMComp=0;
Do[
	Print["Scalar term ", k, "/", NScalarTerms, ": ", Short[ScalarEOMTerms[[k]],1]];
	$TermStart=AbsoluteTime[];
	termComp=EuclideanSplinter[ScalarEOMTerms[[k]]];
	$TermEnd=AbsoluteTime[];
	Print["  -> done in ", NumberForm[$TermEnd-$TermStart,{4,1}], " s"];
	EuclideanScalarEOMComp+=termComp;
,{k,NScalarTerms}];

Comment@"Simplify the summed scalar component equation.";
Print["Simplifying scalar component equation..."];
EuclideanScalarEOMComp=Map[Simplify,EuclideanScalarEOMComp];
Print["Scalar simplification complete."];

Comment@"Display the scalar component equation.";
AppendText@"Euclidean scalar field equation (component form):";
Print["KG_E = ", InputForm@EuclideanScalarEOMComp];
AppendMath[EuclideanScalarEOMComp==0];

Subsection@"Kinetic scalar X_E";
AppendCell["Subsection","Kinetic scalar X_E"];

Comment@"Compute the component-level value of the kinetic scalar X_E = g_E^{ab} nabla_a phi nabla_b phi.";
Print["Computing kinetic scalar X_E components..."];
KineticXEComp=EuclideanSplinter[KineticXE];
KineticXEComp//=Simplify;
Print["X_E = ", InputForm@KineticXEComp];
AppendText@"Kinetic scalar (component form):";
AppendMath[KineticXEComp];

(*=================================*)
(*  Build equation association     *)
(*=================================*)

Section@"Equation export";
AppendCell["Section","Equation export"];

Comment@"Build association for the Euclidean metric EOM component equations.";
$CoordNames={"CTau","CRho","CTheta","CPhi"};
EuclideanMetricEOMAssoc=Association[];
Do[
	comp=EuclideanMetricEOMComp[[i,j]];
	If[comp=!=0,
		EuclideanMetricEOMAssoc[{$CoordNames[[i]],$CoordNames[[j]]}]=comp;
	];
,{i,4},{j,i,4}];

Comment@"Build association for the Euclidean scalar EOM.";
EuclideanScalarEOMAssoc=Association[{}->EuclideanScalarEOMComp];

Comment@"Assemble the top-level association.";
SphericalEquations=Association[
	"EuclideanMetricEOM"->EuclideanMetricEOMAssoc,
	"EuclideanScalarEOM"->EuclideanScalarEOMAssoc,
	"KineticXEScalar"->KineticXEComp
];

Comment@"Export the equation association to Spherical.mx.";
$MXPath=FileNameJoin@{$ThisDirectory,"Spherical.mx"};
Print["Exporting equations to ", $MXPath];
Export[$MXPath,SphericalEquations];
Print["Equation export complete. Keys: ", Keys@SphericalEquations];

(*=================================*)
(*  Export final PDF               *)
(*=================================*)

Section@"Export to PDF";
Comment@"Exporting all results to a single PDF file.";
Print["Exporting to ", $PDFPath];
UsingFrontEnd[
	nb=NotebookPut@Notebook[$Cells];
	NotebookPrint[nb,$PDFPath];
	NotebookClose[nb];
];
Print["Export complete."];

Comment@"Clean up temporary cells file.";
DeleteFile[$CellsPath];

Supercomment@"(This is the end of the Euclidean script.)";
Quit[];
