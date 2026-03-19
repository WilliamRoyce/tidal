(*  SpecialFunctions  *)
(*====================*)

DisplayRule~SetAttributes~HoldAll;
DisplayRule[InputExpr_,InputRule_]:=Module[{Expr=Evaluate@InputExpr,EqnLabelValue=ToString@Defer@InputRule},

	EqnLabelValue//=StringDelete[#,"Defer["]&;
	EqnLabelValue//=StringDelete[#,"]"]&;

	Expr=Expr/.Evaluate@InputRule;
	Expr//=ToCanonical;
	Expr//=ContractMetric;
	Expr//=ScreenDollarIndices;
	Expr//=CollectTensors;
	Expr//=ScreenDollarIndices;
	Expr//=FullSimplify;
	DisplayExpression[(InputExpr->Expr),EqnLabel->EqnLabelValue];
Expr];

Recanonicalize[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ToCanonical;
	Expr//=ContractMetric;
	Expr//=ScreenDollarIndices;
	Expr//=CollectTensors;
	Expr//=ScreenDollarIndices;
	Expr
];

SoftDemonstrateComponents[InputLeft_,InputRight_]:=Module[{},
	(InputLeft->MatrixForm@InputRight)//DisplayExpression;
];

ComponentAnalysis::ExtrinsicCurvatureNotEqual="There is a fatal structural mismatch (not just an overall sign difference) between your four-dimensionally-derived derived and three-dimensionally-derived extrinsic curvature components.";
ComponentAnalysis::ValuesNotEqual="Gauss-Codazzi and induced values disagree. You should expect complete agreement only for the downstairs components of the extrinsic curvature, due to abuse of the four-dimensional metric inverse, so this is not necessarily a problem. Even in the case of the extrinsic curvature with downstairs components, note that the unit-timelike normal in four dimensions is only defined up to a sign, and so the overall sign of the extrinsic curvature can differ.";

CancelValues[InputLeft_,InputRight_]:=Module[{VanishingExpr=InputLeft-InputRight},
	VanishingExpr//=List;	
	VanishingExpr//=FullSimplify;
	VanishingExpr//=Flatten;
	VanishingExpr//=DeleteDuplicates;
	If[!(VanishingExpr==={0}),
		Return@Message[ComponentAnalysis::ValuesNotEqual];
	];
];

DemonstrateComponents[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=ToBasis[$Chart];
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=ToValues;
	Expr//=MatrixForm;	
	(InputExpr->Expr)//DisplayExpression;
];

GetComponents[InputExpr_]:=Module[{Expr=InputExpr,ListExpr},
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=ToValues;
	ListExpr=Expr;
	Expr//=MatrixForm;	
ListExpr];

ObtainComponents[InputExpr_]:=Module[{Expr=InputExpr},
	Expr//=SeparateMetric[G];
	Expr//=ToBasis[$Chart];
	Expr//=TraceBasisDummy;
	Expr//=ToValues;
	Expr//=ComponentArray;
	Expr//=ToValues;
Expr];

CreateVector[LHS_,RHS_]:=Module[{TensorHead=Head@LHS,
		Indices=First@(List@@LHS),Expr},
	Block[{Print=Null},
		ComponentValue[ComponentArray@Evaluate@LHS,Evaluate@RHS];
	];
	Expr=TensorHead[-Indices];
	Block[{Print=Null},
		ComponentValue[ComponentArray@Expr,ObtainComponents@Expr];
	];
	DisplayExpression@(Grid[{{LHS->MatrixForm@GetComponents@LHS,
		Expr->MatrixForm@GetComponents@Expr}},Spacings->5]);
];

Splinter[InputExpr_]:=Module[{
		Expr=InputExpr},

	Expr//=SeparateMetric[G];
	Expr//=ToBasis[$Chart];
	Expr//=ComponentArray;
	Expr//=ToValues;
	Expr//=TraceBasisDummy;
	Expr//=ToValues;
	Expr=Map[Simplify,Expr,{2}];
	(*Expr=Map[FullSimplify,Expr,{2}];*)
Expr];

GeneralCreateComponents[InputTensor_,FullyExpandedExpression_]:=Module[{
		LeftExpr=InputTensor,
		RightExpr=FullyExpandedExpression},
	LeftExpr//=ToBasis[$Chart];
	LeftExpr//=ComponentArray;
	Block[{Print=Null},
		ComponentValue[Evaluate@LeftExpr,Evaluate@RightExpr];
	];
];

GeneralDenoteComponents[InputTensor_,FullyExpandedExpression_]:=Module[{
		Expr=FullyExpandedExpression},
	Expr//=MatrixForm;	
	(InputTensor->Expr)//DisplayExpression;
];
