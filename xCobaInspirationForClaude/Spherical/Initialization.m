(*==================*)
(*  Initialization  *)
(*==================*)

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
