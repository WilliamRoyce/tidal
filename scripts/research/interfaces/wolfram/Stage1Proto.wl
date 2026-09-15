(* ::Package:: *)
(* cspell:words WAVEOP *)

(* ==============================================================================
   MODULE:   R1Proto` — research prototype of the committed Stage-1 package (R-1, #566)

   PURPOSE:  Take ONE theory as data (an Association imported from WXF) and derive its
             PSALTer spectrum, with no generated code anywhere:
               validate names -> parse each operator string HELD -> whitelist its symbols
               -> declare couplings and fields in Global` -> build L = Sum c_i O_i
               -> ParticleSpectrum -> export spectrum.wxf.

   INPUT:    <|"fingerprint" -> hex, "fields" -> <|name -> <|"indices" -> {"-a",...},
             "print_as" -> s, "print_source_as" -> s|>|>, "lagrangian" -> <|coupling ->
             operator string|>, "derivation" -> <|"max_laurent_depth" -> n|>|>

   WHY THESE RULES (evidence: docs/cosmology/r1_planning_record.md §2.3, §2.5, App. B2):
     * User text is parsed with ToExpression[..., HoldComplete] and never evaluated before
       every symbol is on a whitelist, so `m_phi` (Pattern), `Print[..]`, `Quit[]` or a
       coupling inside an operator are refused without running.
     * Field and coupling symbols are created in Global`. PSALTer builds a context from
       ToString[field] (DefField.m:37) and looks it up (ValidateLagrangian.m:37); a symbol
       born in a private context fails that lookup and mis-creates its source tensor
       (DefField.m:52).
     * ParticleSpectrum throws uncaught (ParticleSpectrum.m:23); it is wrapped in Catch so
       a PSALTer refusal becomes a Failure, not an exit-0 abort (#561).
     * The theory name is derived from the fingerprint, so it cannot collide with a symbol
       PSALTer exports (UpdateTheoryAssociation.m:8).

   CONVENTIONS: PSALTer's — signature (+,-,-,-), epsilon_0123 = +1 (PSALTer.m:101).
   Not production code: the production package is I-S1B's.
   ============================================================================== *)

BeginPackage["R1Proto`", {"xAct`xCore`", "xAct`xPerm`", "xAct`xTensor`", "xAct`xCoba`",
  "xAct`xTras`", "xAct`SymManipulator`", "xAct`PSALTer`"}];

RunStage1::usage = "RunStage1[theory, outDir] validates a theory Association, runs ParticleSpectrum and writes spectrum.wxf into outDir. RunStage1[theory, outDir, \"ValidateOnly\" -> True] stops after validation. Returns an Association on success, a Failure otherwise.";
ValidateTheory::usage = "ValidateTheory[theory] parses and whitelists a theory Association without evaluating user text. Returns an Association of held parts or a Failure.";
$R1ProtoSchema::usage = "$R1ProtoSchema is the fingerprint/manifest schema this prototype writes.";

Begin["`Private`"];

$R1ProtoSchema = "r1-proto-1";
$conventions = <|"signature" -> {1, -1, -1, -1}, "epsilon0123" -> 1|>;
$parseContextPath = {"xAct`PSALTer`", "xAct`xTensor`", "xAct`xCore`", "System`"};

fail[tag_String, msg_String] := Failure[tag, <|"MessageTemplate" -> msg|>];

(* ---- names ---------------------------------------------------------------- *)

(* Parse a bare name held, in Global`, with only PSALTer/xAct/System visible. *)
heldName[name_String] := Block[{$Context = "Global`", $ContextPath = $parseContextPath},
  ToExpression[name, InputForm, HoldComplete]];

contextOf[HoldComplete[s_Symbol]] := Context[Unevaluated[s]];
contextOf[_] := None;

(* A coupling or field name: letters and digits, starting with a letter (no underscore:
   `m_phi` would parse as Pattern[m, Blank[phi]]), and it must be a fresh Global` symbol,
   not one PSALTer/xAct already owns (e.g. G, CD, V). *)
checkName[kind_String, name_String] := Which[
  !StringMatchQ[name, RegularExpression["[A-Za-z][A-Za-z0-9]*"]],
    fail["BadName", kind <> " name " <> ToString[name, InputForm] <>
      " must be letters and digits starting with a letter (no underscores: m_phi parses as a pattern)."],
  contextOf[heldName[name]] =!= "Global`",
    fail["NameCollision", kind <> " name " <> name <> " is already a symbol in " <>
      ToString[contextOf[heldName[name]]] <> "; choose another name."],
  True, heldName[name]];
checkName[kind_String, name_] :=
  fail["BadName", kind <> " name must be a string, got " <> ToString[name, InputForm]];

(* ---- operators ------------------------------------------------------------ *)

indexSymbols := indexSymbols = Symbol["xAct`PSALTer`" <> #] & /@ CharacterRange["a", "z"];

(* Every symbol in a held expression, heads included, as HoldComplete[sym] (unevaluated). *)
heldSymbols[h_HoldComplete] := Cases[h, s_Symbol :> HoldComplete[s], {1, Infinity}, Heads -> True];

parseOperator[coupling_String, str_String, fieldHeld_List, couplingHeld_List] := Module[
  {held, syms, allowed, bad, atoms},
  If[!SyntaxQ[str],
    Return[fail["SyntaxError", "operator for " <> coupling <> " is not one syntactically valid expression (error near character " <>
      ToString[SyntaxLength[str]] <> ")."]]];
  held = Block[{$Context = "Global`", $ContextPath = $parseContextPath},
    ToExpression[str, InputForm, HoldComplete]];
  If[!MatchQ[held, HoldComplete[_]],
    Return[fail["NotOneExpression", "operator for " <> coupling <> " must be exactly one expression."]]];
  syms = heldSymbols[held];
  (* the outer HoldComplete head is itself collected once; any further one came from the user *)
  If[Count[syms, HoldComplete[HoldComplete]] > 1,
    Return[fail["ForbiddenSymbol", "operator for " <> coupling <> " contains HoldComplete."]]];
  syms = DeleteDuplicates@DeleteCases[syms, HoldComplete[HoldComplete]];
  allowed = Join[HoldComplete /@ {Plus, Times, Power, xAct`PSALTer`CD}, HoldComplete /@ indexSymbols, fieldHeld];
  bad = Complement[syms, allowed];
  If[IntersectingQ[bad, couplingHeld],
    Return[fail["CouplingInOperator", "operator for " <> coupling <> " contains the coupling(s) " <>
      StringRiffle[ToString /@ (Intersection[bad, couplingHeld] /. HoldComplete[s_] :> SymbolName[Unevaluated[s]]), ", "] <>
      "; each coupling must multiply its operator from outside (S = Sum c_i O_i), so the action stays linear in the couplings."]]];
  If[bad =!= {},
    Return[fail["ForbiddenSymbol", "operator for " <> coupling <> " uses symbol(s) that are not declared fields, CD, indices a-z or + * ^: " <>
      StringRiffle[(bad /. HoldComplete[s_] :> Context[Unevaluated[s]] <> SymbolName[Unevaluated[s]]), ", "] <> "."]]];
  atoms = Cases[held, x : (_Real | _String | _Complex) :> HoldComplete[x], {1, Infinity}, Heads -> True];
  If[atoms =!= {},
    Return[fail["InexactOrStringAtom", "operator for " <> coupling <> " contains non-integer atoms (use exact rationals like 1/2): " <>
      ToString[atoms, InputForm] <> "."]]];
  held];

(* ---- validation ------------------------------------------------------------ *)

ValidateTheory[theory_Association] := Catch[Module[
  {fields, lag, deriv, fieldNames, couplingNames, fieldHeld, couplingHeld, indices, ops, depth, fp},
  If[Complement[Keys[theory], {"fingerprint", "fields", "lagrangian", "derivation"}] =!= {},
    Throw[fail["UnknownKeys", "unknown theory keys: " <> ToString[Complement[Keys[theory], {"fingerprint", "fields", "lagrangian", "derivation"}]]], $tag]];
  {fields, lag, deriv, fp} = Lookup[theory, {"fields", "lagrangian", "derivation", "fingerprint"}, Missing[]];
  If[!AssociationQ[fields] || Length[fields] == 0, Throw[fail["NoFields", "theory needs a non-empty fields mapping."], $tag]];
  If[!AssociationQ[lag] || Length[lag] == 0, Throw[fail["NoLagrangian", "theory needs a non-empty lagrangian mapping."], $tag]];
  If[!StringQ[fp] || !StringMatchQ[fp, RegularExpression["[0-9a-f]{64}"]], Throw[fail["NoFingerprint", "theory is missing its fingerprint."], $tag]];
  If[Names["xAct`PSALTer`a"] === {}, Throw[fail["NoIndices", "PSALTer's abstract indices are not loaded (Needs[\"xAct`PSALTer`\"] first)."], $tag]];
  fieldNames = Keys[fields]; couplingNames = Keys[lag];
  fieldHeld = checkName["field", #] & /@ fieldNames;
  couplingHeld = checkName["coupling", #] & /@ couplingNames;
  Scan[If[FailureQ[#], Throw[#, $tag]] &, Join[fieldHeld, couplingHeld]];
  If[IntersectingQ[fieldNames, couplingNames], Throw[fail["NameCollision", "a name is used as both field and coupling."], $tag]];
  indices = Map[
    Function[f, With[{spec = fields[f]},
      If[!AssociationQ[spec] || !ListQ[spec["indices"]] ||
          !AllTrue[spec["indices"], StringQ[#] && StringMatchQ[#, RegularExpression["-?[a-z]"]] &],
        Throw[fail["BadIndices", "field " <> f <> " needs indices like [-a] or [-a, -b]."], $tag]];
      Map[If[StringStartsQ[#, "-"], -Symbol["xAct`PSALTer`" <> StringDrop[#, 1]], Symbol["xAct`PSALTer`" <> #]] &, spec["indices"]]]],
    fieldNames];
  ops = KeyValueMap[parseOperator[#1, #2, fieldHeld, couplingHeld] &, lag];
  Scan[If[FailureQ[#], Throw[#, $tag]] &, ops];
  ops = MapThread[If[MatchQ[#1, HoldComplete[0]], Throw[fail["ZeroOperator", "operator for " <> #2 <> " is zero."], $tag], #1] &, {ops, couplingNames}];
  depth = Lookup[Replace[deriv, Except[_Association] -> <||>], "max_laurent_depth", 1];
  If[!MemberQ[{1, 2, 3}, depth], Throw[fail["BadDerivation", "max_laurent_depth must be 1, 2 or 3 (PSALTer.m:82)."], $tag]];
  <|"fingerprint" -> fp, "theoryName" -> "Theory" <> StringTake[fp, 12],
    "fieldNames" -> fieldNames, "fieldHeld" -> fieldHeld, "indices" -> indices,
    "fieldSpecs" -> fields, "couplingNames" -> couplingNames, "couplingHeld" -> couplingHeld,
    "operators" -> ops, "maxLaurentDepth" -> depth|>
  ], $tag];
ValidateTheory[other_] := fail["NotAnAssociation", "theory must be an Association, got " <> ToString[Head[other]]];

(* ---- derivation ------------------------------------------------------------ *)

Options[RunStage1] = {"ValidateOnly" -> False};

RunStage1[theory_, outDir_String, OptionsPattern[]] := Module[
  {v, lagrangian, timing, status, assoc, waveOp, name, out},
  v = ValidateTheory[theory];
  If[FailureQ[v], Return[v]];
  Print["STAGE1_VALIDATED couplings=", v["couplingNames"], " fields=", v["fieldNames"]];
  If[TrueQ[OptionValue["ValidateOnly"]], Return[<|"validated" -> True, "couplings" -> v["couplingNames"]|>]];
  name = v["theoryName"];
  (* declare in Global` — see header: PSALTer derives contexts from ToString[symbol] *)
  Scan[Function[h, ReleaseHold[h /. HoldComplete[s_] :> Hold[xAct`xTensor`DefConstantSymbol[s]]]], v["couplingHeld"]];
  MapThread[
    Function[{h, idx, f}, With[{spec = v["fieldSpecs"][f]},
      ReleaseHold[h /. HoldComplete[s_] :> Hold[xAct`PSALTer`DefField[s @@ idx,
        PrintAs -> Lookup[spec, "print_as", f], PrintSourceAs -> Lookup[spec, "print_source_as", "j"]]]]]],
    {v["fieldHeld"], v["indices"], v["fieldNames"]}];
  lagrangian = Total@MapThread[ReleaseHold[#1] * ReleaseHold[#2] &, {v["couplingHeld"], v["operators"]}];
  If[PossibleZeroQ[lagrangian], Return[fail["ZeroLagrangian", "the Lagrangian evaluates to zero."]]];
  {timing, status} = AbsoluteTiming@Catch[
    xAct`PSALTer`ParticleSpectrum[lagrangian, TheoryName -> name,
      MaxLaurentDepth -> v["maxLaurentDepth"], ShowPropagator -> False];
    "completed"];
  If[status =!= "completed", Return[fail["PSALTerRefused", "ParticleSpectrum threw (see its message above): " <> ToString[status, InputForm]]]];
  assoc = Symbol["Global`" <> name];
  If[!AssociationQ[assoc] || !KeyExistsQ[assoc, xAct`PSALTer`WaveOperator],
    Return[fail["NoHarvest", "ParticleSpectrum completed but Global`" <> name <> " holds no WaveOperator."]]];
  waveOp = assoc[xAct`PSALTer`WaveOperator];
  out = <|"schema" -> $R1ProtoSchema, "fingerprint" -> v["fingerprint"], "theory_name" -> name,
    "conventions" -> $conventions, "couplings" -> Sort[v["couplingNames"]],
    "psalter_version" -> xAct`PSALTer`Private`$Version[[1]], "wolfram_version" -> $VersionNumber,
    "particle_spectrum_wall_s" -> Round[timing, 0.01],
    "wave_operator" -> waveOp, "pseudo_determinant" -> Lookup[assoc, xAct`PSALTer`PseudoDeterminant, Missing[]]|>;
  Export[FileNameJoin[{outDir, "spectrum.wxf"}], out, "WXF"];
  Print["STAGE1_THEORY_NAME=", name];
  Print["STAGE1_WALL_S=", Round[timing, 0.01]];
  Print["STAGE1_SECTORS=", Length[waveOp]];
  Print["STAGE1_WAVEOP_SHA256=", Hash[BinarySerialize[waveOp], "SHA256", "HexString"]];
  Print["STAGE1_PSALTER_VERSION=", xAct`PSALTer`Private`$Version[[1]]];
  <|"theoryName" -> name, "wallSeconds" -> timing|>];

End[];
EndPackage[];
