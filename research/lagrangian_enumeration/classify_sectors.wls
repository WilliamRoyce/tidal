#!/usr/bin/env wolframscript
(* ::Script:: *)
(*
  Classify enumerated sectors by physical constraints.
  Reads enumeration_results.json and adds: ghost status, derivative order,
  perturbation order, topological flag, IBP partner identification.

  Runs on the RESULTS of the main enumeration — no re-enumeration needed.
*)

Print["=== Sector Classification ===\n"];

(* Read existing results *)
resultsFile = FileNameJoin[{DirectoryName[$InputFileName], "enumeration_results.json"}];
counts = Import[resultsFile, "RawJSON"];

Print["Loaded ", Length[counts], " sectors from enumeration_results.json"];

(* === Helper: count building blocks in sector name === *)
CountBlock[name_String, block_String] := Module[{clean, parts},
  clean = StringReplace[name, {"eps " -> "", "CS: " -> ""}];
  parts = StringTrim /@ StringSplit[clean, " x "];
  Count[parts, block]
];

(* === Ghost classification === *)
(* Safe building blocks: {R, T, F, DT} — all O(d^1 Gamma) or O(d^1 A) *)
(* Ghost-prone: {DR, DF} — O(d^2 Gamma) or O(d^2 A) *)
GhostStatus[name_String] := Module[{nDR, nDF},
  nDR = CountBlock[name, "DR"];
  nDF = CountBlock[name, "DF"];
  Which[
    nDR + nDF == 0, "safe",       (* Yang-Mills type *)
    nDR + nDF == 1, "borderline", (* linear in higher derivative — like R in GR *)
    True, "ghost_prone"           (* quadratic+ in higher derivatives *)
  ]
];

(* === Derivative order in fundamental fields === *)
(* In first-order (Palatini) PGT: Gamma independent of g *)
(* R ~ dGamma (1), T ~ Gamma (0), DT ~ dGamma (1), DR ~ d^2Gamma (2) *)
(* F ~ dA (1), DF ~ d^2A (2) *)
DerivOrderGamma[name_String] := Module[{nR, nDR, nDT, nT},
  nR = CountBlock[name, "R"] - CountBlock[name, "DR"];
  nDR = CountBlock[name, "DR"];
  nDT = CountBlock[name, "DT"];
  nT = CountBlock[name, "T"] - CountBlock[name, "DT"];
  (* T contributes 0 derivatives, R and DT contribute 1, DR contributes 2 *)
  nR * 1 + nDT * 1 + nDR * 2
];

DerivOrderA[name_String] := Module[{nF, nDF},
  nF = CountBlock[name, "F"] - CountBlock[name, "DF"];
  nDF = CountBlock[name, "DF"];
  nF * 1 + nDF * 2
];

(* === Perturbation order (B0!=0, Tbar=0, Rbar=0) === *)
EpsOrder[name_String] := Module[{totalR, totalT, totalF},
  totalR = CountBlock[name, "R"];  (* includes DR *)
  totalT = CountBlock[name, "T"];  (* includes DT *)
  totalF = CountBlock[name, "F"];  (* includes DF *)
  totalR + totalT + Max[0, totalF - 1]
];

(* === Topological identification === *)
IsTopological[name_String] := (name === "eps F x F");
(* F.Ftilde is the only total derivative. Torsionful GB and Pontryagin are NOT topological. *)

(* === IBP pair identification === *)
(* Sectors are IBP partners if one is obtained by moving D prefix between blocks *)
FindIBPPartner[name_String, allNames_List] := Module[
  {clean, hasEps, parts, i, j, newParts, partnerName, perms},
  hasEps = StringStartsQ[name, "eps "];
  clean = StringReplace[name, "eps " -> ""];
  parts = StringTrim /@ StringSplit[clean, " x "];

  Do[
    If[StringStartsQ[parts[[i]], "D"],
      Do[
        If[i != j && !StringStartsQ[parts[[j]], "D"],
          newParts = parts;
          newParts[[i]] = StringDrop[newParts[[i]], 1]; (* remove D *)
          newParts[[j]] = "D" <> newParts[[j]];        (* add D *)
          (* Try ALL permutations of the new blocks (naming order varies) *)
          perms = Permutations[newParts];
          Do[
            partnerName = If[hasEps, "eps ", ""] <> StringRiffle[perm, " x "];
            If[MemberQ[allNames, partnerName] && partnerName =!= name,
              Return[partnerName, Module]
            ],
            {perm, perms}
          ];
        ],
        {j, Length[parts]}
      ]
    ],
    {i, Length[parts]}
  ];
  None
];

IsIBPSafeRep[name_String, partner_] := Module[{},
  If[partner === None, Return[True]];
  (* Safe rep = fewer DF+DR blocks than partner *)
  (CountBlock[name, "DF"] + CountBlock[name, "DR"]) <=
  (CountBlock[partner, "DF"] + CountBlock[partner, "DR"])
];

(* === Mixing classification === *)
MixingType[name_String, count_Integer] := Module[
  {hasR, hasT, hasF, nF},
  hasR = CountBlock[name, "R"] > 0;
  hasT = CountBlock[name, "T"] > 0;
  hasF = CountBlock[name, "F"] > 0;
  nF = CountBlock[name, "F"];
  Which[
    count == 0, "vanishing",
    hasR && hasT && hasF, "h-t mixing via B0",
    hasR && hasF && !hasT, If[nF >= 2, "h-f (Gertsenshtein)", "h-f mixing"],
    hasT && hasF && !hasR, If[nF >= 2, "t-f mixing via B0", "t-f mixing"],
    hasR && hasT && !hasF, "h-t mixing (B0-independent)",
    hasR && !hasT && !hasF, "curvature self-interaction",
    hasT && !hasR && !hasF, "torsion self/kinetic",
    hasF && !hasR && !hasT, "EM self/kinetic",
    True, "unclassified"
  ]
];

(* === Apply all classifications === *)
allNames = Keys[counts];
fullClassified = Association[];

Print["\n", StringPadRight["Sector", 25], " | ",
  StringPadRight["N", 5], " | ",
  StringPadRight["Ghost", 12], " | ",
  StringPadRight["eps", 3], " | ",
  StringPadRight["Top", 3], " | ",
  StringPadRight["dGam", 4], " | ",
  StringPadRight["dA", 3], " | ",
  StringPadRight["IBP partner", 22], " | ",
  "Safe"];
Print[StringJoin[Table["-", 110]]];

Do[
  Block[{name = kv[[1]], count = kv[[2]],
    ghost, epsOrd, topo, ibpP, safeRep, dG, dA, mix},
    ghost = GhostStatus[name];
    epsOrd = EpsOrder[name];
    topo = IsTopological[name];
    ibpP = FindIBPPartner[name, allNames];
    safeRep = IsIBPSafeRep[name, ibpP];
    dG = DerivOrderGamma[name];
    dA = DerivOrderA[name];
    mix = MixingType[name, count];

    If[count > 0,
      Print[
        StringPadRight[name, 25], " | ",
        StringPadRight[ToString[count], 5], " | ",
        StringPadRight[ghost, 12], " | ",
        StringPadRight[ToString[epsOrd], 3], " | ",
        StringPadRight[If[topo, "Y", "-"], 3], " | ",
        StringPadRight[ToString[dG], 4], " | ",
        StringPadRight[ToString[dA], 3], " | ",
        StringPadRight[If[ibpP === None, "-", ibpP], 22], " | ",
        If[safeRep, "Y", "-"]
      ];
    ];

    fullClassified[name] = <|
      "count" -> count,
      "mixing" -> mix,
      "ghost_status" -> ghost,
      "eps_order" -> epsOrd,
      "topological" -> topo,
      "ibp_partner" -> If[ibpP === None, Null, ibpP],
      "ibp_safe_rep" -> safeRep,
      "deriv_order_Gamma" -> dG,
      "deriv_order_A" -> dA
    |>;
  ],
  {kv, Normal[counts]}
];

(* === Summary === *)
Print["\n=== PHYSICAL VIABILITY SUMMARY ==="];

safePhys = Select[fullClassified,
  #["count"] > 0 && #["ghost_status"] === "safe" &&
  !#["topological"] && #["ibp_safe_rep"] &];
borderPhys = Select[fullClassified,
  #["count"] > 0 && #["ghost_status"] === "borderline" &&
  !#["topological"] && #["ibp_safe_rep"] &];
ghostSectors = Select[fullClassified,
  #["count"] > 0 && #["ghost_status"] === "ghost_prone" &];
topoSectors = Select[fullClassified, #["topological"] &];
ibpRedundant = Select[fullClassified,
  #["count"] > 0 && !#["ibp_safe_rep"] &];

Print["  Safe (Yang-Mills type): ", Length[safePhys], " sectors, ",
  Total[#["count"] & /@ Values[safePhys]], " contractions"];
Print["  Borderline (linear in DR/DF): ", Length[borderPhys], " sectors, ",
  Total[#["count"] & /@ Values[borderPhys]], " contractions"];
Print["  Ghost-prone (quadratic DR/DF): ", Length[ghostSectors], " sectors, ",
  Total[#["count"] & /@ Values[ghostSectors]], " contractions"];
Print["  Topological (total derivative): ", Length[topoSectors], " sectors"];
Print["  IBP-redundant (partner counted): ", Length[ibpRedundant], " sectors"];

(* Core quadratic count *)
coreEven = Select[safePhys,
  #["eps_order"] <= 2 && #["deriv_order_Gamma"] <= 2 && #["deriv_order_A"] <= 1 &];
Print["\n  Core parity-even (R^2+T^2+F^2+RxF): ",
  Total[#["count"] & /@ Select[Values[coreEven], !StringStartsQ[#["mixing"], "eps"] &]]];

(* === Export === *)
outDir = DirectoryName[$InputFileName];

Export[FileNameJoin[{outDir, "enumeration_classified.json"}], fullClassified];
Export[FileNameJoin[{outDir, "enumeration_physical.json"}], safePhys];

Print["\nExported: enumeration_classified.json (full), enumeration_physical.json (safe only)"];
Print["=== DONE ==="];
