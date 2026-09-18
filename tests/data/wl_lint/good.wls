f[x_, y_: 1] := (x
   + y);
h[e_, cands_List: {1, -1}] := e;
f[1];
f[1, 2];
h[3];
Needs["xAct`xTensor`"];
Check[Needs["xAct`xPand`"], failed = True];
Print[System`$Version]
Print[xAct`xPand`$Version]
