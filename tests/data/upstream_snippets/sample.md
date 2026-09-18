# A report with two snippets

<!-- setup: DefManifold[M, 4, {a, b, c, d}]; -->
```wolfram
Needs["xAct`xTensor`"]; $RiemannSign    (* 1 *)
q[a, b] T[-a]                          (* q[a, b]*T[-a] *)
```

Some text.

```wolfram
ToCanonical[x]      (* -K[-b, a, -a], unreduced *)
ToCanonical[y]      (* still unreduced *)
f[]
```
