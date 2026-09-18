<!-- To file at https://github.com/THelpin/xBrauer_Bundle/issues/new . The first line after
this comment is the title; everything after it is the body. Every snippet was run as written by
check_snippets.py; the claims table is claims.md. -->
**Title:** After loading xBrauer, `ContractMetric` stops contracting an induced metric into a tensor without `Master`, and `SeparateMetric` stops separating a derivative index

**Environment:** Wolfram Engine 14.3.0; xAct 1.3.0 (xTensor 1.3.0); xBrauer 1.1.0 (2023-11-29), `xBrauer_Bundle` at `48be67e1`, unmodified.

We met this through xPand, whose perturbation fields are tensors on a slice with no `Master`: once xBrauer is loaded (for instance through xMAG), xPand's splits keep terms that should vanish by transversality, with no message. The reproduction needs only xTensor:

```wolfram
Needs["xAct`xTensor`"]
DefManifold[M, 4, {a, b, c, d}];
DefMetric[-1, g[-a, -b], CD, {";", "\[Del]"}];
DefTensor[v[a], M];
AutomaticRules[v, MakeRule[{v[-a] v[a], -1}], Verbose -> False];
DefMetric[1, q[-a, -b], cdq, {"|", "D"}, InducedFrom -> {g, v}];
DefTensor[T[-a], M, OrthogonalTo -> {v[a]}];
DefTensor[S[-a], M, OrthogonalTo -> {v[a]}, Master -> q];
DefTensor[f[], M];
ContractMetric[q[a, b] T[-a]]    (* T[b] *)
SeparateMetric[g][cdq[a]@f[], a]    (* separated with q, a dummy index appears *)
Needs["xAct`xBrauer`"]
ContractMetric[q[a, b] T[-a]]    (* q[a, b]*T[-a] *)
ContractMetric[q[a, b] T[-a], q]    (* q[a, b]*T[-a] *)
ContractMetric[q[a, b] S[-a]]    (* S[b] *)
ContractMetric[g[a, b] T[-a]]    (* T[b] *)
SeparateMetric[g][cdq[a]@f[], a]    (* cdq[a][f[]] *)
```

### 1. The contraction

`xBrauer.m:1938-1940` issue three `Unset`s on xTensor's `ContractMetric1` rules, and `:1947-1974` install twelve replacements with an extra condition — for a tensor slot

```wolfram
MetricOfTensor[tensor] === metric || (MetricOfTensor[tensor] === Null && FirstMetricQ[metric])
```

and the same with `MasterOf[covd]` for a derivative slot. For a tensor without `Master`, `MetricOfTensor[tensor]` is `Null` (`:315`), so the condition reduces to `FirstMetricQ[q]`, which is `False` for an induced metric: it is never the bundle's first. The contraction is refused even when the metric is named. A tensor with `Master -> q` passes through `MetricOfTensor === q`, and contraction with the first metric `g` is unaffected.

One thing that can hide this: a tensor declared with `ProjectedWith -> {q[a, -b]}` gets a projector rule of its own that contracts `q[a, b] T[-a]` as soon as the product is built, before `ContractMetric` runs, so it looks unaffected.

Suggested (not tested): let the fallback accept an induced metric as well, `(MetricOfTensor[tensor] === Null && (FirstMetricQ[metric] || InducedMetricQ[metric]))`, with the same change for `MasterOf[covd]`.

### 2. The separation

`xBrauer.m:1866` removes xTensor's catch-all `SeparateMetric1` definition (`xTensor.m:8682`) and dispatches induced and frozen metrics to `SeparateMetric2`, which `:1874` replaces; its derivative branch separates only when `MasterOf[covd] === metric`. So `SeparateMetric[g]` on an index of the induced derivative `cdq` (whose `MasterOf` is `q`) now returns the expression unchanged, where before it separated that index with `q`. xPand's `IndicesDown` is built on `SeparateMetric` and changes with it.

**Not related:** #2 concerns the indexed form `MetricOfTensor[tensor[inds]]`; the condition above uses the head form.
