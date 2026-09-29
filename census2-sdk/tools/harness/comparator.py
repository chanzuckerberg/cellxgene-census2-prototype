"""The comparator (task 2.3): decide whether two canonical results are the same answer.

Both sides are in the canonical form of `canonical.py`: rows ordered by identity,
categoricals decoded, X in the sorted coordinate space. The comparator normalizes
nothing, so every difference it reports is a real difference.

The rules:

1. Row identity comes first. `obs` rows are keyed by `(dataset_id, observation_joinid)`
   and `var` rows by `feature_id`. A missing, extra, or duplicated key is a row-set
   finding that names the first affected identity. A right side whose keys are a strict
   subset of the left is its own kind, `strict-subset`, because that is the shape an
   exclusion bug produces.
2. Values are compared only where identities align. A mismatch names the first cell
   that differs: its identity, its column, and both values.
3. Floats compare exactly, including the raw counts in X. The one exception is the
   layers in `FLOAT_TOLERANT` (`normalized`), which compare within `rtol=1e-6`,
   because the operation order of derived floating-point arithmetic can differ
   between stacks.
4. A column present on one side only is a finding.
5. Two integer columns compare exactly. A float column, or an int/float pair,
   compares NaN-aware (NaN equals NaN) in float64, and the dtype disagreement is
   tolerated when the values compare equal. The Census 1 SDK's dtype for a nullable
   integer column (the spatial `array_row`, `array_col`, `in_tissue`) depends on
   whether the selected rows include a null (evidence 18.1). A kind disagreement
   (string vs numeric) is a finding.
6. An axis that neither side holds, such as `obs` in a presence result, is not
   compared. An axis that either side holds must carry its identity columns on both
   sides, or the axis is a finding.

`compare` returns a list of findings. An empty list means the same answer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from harness.canonical import Arrays, decode_strings

FLOAT_TOLERANT = {"normalized": 1e-6}  # layer name -> rtol; everything else is exact


@dataclass(frozen=True)
class Finding:
    kind: str  # row-set | strict-subset | column-set | dtype | shape | value
    where: str  # the array or axis the finding is about
    detail: str  # names the first affected identity or column and both values

    def __str__(self) -> str:
        return f"[{self.kind}] {self.where}: {self.detail}"


def _column(arrays: Arrays, name: str) -> np.ndarray | None:
    """A column as comparable values: decoded strings, or the numeric array itself."""
    if f"{name}.offsets" in arrays:
        strings = decode_strings(arrays[f"{name}.offsets"], arrays[f"{name}.bytes"])
        return np.asarray(strings, dtype=object)
    return arrays.get(name)


def _identity(arrays: Arrays, axis: str) -> np.ndarray | None:
    if axis == "obs":
        ds, oj = (
            _column(arrays, "obs.dataset_id"),
            _column(arrays, "obs.observation_joinid"),
        )
        if ds is None or oj is None:
            return None
        return np.asarray([f"{a}|{b}" for a, b in zip(ds, oj, strict=True)])
    fid = _column(arrays, "var.feature_id")
    return None if fid is None else np.asarray(fid, dtype=str)


def _column_names(arrays: Arrays, axis: str) -> set[str]:
    out = set()
    for key in arrays:
        if key.startswith(f"{axis}."):
            out.add(key.removesuffix(".offsets").removesuffix(".bytes"))
    return out


def _compare_axis(axis: str, left: Arrays, right: Arrays) -> list[Finding]:
    findings: list[Finding] = []

    lid, rid = _identity(left, axis), _identity(right, axis)
    if lid is None or rid is None:
        if (lid is None) != (rid is None):
            findings.append(
                Finding("column-set", axis, "identity columns present on one side only")
            )
        elif _column_names(left, axis) or _column_names(right, axis):
            findings.append(
                Finding("column-set", axis, "identity columns missing on both sides")
            )
        return findings

    # 1. row sets, by identity multiset
    if not np.array_equal(lid, rid):
        lset, rset = set(lid), set(rid)
        missing, extra = sorted(lset - rset), sorted(rset - lset)
        if missing and not extra and len(rid) < len(lid):
            findings.append(
                Finding(
                    "strict-subset",
                    axis,
                    f"right is a strict subset: {len(lid) - len(rid)} rows missing, "
                    f"first missing identity {missing[0]!r}",
                )
            )
        elif missing or extra:
            what = (
                f"first missing {missing[0]!r}"
                if missing
                else f"first extra {extra[0]!r}"
            )
            findings.append(
                Finding(
                    "row-set",
                    axis,
                    f"{len(missing)} missing / {len(extra)} extra rows; {what}",
                )
            )
        else:  # same key set, different multiplicity; canonical order makes this duplication
            keys_l, counts_l = np.unique(lid, return_counts=True)
            _, counts_r = np.unique(rid, return_counts=True)
            diff = np.nonzero(counts_l != counts_r)[0]
            key = keys_l[diff[0]] if len(diff) else lid[0]
            findings.append(Finding("row-set", axis, f"duplicated identity {key!r}"))
        return findings  # values are meaningless while rows differ

    # 2. column sets
    lcols, rcols = _column_names(left, axis), _column_names(right, axis)
    for col in sorted(lcols - rcols):
        findings.append(Finding("column-set", col, "present on left only"))
    for col in sorted(rcols - lcols):
        findings.append(Finding("column-set", col, "present on right only"))

    # 3. values, the first differing cell named by identity
    for col in sorted(lcols & rcols):
        lv, rv = _column(left, col), _column(right, col)
        assert lv is not None and rv is not None
        lk, rk = lv.dtype.kind, rv.dtype.kind
        if {lk, rk} <= set("iu"):  # rule 5: integers compare exactly
            neq = lv != rv
        elif {lk, rk} <= set("iuf"):  # rule 5: floats, NaN-aware, dtype-tolerant
            l64, r64 = lv.astype(np.float64), rv.astype(np.float64)
            neq = ~((l64 == r64) | (np.isnan(l64) & np.isnan(r64)))
        elif lv.dtype != rv.dtype and (lk not in "OU" or rk not in "OU"):
            findings.append(Finding("dtype", col, f"{lv.dtype} vs {rv.dtype}"))
            continue
        else:
            neq = lv != rv
        if neq.any():
            i = int(np.argmax(neq))
            findings.append(
                Finding(
                    "value",
                    col,
                    f"identity {lid[i]!r}: left {lv[i]!r} vs right {rv[i]!r}",
                )
            )
    return findings


def _compare_layer(
    name: str,
    left: Arrays,
    right: Arrays,
    obs_id: np.ndarray | None,
    var_id: np.ndarray | None,
) -> list[Finding]:
    prefix = f"X.{name}."
    lshape, rshape = left[prefix + "shape"], right[prefix + "shape"]
    if not np.array_equal(lshape, rshape):
        return [Finding("shape", f"X.{name}", f"{tuple(lshape)} vs {tuple(rshape)}")]
    for part in ("indptr", "indices"):
        if not np.array_equal(left[prefix + part], right[prefix + part]):
            diff = int(np.argmax(left[prefix + part] != right[prefix + part]))
            return [
                Finding(
                    "value",
                    f"X.{name}.{part}",
                    f"sparsity structure differs at position {diff}",
                )
            ]
    lv, rv = left[prefix + "data"], right[prefix + "data"]
    rtol = FLOAT_TOLERANT.get(name)
    equal = np.allclose(lv, rv, rtol=rtol, atol=0) if rtol else np.array_equal(lv, rv)
    if equal:
        return []
    neq = ~np.isclose(lv, rv, rtol=rtol, atol=0) if rtol else lv != rv
    k = int(np.argmax(neq))
    row = int(np.searchsorted(left[prefix + "indptr"], k, side="right") - 1)
    col = int(left[prefix + "indices"][k])
    cell = obs_id[row] if obs_id is not None and row < len(obs_id) else f"row {row}"
    gene = var_id[col] if var_id is not None and col < len(var_id) else f"col {col}"
    return [
        Finding(
            "value",
            f"X.{name}",
            f"cell {cell!r} x gene {gene!r}: left {lv[k]!r} vs right {rv[k]!r}",
        )
    ]


def compare(left: Arrays, right: Arrays) -> list[Finding]:
    findings: list[Finding] = []
    findings += _compare_axis("obs", left, right)
    findings += _compare_axis("var", left, right)

    row_findings = {f.kind for f in findings} & {"row-set", "strict-subset"}
    lnames = {k.split(".")[1] for k in left if k.startswith("X.")}
    rnames = {k.split(".")[1] for k in right if k.startswith("X.")}
    for name in sorted(lnames ^ rnames):
        side = "present on left only" if name in lnames else "present on right only"
        findings.append(Finding("column-set", f"X.{name}", side))
    if not row_findings:  # X positions are meaningless while row sets differ
        obs_id, var_id = _identity(left, "obs"), _identity(left, "var")
        for name in sorted(lnames & rnames):
            findings += _compare_layer(name, left, right, obs_id, var_id)
    return findings
