"""Compare the digests that `digest.py` wrote for three configurations.

Usage: python compare_digests.py MAIN.json PER_DATASET.json CONSOLIDATED.json [OUT.json]

A case is the same when every digest but `stores_opened` agrees. For a case that
differs, the output names the parts that differ.
"""

import json
import pathlib
import sys

names = ("main", "per_dataset", "consolidated")
runs = {
    k: json.loads(pathlib.Path(p).read_text())
    for k, p in zip(names, sys.argv[1:4], strict=True)
}
out = {}
for case in runs["main"]:
    vals = {k: runs[k].get(case) for k in names}
    strip = {
        k: None if v is None else {x: y for x, y in v.items() if x != "stores_opened"}
        for k, v in vals.items()
    }
    same = strip["main"] == strip["per_dataset"] == strip["consolidated"]
    diff = []
    if not same and all(strip.values()):
        for part in ("n_obs", "X", "normalized"):
            if len({json.dumps(s[part]) for s in strip.values()}) > 1:
                diff.append(part)
        for part in ("obs", "var"):
            cols = set().union(*(s[part] for s in strip.values()))
            for col in sorted(cols):
                if len({s[part].get(col) for s in strip.values()}) > 1:
                    diff.append(f"{part}.{col}")
    out[case] = {
        "same": same,
        "n_obs": None if vals["main"] is None else vals["main"]["n_obs"],
        "stores_opened": {k: (v or {}).get("stores_opened") for k, v in vals.items()},
        "differs": diff,
    }
    print(case, same, out[case]["n_obs"], out[case]["stores_opened"], diff)
if len(sys.argv) > 4:
    pathlib.Path(sys.argv[4]).write_text(json.dumps(out, indent=1))
