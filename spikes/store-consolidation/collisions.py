"""Count the author cell names (the obs index) that occur in more than one store.

Usage: python collisions.py CATALOG
"""

import collections
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from anndata.io import read_elem
from census2_sdk import catalog as _catalog


def names(uri: str) -> set[str]:
    obs = _catalog.open_store(uri)["obs"]
    return {str(v) for v in np.asarray(read_elem(obs[str(obs.attrs["_index"])]))}


rows = _catalog.load(sys.argv[1])
uris = [str(u) for u in rows["store_uri"]]
with ThreadPoolExecutor(16) as pool:
    sets = list(pool.map(names, uris))
count: collections.Counter[str] = collections.Counter()
for s in sets:
    count.update(s)
shared = {n for n, c in count.items() if c > 1}
stores = sum(1 for s in sets if s & shared)
print({"names_in_more_than_one_store": len(shared), "stores_with_one": stores})
