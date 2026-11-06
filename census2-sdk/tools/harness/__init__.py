"""The validation harness: the query corpus, the canonical result form,
the comparator, and the oracle runner.

Census 1 is the oracle. `oracle_runner.py` records its answer to every corpus query in
the canonical form. A stack under test answers the same queries, and `comparator.compare`
reports every difference between the two answers. `sdk_adapter.py` is the SDK's
adapter to this harness.
"""
