"""Local type stub for the subset of admet-ai consumed by DrugOS.

The upstream distribution ships no type annotations; documented in
doc/06-technology-stack.md per the release quality gate (G2).
"""

from typing import Any

class ADMETModel:
    def __init__(self) -> None: ...
    def predict(self, smiles: list[str], verbose: bool = False) -> Any: ...
