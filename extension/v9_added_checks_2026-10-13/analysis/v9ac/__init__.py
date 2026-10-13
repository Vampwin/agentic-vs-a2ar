"""ADR-AIIM-427/428: four added checks and the evaluator wrappers that carry them. Outside 10_new_system: the frozen evaluators are imported and
called unchanged, and the added checks are appended to their output. The wrappers are importable as `v9ac.G2_plus` etc. so that the frozen
contract engine can load them from a copied contract (`evaluator: v9ac.G2_plus`)."""
