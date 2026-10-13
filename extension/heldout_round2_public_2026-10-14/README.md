# Second blind round of fault cases (N01 to N07), 2026-10-13 and 2026-10-14

Seven fault cases were written in a new session of a different model, without access to the project, and evaluated once against the frozen evaluators of `pipeline/` under rules registered before the case file was opened (Section 3.1 of the article; Supplementary Section 2).

| file | content |
|---|---|
| `PUBLIC_PREREG_ADR452.md` | the registration and addenda 2 to 4 |
| `evaluator_fingerprint.txt` | SHA-256 of the evaluator, contract and engine files of the primary arm and of the added-check files of the secondary arm, computed before the case file was opened; the primary-arm files equal those under `pipeline/` |
| `results_public.json` | gate, check, declared expectation, states and classification of each case |

**Withheld**, for the reason given for the first round (the cases and the workspaces would allow the faults to be reconstructed): the case file, the faulted workspaces, the mapping, the two scripts and the raw results with their evidence. They are available from the corresponding author on request. SHA-256 of what is withheld, so that anything later supplied can be checked against it:

- case file as registered (final version): `e2289991eef00bf8955055838a2151eb3b5fedc1f930f46774bc13ac522cfdc0`
- mapping of the cases to the inputs of the evaluators (MAPPING_ADR452.md): `9da7253d98a558ff5b06e8f078f40d09c23f89c60602e50f0b53ee06d7071b36`
- primary-arm script (analysis/run_round2.py): `31496bc97d99ec9a046dec32c591913b5d4221fb485e7d138e3a60d9b009477e`
- secondary-arm script (analysis/run_round2_secondary.py): `b4003014c8268153ef915c0fa1c9d9d992525a0e256966cc5411d3c30085cad6`
- raw primary result with evidence (results_primary.json): `57c546fd4e7022ed7214efdcdb7dc2044d049eaa811c844534598f4c56bf7ce8`
- raw secondary result (results_secondary.json): `3b1bf02e810251fd0fb0b60dc9873e152932e37f799148a8e62162917eed1424`
- addendum 1 (ADDENDUM_1_ADR452.md): `f2b493819798732cee475b5f40212fbc7812432939dd71e23d47e12d1b5a5421`
- first version of the case file (replaced in two cases before the run): `9498e333f55a870c14b41eb320d7274de9a0a88fce2c5d33bc28c8d56fe850cb`

The run was made once by the Claude agent session that wrote the evaluators, so it is not independent; see the article.
