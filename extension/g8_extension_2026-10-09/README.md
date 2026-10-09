# g8_extension_2026-10-09 — G8 on the v8 working documents, and an audit of the numbers added since the freeze

`10_new_system/scripts/g8_number_traceability.py` was run **unmodified** on the accepted-changes text of
`KDSS_A2AR_AIIM_v8_WORKING.docx` and `Supplementary_Information_v8_WORKING.docx` (`inputs/` in the development tree, which holds the manuscript text and is **not deposited**; built by
`analysis/build_g8_inputs.py`: the manuscript as markdown, the supplement as a docx with tracked deletions removed,
because G8's docx reader keeps deleted text). Its output is
`verification_summary_2026-08-21/g8_number_traceability_v8_accepted_2026-10-09.json`; the original
`g8_number_traceability.json` (the v6-markdown run of 2026-10-05) was overwritten by the run and restored
byte-identical (`661e3bc6…`).

`analysis/new_number_audit.py` is a second, weaker audit written because G8 anchors a fixed list of v6-era numbers to
literal strings and cannot see a number added since. It checks every numeric token inside a tracked insertion against
the numbers in the post-freeze run artifacts (tier A), then the decision-log entries from ADR-AIIM-311 (tier B); what
matches neither is listed for a person (tier C). It detects typos and drift, not mismatched quantities, and it is not G8.
Result and caveats: ADR-AIIM-338.
