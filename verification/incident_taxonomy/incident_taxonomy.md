# Incident taxonomy

Referenced from Section 3.6 of `manuscript_AIIM_v6.md`. Fifteen incidents
recorded in the project decision log (Supplementary S5) during construction
and execution of the pipeline, classified by the layer each occurred in, what
detected it, and how it resolved. Three are scientific, six are
infrastructure, six are agent-reasoning errors.

| Layer | Incident | Detector | Resolution |
|---|---|---|---|
| Scientific | Column-selection defect silently validating the wrong field | Contract evaluator | Corrected; no synthetic fixture had triggered it |
| Scientific | Pose-RMSD computation admitting a mispositioned pose | Contract evaluator | Corrected |
| Scientific | Net system charge +0.996 in the membrane build | Pre-run integrity check | One Cl⁻ added; a retained crystallographic Na⁺ had not been counted |
| Infrastructure | Per-compound memory growth toward exhaustion | Resident-set comparison across restarts | 30-minute forced-restart wrapper |
| Infrastructure | Host restarts terminating all docking processes | Process check plus boot-time inspection | Idempotent resume; zero rows lost |
| Infrastructure | CPU oversubscription — three chunks each claiming all cores | Load average 38–52 on 20 cores | Explicit thread cap; load fell to 21.6 |
| Infrastructure | An edited shell script not reloaded by an already-running loop | Inspection of the running process command line | Wrapper processes restarted explicitly |
| Infrastructure | Three glycosides never converging | Identical progress counts across restart cycles | Documented exclusion with recorded reason |
| Infrastructure | Stray directory created by a host path used inside the Linux guest | Directory listing | Removed; path resolution now probes for a known file |
| Agent reasoning | Literature precedent asserted without a search, resting on non-human-receptor data | Cross-check against the project's own frozen labels | Claim withdrawn |
| Agent reasoning | Process loss attributed to system sleep | Operating-system event log inspected | Corrected to user-initiated restarts; the proposed remediation was withdrawn as ineffective |
| Agent reasoning | Journal scope judged from the scope statement alone | The journal's published record searched | Reversed |
| Agent reasoning | "This force field reports no parameter-penalty score" | Tool output read directly | Corrected; one ligand torsion scores 87.0 |
| Agent reasoning | "The simulation package is not installed" | Only one of four guest environments had been checked | Corrected |
| Agent reasoning | A membrane-build failure attributed to the packing stage on the strength of a clash count measured on the *post*-parameterisation file, not the packing output | The two files measured separately | Attribution corrected — the packing output held 28 sub-0.7 Å contacts, the parameterised file 119; most were introduced when hydrogens and virtual sites were built onto it |

## Two patterns

**Every agent-reasoning error was caught by checking a primary source the
agent had not consulted** — an event log, a tool's own output, a journal's
back catalogue, the project's frozen labels — and none was caught by the
agent re-reading its own reasoning.

The last two scientific rows were caught by checks that exist only because
the workflow requires artifacts to be verified before use; both would
otherwise have silently degraded a downstream result rather than failing
visibly.

These are reported because a system that reports only its successes cannot
substantiate a claim about error detection.
