# Titans Drug Discovery — v0.6

A local research workbench for gene-signature matching. Inputs are public human
microglia differential-expression results and LINCS perturbation signatures.
LINCS cell lines are not microglia assays. Outputs are experimental hypotheses,
not medicines ranked for clinical use.

## Correctness and interpretation

- GMT parsing preserves empty descriptions and rejects duplicate terms.
- The context algorithm is **directional_set_ks/v2**: signed, unweighted set
  enrichment against disease-weight ranks, evaluated at both sides of hits and
  at tie-group boundaries. Reversal is positive; mimicry is negative. It is not
  Broad CLUE WTCS, NCS or tau. `weighted_connectivity` is a historical API name.
- Consensus ranking and its one-sided permutation test both use cosine reversal.
  Set reversal, XSum and mean rank are secondary descriptive diagnostics; binary
  drug signs are never turned into an alphabetical or query-dependent ranking.
- The default disease signature uses the Mathys mixed-model consensus (`mix`).
  `--model ind` retains the original cell-level selection for sensitivity checks.
  Neither is an independent cohort or newly computed donor-level pseudobulk model.
- Context nulls permute only represented query genes. The default coverage floor
  remains three genes; inspect `--min-hits 5` as a declared sensitivity analysis,
  not a way to select whichever FDR looks best.
- Cell-balanced scores, leave-one-cell-out diagnostics and best-context scores
  are descriptive. FDR tests only the overall mean across coverage-eligible
  contexts. Contexts are not independent donor replicates.
- `control_panel_ad.csv` replaces the misleading gold-standard designation.
  Its labels remain unadjudicated. AUROC is null; a separately labeled legacy
  diagnostic is retained. The DAM rescue panel is a different contract: proposed
  roles are unobserved until an iPSC-microglia assay receipt exists.
- Packets retain actual metrics, identifiers and mechanism flags. Missing values
  stay null. Unknown mechanisms are not negative toxicity findings. Unnamed
  perturbations stay in results but cannot become experiment packets.
- Result and packet schemas are v2. Complete rankings, run configuration, source
  hashes and runtime/code receipts accompany each output. Sources with recorded
  hashes are verified on read; sources without upstream hashes are labeled.
  Ambiguous HGNC aliases are unresolved rather than arbitrarily assigned.

## Run

The public inputs must already exist in `data/cache`. No patient data or network
service is required. Python 3.12 was verified using `requirements-verified.txt`;
the package supports Python >=3.10. The optional Vault identity overlay needs
NetworkX and Plotly; HGNC identities work independently.

```bash
cd /path/to/forge-workspace/microglia_repurposing
.venv/bin/python run_context_engine.py --all-harmonizome --n-perm 2000
.venv/bin/python run_engine.py --all-lincs-consensus --n-perm 2000
.venv/bin/python run_upgrade_review.py --n-perm 2000
.venv/bin/python run_rescue_panel.py
.venv/bin/python -m pytest tests -q -W error::RuntimeWarning
```

Both CLIs accept `--model mix|ind`, `--min-hits`, `--seed`, `--n-perm` and
`--output`. Without a full-panel flag they run the exploratory controls.

`run_upgrade_review.py` creates three context analyses (mixed, original,
stricter overlap), two consensus analyses, complete context scores, sensitivity
comparisons and a portable `data/review/index.html` page. Open that file locally.
The page supports searching, sorting, candidate details and review-note export.
Notes do not automatically update scores or validation. Export notes before
leaving the page; they are held only in the current page session.

## Performance

Normalized context matrices use a source/mapping/code-keyed, hash-verified NPZ
cache with pickle disabled. Context lookup and compound metadata are indexed
once. Cold and warm timings depend on input size and local hardware; timing
receipts are recorded by the review command. No measured or inferred L1000
continuous data is fabricated from binary inputs.

## Remaining research work

The DAM rescue protocol is frozen and unobserved. It does not yet supply
adjudicated iPSC-microglia labels, viability results or clinical validation.
SEA-AD transfer is an external-cohort evaluation, not independent replication.
Next expansion order: run the declared DAM assay with receipts; pin remaining
holdout identities; assess adjacent brain cells; then establish a separate
muscle program with strength, mobility and daily function endpoints. Tissue RNA
reversal and circulating protein changes must remain distinct evidence types.
Continuous LINCS ingestion needs real ranked profiles, assay-quality metadata
and a validated reference-method comparison.

## DAM-state rescue panel (v0.6)

`data/rescue/` is a wet-lab contract, not another ranking. The query is a frozen
Keren-Shaul core DAM vs homeostatic gene set, pinned to HGNC. The matrix is
human iPSC-microglia. Viability is a required gate. LINCS scores cannot
adjudicate these labels.

```bash
cd /path/to/forge-workspace/microglia_repurposing
.venv/bin/python run_rescue_panel.py
```

Donepezil and memantine are proposed negatives. Pexidartinib is a depletion
confounder. Vorinostat is a cytotoxicity confounder. Minocycline and
pioglitazone stay legacy-not-rescue. Sirolimus is a holdout hypothesis with a
pinned InChIKey; GW3965 and MCC950 remain identity-unresolved and cannot become
eligible controls. Vehicle, LPS/IFN-γ, Aβ oligomers and the viability readout
are assay controls, not ranking labels. Rescue AUROC stays null until viable
observations exist for both a rescued and a not-rescued ranking control.

Open **DAM rescue panel** in the dashboard to inspect the freeze. Observations
require an assay receipt. Changing the gene set or matrix requires a new
protocol hash; do not edit the frozen files to look prespecified.

## Desktop dashboard (v0.6)

Double-click **Launch Titans Drug Discovery.command** on the Desktop. The launcher
starts a loopback-only server on port 8879, verifies app/version/workspace identity,
and reuses a matching running instance. It checks a code fingerprint and can restart
an older instance only after verifying the listener's user and exact project command.
It does not kill another app occupying the port.

The dashboard reads verified analyses from `data/review/`, fetches candidates in
pages of 30 and loads evidence details on selection. Reviews persist separately in
`data/reviews.sqlite3`, with an append-only edit history. Observed results require
an outcome and an assay receipt/reference. Reviews do not modify analysis scores
or confer biological validation. Export reviews from the dashboard for backup.

Run `run_upgrade_review.py` to regenerate the canonical analyses, then click
**Refresh analyses**. The replacement is validated before it becomes active;
failed refreshes leave the previous workspace available. Analyses are archived
in SQLite on successful load. Reviews remain attached to their original analysis
hash, and the saved-review queue can open the original snapshot. Analyses that
were overwritten before v0.5 first loaded them cannot be recovered automatically.
Browser drafts are separate from saved reviews and survive refresh/reload; saving
adds the note to the SQLite review history. Browser drafts are not shared between browsers.

## External cohort evaluation

Open **Replication & benchmarks** to compare the frozen Mathys development run
with six declared SEA-AD MTG Micro-PVM states. All six published aggregate tables
are retained with download receipts and SHA-256 checksums. The source tables have
2024 modification dates; they are not the superseding 2026 individual-cell release.
Only public aggregate differential-expression tables were downloaded.

The initial completed evaluation produced **four usable signatures and two
insufficient signatures**. Both positive agreement and opposing directions occur;
missing coverage and exact neutral scores remain separate categories. Each state
uses its own FDR family. Searching across states does not supply joint significance.

`data/replication/protocol.json` was frozen after file inventory/header inspection
and before examining external gene coefficients. It pins the baseline run, scoring
code, input files, 2,000 permutations, seed 0 and minimum overlap 3. Published
progression p-values are BH-adjusted over each table's finite p-value family;
q <= 0.05 genes are ordered by q, absolute coefficient and publisher row ID, then
capped at 100 HGNC-resolved genes per direction. The unnamed CSV index contains
gene symbols; the publisher's `gene_id` column is a numeric row identifier.
The adapter's implementation receipt is separate from the frozen statistical protocol.
Sparse states are recorded without lowering thresholds. Frozen code/input drift
stops a rerun; amendments require a distinct protocol and evaluation, not editing
the original freeze to make a changed analysis appear prespecified.

```bash
cd /path/to/forge-workspace/microglia_repurposing
.venv/bin/python fetch_replication_sources.py  # public aggregate downloads, cached thereafter
.venv/bin/python run_replication.py            # offline evaluation of every frozen state
```

Signatures, full analysis runs, implementation receipts and complete reports are
stored by content hash under `data/replication/`. The report pointer is published
only after all declared states finish. A rerun does not overwrite the original
report or baseline scores. The dashboard compares exact perturbation IDs across
the full union, with nulls for missing coverage.

**Interpretation:** SEA-AD recruitment (ACT/UW ADRC) differs from ROSMAP, but a
donor-level crosswalk is unavailable locally. Donor non-overlap is unverified.
Brain region, Micro-PVM state and continuous progression contrast also differ
from the original Mathys comparison. These results are external cohort transfer
evaluations, not established independent replication, experimental rescue or
clinical validation. All SEA-AD states share a cohort and the LINCS panel.

The benchmark workspace contains six legacy compounds, three with linked primary
evidence dossiers. Clinical outcomes and microglia depletion findings are distinct
from rescue-assay labels. **Zero rescue controls are human-adjudicated; rescue
AUROC remains null.** Existing compounds have already been examined and are not
a blinded holdout. Proposed assay controls are clearly labeled as unobserved.

Source records: [SEA-AD resources](https://brain-map.org/consortia/sea-ad/our-resources),
[public aggregate files](https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html),
[ROSMAP recruitment](https://www.radc.rush.edu/molecular_networks/datasets.html).
The next scientific step is donor/contrast adjudication and reference-backed
rescue assays in the declared cell state; a positive overlap alone is insufficient.


## DAM rescue revision 2: evidence requirements

The active DAM workspace reads `data/rescue/v2/`. The original four v1 JSON
artifacts remain unchanged. The legacy freeze script now refuses to overwrite
an existing freeze or observation log. Run `.venv/bin/python run_rescue_panel.py`
to regenerate the active report; the dashboard reads the same evaluator.

This revision is an engineering contract awaiting assay qualification. The
Keren-Shaul citation is corrected to https://pubmed.ncbi.nlm.nih.gov/28602351/
(DOI 10.1016/j.cell.2017.05.018). Each proposed gene now carries an explicit
pending source-locator, direction, species, and human-transfer review. Correcting
the bibliography does not validate the manually assembled human signature.

A proposed functional endpoint, amyloid clearance per viable cell, is required
alongside signature reversal, homeostatic restoration, and quantitative viability.
Aβ and LPS/IFN-γ have separate induction specifications. Methods, control-derived
thresholds, sample-size justification, dose/time plan, and adjudication plan are
still unset. Freeze these in a new protocol revision before confirmatory work.
The software minimum of two donors and two replicates is only a safeguard against
single-sample claims, not a scientifically justified sample-size calculation.

Observation intake retains incomplete records for review, but incomplete records
cannot qualify. It records donor, replicate, batch, arm, dose, time, induction
shift, viability-loss fraction, signature score, homeostatic restoration, and
functional effect. Receipt files must reside under `data/rescue/` and match their
SHA-256. Adjudication requires a reviewer, evidence receipt, decision, and the
`observation_identity(record)` hash binding that decision to the exact record.
Hashes provide integrity and linkage, not proof of experimental truth or reviewer
identity; independent scientific review remains necessary.

All observations contribute to eligibility. A later success cannot hide an earlier
failure. Conflicting outcomes and multiple conditions stay blocked pending a
separate condition-specific analysis. Duplicate observation IDs, wrong cell models,
nonfinite measurements, and malformed resealed logs are rejected. Negative labels
cannot contradict a complete set of passing rescue measurements.

Prediction AUROC is deliberately unavailable. A separate prospectively registered
prediction table and independent outcome adjudication must be implemented before
reporting predictive performance. Assay reversal scores are never used as a
substitute, and development and holdout compounds are not pooled.

### Next scientific work

1. Review the 22 gene directions against primary source tables and human evidence.
2. Qualify each induction arm and functional readout; define quantitative thresholds
   and a donor/replicate analysis plan from assay performance and power requirements.
3. Freeze a new qualified protocol and predictions before confirmatory observations.
4. Collect blinded, source-backed observations and independently adjudicate outcomes.

The neighboring `../osteoarthritis_repurposing` workbench has separate cartilage,
synovium, and chondrocyte-rescue contracts. Its evidence and hashes have not been
changed by this upgrade. Its rescue evaluator needs a separate corresponding
review; DAM qualification does not establish cartilage rescue.
