# ninetails 1.0.84 (development)

## cDNA pipeline: independent validation of orientation calls

The cDNA pipeline classifies each read as `polyA`, `polyT` or `unidentified` by
matching Dorado's SSP and VNP primers against both ends of the basecalled
sequence. That call is made entirely from sequence, so nothing in the pipeline
previously contradicted it when it was wrong.

* **`infer_cdna_layout()`** reads the orientation off the *signal* instead, by
  comparing the sizes of the pre-tail and post-tail regions and taking the larger
  side to be the transcript body. It is deliberately strict — a region below
  `min_region_samples` (default 100), or a larger/smaller ratio below
  `strictness_ratio` (default 1.5), returns `"ambiguous"` rather than committing
  to an orientation.
* **`cdna_layout_agreement_marker()`** pairs the sequence-based call with the
  layout-based one and collapses the result to a single character, so
  disagreements are visible at a glance in a table of thousands of reads.
* **`launch_cdna_signal_browser()`** puts the two side by side. It takes a Dorado
  summary and the output of `detect_orientation_multiple()`, joins them on
  `read_id`, and reports how many reads carried no matching `tail_type` rather
  than silently dropping them.

The two methods are independent by construction — one reads primers, the other
reads region geometry — so agreement is evidence and disagreement is a read worth
looking at.

## Dorado 1.4.0 compatibility

Dorado 1.4.0 renamed the POD5 column in its summary output from `filename` to
`input_filename`. `preprocess_inputs()` and `launch_cdna_signal_browser()` now
accept either spelling, normalising to `filename` internally.

## Bug fixes

* **Read IDs containing underscores are no longer truncated.** Chunk names have
  the form `<read_id>_<chunk_index>`, and the read ID was recovered by splitting
  on `_`, which cut UUID-style identifiers short and orphaned their predictions.
  Only the trailing `_<digits>` is stripped now.
* **`detect_orientation_multiple()` no longer overwrites its input files.** The
  classified table is written to a sibling `<name>_classified.tsv`, leaving the
  original sequence files untouched. Overwriting them risked corrupting user
  input on a partial failure, and made re-runs non-idempotent — the second pass
  would re-classify an already-tagged file.
* **An empty signal vector no longer aborts feature extraction.**
  `filter_signal_by_threshold()` returns no pseudomoves for a zero-length signal
  rather than iterating its z-score loop backwards.
* `split_tail_centered_dorado()` pads chunks that would start before the
  beginning of the signal with random draws from the five most frequent observed
  values, rather than a fixed constant, keeping chunk length at 100 without
  introducing a flat artificial segment for the CNN to learn from.

## Documentation

* The package now has a changelog, rendered on the pkgdown site.
* The pkgdown reference index groups all exported functions by pipeline and role,
  with the Guppy pipeline marked legacy throughout.
* Vignettes covering detection, postprocessing, plotting, signal inspection, the
  Shiny dashboard and tailfindr compatibility.

---

# ninetails 1.0.78

*Released 2026-05-07.*

Dorado-compatible DRS processing with a Shiny application for signal browsing.

* `launch_signal_browser()` for Dorado DRS (POD5) data and
  `launch_signal_browser_guppy()` for Guppy legacy (fast5) data, with non-A
  modification overlays on the raw signal trace.
* The cDNA pipeline is present but **under construction** and should not be used
  for analysis.

---

# ninetails 1.0.4

*Released 2025-09-28.*

* **Full POD5 and Dorado support.** `check_tails_dorado_DRS()` processes direct
  RNA sequencing data basecalled with Dorado ≥ 1.0.0, reading poly(A) coordinates
  from the Dorado summary rather than Nanopolish.
* **The fast5/Guppy pipeline is preserved as legacy mode** via
  `check_tails_guppy()`, and is no longer actively developed.
* Additional pipelines announced as forthcoming.

---

# ninetails 1.0.3

*Released 2024-12-05.*

* Fixed the export of the tailfindr compatibility function.

---

# ninetails 1.0.2

*Released 2024-04-02.*

* **Compatibility with tailfindr output.** `check_polya_length_filetype()`
  detects whether a poly(A) length file came from Nanopolish or tailfindr and
  converts the latter via `convert_tailfindr_output()`. Only DRS data are
  accepted; tailfindr cDNA output is rejected explicitly.

The `v.1.0.2_manuscript` tag (2024-08-13) marks the version accompanying the
manuscript.

---

# ninetails 1.0.0

*Released 2023-11-23.*

* First stable release. Features added, code optimised, typos fixed.

---

# ninetails 0.9.0

*Released 2023-10-03.*

* Performance fixes.

---

# ninetails 0.7.0

*Released 2023-01-18.*

* Performance fixes.
* Fixed potential issues arising from input file incompatibility.
* Clarified comments and descriptions throughout.

---

# ninetails 0.4.1

*Released 2022-10-13.*

* Added statistical functions, data visualisation, and postprocessing features.
* Fixed a bug that diverted `stdout`.

---

# ninetails 0.4.0

*Tagged 2022; no accompanying release notes.*

---

# ninetails 0.3.1

*Released 2022-09-02.*

* Windows compatibility.

---

# ninetails 0.3.0

*Released 2022-08-23.*

* Fixed R dependency issues.
* Fixed the Nanopolish input argument.

---

# ninetails 0.2.0

*Released 2022-08-23.*

* Prerelease with major updates.

---

*Entries for 1.0.78 and earlier were reconstructed from the published GitHub
release notes at <https://github.com/LRB-IIMCB/ninetails/releases> and are
necessarily coarser than the development section above.*
