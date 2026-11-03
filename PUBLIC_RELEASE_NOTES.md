# Public release notes

Prepared 2026-11-03. This repository is a **public copy** of an author-maintained research repository. Its Git history, other branches and pull requests are not carried over. This file records every way the public copy differs from the original so that no recorded hash or piece of evidence is silently inconsistent.

## 1. Not included

* The Git history of the original repository and its branches, including assistant-authored audit branches and a pull request.
* Eight patient-level audit tables written by `tests/s2prep_metadata_audit.R` (`evidence/s2prep_audit_*.csv`). They are derived from the source dataset's event metadata (identifiers, episode timing, source-row indices) and are not needed for any result described here. Anyone with the source data who runs that script would regenerate them locally; they are not part of this release.
* Pseudonymous patient identifiers with episode numbers and times: two records in `evidence/s2prep_metadata_audit.json` were removed (section 3).
* The source dataset itself, its concentration values, and the protocol preflight file. None has ever been part of this repository.
* Local environment directories (`renv/library`).

## 2. Local-path sanitisation

Four files named a path on the author's machine. They were replaced as below. Behaviour is identical when the environment variables are set; with no variable set, the functions look for a bare file name in the working directory. The original files are preserved privately by the author.

| File | Change | SHA-256 before | SHA-256 after |
|---|---|---|---|
| `R/pkpd02_source_reader.R` | default workbook path replaced by `Sys.getenv("PKPD02_SOURCE_XLSX", "Amikacin_Data_File.xlsx")` | `0a403b42d215ca8d4f7c5636f2cba7e4109eec12b7a271fb77e4caf18639f920` | `d5bc8ec1edd0a9009e05d68b6a0d93049e07d382a6054159a39e677b06ef44ba` |
| `tests/s2prep_metadata_audit.R` | preflight-file path replaced by `Sys.getenv("PKPD02_PREFLIGHT_JSON", "PKPD02_v11_metadata_preflight.json")` | `e77a67e3bd3d54ae80cfcca5a10bec3ce50675a433ac2a88f2eeeaa96f48f99a` | `019bb224c1c3d2947255641e5cd11392271ed907a74114c41947a5f4bbe8c865` |
| `reports/S1_verification_report.md` | a local folder path in one sentence replaced by a placeholder | `6a7777855f2f34c28a47039a007e8b7f2f7dc7033dddd1a85e334957c8588774` | `ccbcc14d40166db777a6a2438a162111ef55ffa8f631118cf176489dcfb9a91b` |
| `archive/s1_failed_2026-10-30/reports/S1_verification_report.md` | same sentence, in the archived copy (the two copies remain byte-identical) | `6a7777855f2f34c28a47039a007e8b7f2f7dc7033dddd1a85e334957c8588774` | `ccbcc14d40166db777a6a2438a162111ef55ffa8f631118cf176489dcfb9a91b` |

## 3. Redaction of `evidence/s2prep_metadata_audit.json`

The audit script wrote a list of two "sample and dose at the same time" records, each holding a pseudonymous patient ID, an episode number, a time and row-order detail from the source dataset. That detail is not needed for any result in this repository, so the public copy replaces the list with an object that keeps only two aggregates computed from it.

| File | SHA-256 as the audit script wrote it | SHA-256 in this release |
|---|---|---|
| `evidence/s2prep_metadata_audit.json` | `c90b30a1d369113fec3b9d771f9b4a9c800e338d03040485e42d3e5847f6a08f` | `dc81e3eb80e0717e86779cb553cecb2bde8b1f53013db0e049bbe98981c335cd` |

* **Removed:** the per-record fields `ID`, `EPISODE`, `TIME`, `evids_in_file_order`, `TAD_of_sample`, `file_order_sample_first`, `TAD_gt0_but_dose_listed_first` (2 records).
* **Retained:** every other value, byte-for-byte (the tool asserts this), plus the aggregates *all records list the sample before the dose in file order* = true and *any record with TAD > 0 but the dose listed first* = false.
* **Not a regenerated result.** The file carries a top-level `public_release_redaction` object stating that it is a redacted copy and giving the hash of the original. The original is preserved privately by the author, and re-running `tests/s2prep_metadata_audit.R` with the source workbook recreates it.
* **Nothing depends on the removed field.** `tests/s2prep_run_all.R` reads only `audit_trail$DV_column_loaded` from this file, and the file's own hash is recorded nowhere else in the release (checked when the release was built).

## 4. Consequences for recorded hashes and evidence

* **S-1 archive manifest.** `archive/s1_failed_2026-10-30/MANIFEST.sha256` has one line changed (the S-1 report). Manifest SHA-256 before `67f940ad0af5595fc29f284cbcf7464253c648427005a459bf5110de54926cf8`, after `10a27cbd7ea53873733ae0af4cc9a58fa77ae3de81fb12f62973ffa93aa92759`. `tests/s1r_run_all.R` verifies the archive and the in-place originals against this manifest; it was re-run on this copy and both checks pass. **Nothing else in the failed-S-1 record changed.**
* **Quoted manifest hash.** `protocol/PKPD02_v12_amendment_DRAFT.md` and its archived copy quote the *pre-sanitisation* manifest hash. They were deliberately not edited (they are dated historical text and the archived copy is itself hash-manifested).
* **`evidence/s1r_results.json` was regenerated** by running S-1r on the sanitised tree, because it records the manifest hash and the hash of `R/pkpd02_source_reader.R`. SHA-256 before `a2866e1d762490e0c2d863dff9fdfea1536db7bb5c5e78175218bcdfd7dcfb5d`, after `4eb0ed25338ee34d3617c47ed63c012712cd9b76bd1e19545b480135c8da51c5`. Its numerical content is unchanged; it differs in `timestamp_utc`, the recorded manifest hash and the recorded reader hash. The S-1r point-level CSV evidence reproduced byte-for-byte.
* **Evidence that still records pre-sanitisation hashes.** These files were produced by scripts that read the source workbook, which cannot be re-run here without the data, so they were **not** regenerated. They remain an accurate record of the run on the original files:

| Evidence file | Records the hash of | Recorded value | Public-copy value |
|---|---|---|---|
| `evidence/s2prep_results.json` | `R/pkpd02_source_reader.R` | `0a403b42d215ca8d…` (pre-sanitisation) | `d5bc8ec1edd0a900…` |
| `evidence/s2prep_metadata_audit.json` | `tests/s2prep_metadata_audit.R` | `e77a67e3bd3d54ae…` (pre-sanitisation) | `019bb224c1c3d294…` |
| `evidence/s2prep_results.json` | `tests/s2prep_metadata_audit.R` | `e77a67e3bd3d54ae…` (pre-sanitisation) | `019bb224c1c3d294…` |

  The public copies of those two files differ from the originals only in the removed path (section 2).
* **No scientific method, threshold, protocol text or result was changed.**

## 5. Data source, licence and remaining rights uncertainty

* **Dataset.** Thomson A, Robb F, Siebinga H (2021), *Amikacin data for population pharmacokinetic modelling*, University of Strathclyde, https://doi.org/10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851. The publisher's page (checked 2026-11-03) labels the data file **CC BY 4.0** and describes it as anonymised patient data extracted from routine clinical files.
* **CC BY 4.0 section 3(a) items.** Creators credited above; the page's site-wide copyright line is "Copyright © 2026 University of Strathclyde, its licensors, and contributors"; licence link https://creativecommons.org/licenses/by/4.0/; **modifications:** no copy of the data is distributed, and the retained summaries (counts and structure statistics) were computed from the file's event and covariate columns, with the concentration column skipped at read time; the authors and publisher do not endorse this work.
* **The additional wording.** The same footer says "All rights are reserved, including those for text and data mining, AI training, and similar technologies," and then "For all open access content, the relevant licensing terms apply." On the evidence checked, this is a site-wide notice that defers to the licence on open-access content, and the data file is labelled CC BY 4.0, whose grant is irrevocable and covers sharing and adapting. No separate data-use policy was found, and the publisher has not been asked to confirm how the notice applies to the CC BY file.
* **What this release therefore contains from the dataset:** no rows, no concentration values and no per-patient values; only aggregate structural counts and the pass/fail outcomes of tests run with synthetic parameters.
* **Remaining uncertainty and narrowest remedy.** If the notice were read as restricting use of the CC BY file, the affected material is the three real-structure evidence files (`evidence/s2prep_metadata_audit.json`, `evidence/s2prep_real_structure_domain_check.json`, `evidence/s2prep_results.json`) and the README rows that cite them. Removing those files would not change any synthetic result.

## 6. Authorship and provenance

The original working repository's history contains commits by the author and by an AI agent (Devin); none of it is carried over. This repository starts from a single commit. Authorship, AI-assistance and independently reproducible evidence are described in the README section *Authorship, AI assistance and provenance*; the notes here do not add claims about who wrote what.

## 7. Release integrity

* `RELEASE_MANIFEST.sha256` lists the SHA-256 of every file except itself. Verify with `sha256sum -c RELEASE_MANIFEST.sha256` (Linux) or `shasum -a 256 -c RELEASE_MANIFEST.sha256` (macOS).
* `.gitattributes` contains `* -text`, so Git never converts line endings. `tests/s1r_run_all.R` verifies archived files by hash, and a CRLF conversion on Windows would otherwise make that check fail for a reason unrelated to the science.

## 8. Non-destructive verification

`tools/verify_public.R` (new in the public copy) runs each gate in its own temporary copy of the repository. It exists because running the S-1 script rewrites the in-place S-1 evidence and S-1r then reports a failure that is only an artefact of run order. The S-1 script, its evidence and the archive are unchanged. See README §4.
