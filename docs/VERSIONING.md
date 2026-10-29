<!-- SPDX-License-Identifier: MIT -->
# Versioning

## The governing rule

**A link printed in a paper must resolve to the exact atlas version that paper used,
forever.** A reader following that citation in 2030 must not silently receive a newer
model with different numbers.

Everything below exists to enforce that one sentence.

## How a release is laid out

```
data/releases.json            stable path, never versioned, lists every release
data/v0.1/                    a release directory -- never modified after publication
  mra_web_model.json          the fitted model, with `version` stamped inside it
  nt_level_supply.csv         published nucleotide rates, positions 1-8213
  aa_level_supply.csv         published amino-acid mutabilities
  aa_paths/<Protein>.csv      codon-to-codon path shards, loaded lazily
  nl43_reference.json         reference sequence, protein assignment, codon frames
  hxb2_align.json             HXB2 lookup (Env and Vpu only; identity elsewhere)
  citation.json               machine-readable citation block
  manifest.json               SHA-256 and gzipped size per asset
data/v1.0/                    the next release: a NEW directory, nothing overwritten
```

## The five mechanisms

1. **Versioned asset paths.** Nothing is ever overwritten in place. A new model is a
   new directory. `scripts/build_web_assets.py --version vX.Y` refuses to write into a
   directory that already contains a manifest.

2. **URL-addressable version.** `?v=0.1` selects a release; absent, the portal serves
   `latest`. An unparseable or unknown version is a loud error, never a silent fallback
   to latest -- a typo in a printed URL must fail visibly. Every shareable link the UI
   generates carries an explicit `?v=`, never the implicit default.

3. **Version stamped inside the bundle.** `extract_mra_web_model.py --version` writes
   `version`, `release_date` and `doi` into `mra_web_model.json`, so a downloaded file
   is self-identifying. `loadModelBundle()` in `js/version.js` compares the embedded
   string against the directory it was fetched from and **throws** if they disagree:
   that mismatch means assets were copied between release directories, which would
   silently serve the wrong model.

4. **Provenance in every export.** Browse downloads and custom atlases carry the atlas
   version, the model name and the retrieval date in a header comment. A user who
   downloaded a custom atlas six months ago must be able to say which model produced it.

5. **Not-latest banner.** Viewing any release other than `latest` shows a persistent
   banner naming the current release, what changed, and a link to it -- while leaving
   the old version fully functional.

## Cutting a new release

```bash
conda activate webportal
python scripts/build_web_assets.py --version v1.0 --doi 10.xxxx/yyyyy \
    --src /path/to/data/models/mra/production
node test/roundtrip_nl43.mjs --version v1.0      # must pass before release
python test/compare_roundtrip.py --version v1.0
```

The build script regenerates `data/releases.json` from the directories that exist,
preserving every earlier entry. It will not touch `data/v0.1/`.

**v0.1 does not go to anyone outside the lab until the NL4-3 round-trip test passes.**
See `docs/validation_report.md`.

## Current status of v0.1

- `doi`: **not yet assigned.** Fill in via `--doi` when minted; until then exports say
  `[not yet assigned]` rather than inventing an identifier.
- `release_date`: stamped by the build script at build time.
- v0.1 is a pre-publication build. The intent is to re-stamp the same posterior as
  `v1.0` at paper release and leave `v0.1` in place, so the version machinery is
  exercised with two real directories before anything is cited.
