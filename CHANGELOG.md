# Changelog

All notable changes to this project should be documented in this file.

The project uses [Semantic Versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

## [Unreleased]

### Added
### Changed
### Fixed

## [2.1.0] - 2026-07-16

### Added
- Added support for retrieving and annotating in-frame deletions, insertions, and deletion-insertion variants from cBioPortal, COSMIC, and ClinVar.
- Added a variant_types argument for selecting the protein variant types to retrieve (default: missense).
- Added parsing and internal representations for supported protein and genomic variant types.
- Added genomic annotation and assembly-conversion support for indels.
- Extended gnomAD annotation to supported insertions, deletions, and deletion-insertion variants.
- Added an indel retrieval tutorial and corresponding example output.

## Changed
- Restricted REVEL annotation to missense variants.
- Updated the documentation to describe indel retrieval and the variant_types argument.

## [2.0.1] - 2026-07-15

### Fixed

- Changed ClinVar data source class to filter out genomic variants using unsupported genome assemblies

## [2.0.0] - 2026-06-16

### Added

- Added a documented release process based on Semantic Versioning and annotated git tags
- Added a Markdown changelog with an `Unreleased` section for collecting upcoming changes
- Overhaul of inner workings of the package to support non-missense mutations (PR #265 #270)

### Changed

- Centralized package versioning so `setup.py` reads the version from `cancermuts.__version__`.

## [0.1.0] - 2018-07-02

### Added

- First preliminary release of cancermuts.
- Logging system and associated partial refactoring.
