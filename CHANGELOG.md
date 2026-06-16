# Changelog

All notable changes to this project should be documented in this file.

The project uses [Semantic Versioning](https://semver.org/): `MAJOR.MINOR.PATCH`.

## [Unreleased]

### Added
### Changed
### Fixed

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
