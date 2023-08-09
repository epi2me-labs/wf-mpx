# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Changed
- Bumped minimum required Nextflow version to 23.04.2.
- Enum choices are enumerated in the `--help` output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.
- Updated GitHub issue templates to force capture of more information.
- Example command to use demo data.
- Spaces in sample aliases are now replaced with underscores.
- Deprecated parameter `process_label`.

## [v0.0.7]
### Changed
- Updated whole workflow to bring up-to-date with recent template changes
### Added
- Configuration for running demo data in AWS

## [v0.0.6]
### Changed
- Updated description in manifest

## [v0.0.5]
### Removed
- conda support
### Changed
- Updating templates
- Documentation
### Fixed
- Fix for error with empty VCF files (@MarieLataretu)

## [v0.0.4]
### Added
- Minimum coverage parameter

## [v0.0.3]
### Changed
- Args parser for fastqingress
- Set out_dir option type to ensure output is written to correct directory on Windows

## [v0.0.2]
### Changed
- Tidied up docs
- Restricting assembly to only reads that map to mpx

## [v0.0.1]

First release
