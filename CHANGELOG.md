# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.6] - 2025-03-12
## Changed
- Make normal inputs optional, allow tumor only mode: https://jira.oicr.on.ca/browse/GRD-911

## [1.0.5] - 2024-12-11
### Added
- provision out unfiltered maf
- Boolean for whether filter maf

### Changed 
- output file name changes
- task name changes
- meta changes

### Removed
- filtering for a panel of genes

## [1.0.4] - 2024-11-14
### Changed
- Update maf-filtering module version.

## [1.0.3] - 2024-06-28
### Added
- Added tumour and normal maf clean up before filtering (convert to numeric, replace NaN with 0)
- Added vidarr labels

### Changed
- Changed some output file names to make them clear 
- Removed some ouptput files to simplify
- Updated subworkflows to new version

## [1.0.2] - 2023-11-17
### Fixed
- Remove extra peroid in file extension

## [1.0.1] - 2023-10-27
### Added
- Release a version atfer updated filter maf.

## [1.0.0] - 2023-10-02
### Added
- Created new workflow