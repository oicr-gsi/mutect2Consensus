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

## 1.0.3 - 2024-06-28
- Added tumour and normal maf clean up before filtering (convert to numeric, replace NaN with 0)
- Changed some output file names to make them clear 
- Removed some ouptput files to simplify
- Updated subworkflows to new version
- Added vidarr labels

## 1.0.0 -2023-10-02
- Created new workflow
