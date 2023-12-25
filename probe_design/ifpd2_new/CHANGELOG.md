# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased

## [1.0.1a1] - 2021-12-07
### Added
- New entry points: `db make`, `db dump`, `db info`, `db merge`, and `db_reindex`.
- Module `io` with parsing functions for output generated via oligo-melting, hush, or OligoArrayAux.
- Assert for compatibility with a given dtype. (only max value currently checked)

### Changed
- New database format.
- `dbchk` entry point changed to `db check`.
- Dropped sequence length column, as it can be inferred from the end/start columns.
- Removed `-X` default from `query`.

## [1.0.0-alpha] - 2021-02-26

[1.0.1a1] https://github.com/ggirelli/ifpd2/releases/tag/v1.0.1a1  
[1.0.0-alpha] https://github.com/ggirelli/ifpd2/releases/tag/v1.0.0-alpha  
