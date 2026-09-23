# ABCD Data-Access Notes

The ABCD application in the manuscript uses restricted subject-level data from
the ABCD Study. Raw ABCD files, subject-level derived matrices, and local
intermediate workspaces are not distributed in this public repository.

The public repository is therefore package-centered: it provides the GG-ComBat
implementation, toy examples, tests, and simulation scripts. Authorized ABCD
users can reproduce the application only after obtaining the relevant ABCD data
release through the approved data-access process and reconstructing the private
input files in their own secure environment.

Examples of private inputs used in the authors' local analysis include imaging
feature tables, scanner/site metadata, demographic covariates, and derived
harmonized matrices. These files are intentionally excluded from Git tracking.

This design follows the data-use boundary: method code and simulation settings
are public, while restricted subject-level ABCD data and private analysis
artifacts remain outside the repository.
