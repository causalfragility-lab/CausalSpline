## Response to CRAN reviewer (Uwe Ligges, 2026-03-17)
* Fixed: README.md now links to LICENSE.md (was: LICENSE)
  LICENSE.md is included in the package

## R CMD check results
0 errors | 0 warnings | 1 note

## Notes
* "New submission" note is expected for first CRAN submission
* "Unable to verify current time" is a Windows network/firewall 
  issue unrelated to the package itself
* Flagged words (Hirano, Imbens, IPW, ignorability,
  unconfoundedness) are standard causal inference terminology

## Test environments
* Windows 11 local: R 4.5.1 (2025-06-13 ucrt)
* win-builder: R-devel Windows Server 2022: 1 NOTE only
* GitHub Actions Ubuntu latest (release): OK
* GitHub Actions Ubuntu latest (devel): OK
* GitHub Actions Ubuntu latest (oldrel-1): OK
* GitHub Actions macOS latest (release): OK
* GitHub Actions Windows latest (release): OK

## Downstream dependencies
* None — this is a new package
