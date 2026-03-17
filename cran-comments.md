## Response to CRAN reviewer (Uwe Ligges, 3rd submission 2026-03-17)

* Fixed: Removed the LICENSE.md file URI from README.md entirely.
  LICENSE.md is excluded from the package tarball via .Rbuildignore,
  which caused the invalid file URI note in previous submissions.
  The license is already declared as GPL (>= 3) in DESCRIPTION.
  No file URI link is needed or included.

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
