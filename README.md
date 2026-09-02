# Country-Wide Seagrass Above-Ground Carbon Estimation

R and Google Earth Engine scripts for:

Hafizt, M., Roelfsema, C., Phinn, S., Lyons, M., Wicaksono, P., Salsabila, H.N.,
McMahon, K., Hernawan, U.E., Irawan, A., Adi, N.S., & Rahmawati, S. (in
preparation). **From Area-Based to Pixel-Based: A Country-Wide Approach for
Seagrass Above-Ground Carbon Estimation in Tropical Coastal Ecosystems**.

Companion repositories for the earlier chapters in this thesis:
- Chapter 4, Study 1 (GSE spectral percent-cover):
  [seagrass-spc-gse-indonesia](https://github.com/Emha89/seagrass-spc-gse-indonesia)
- Chapter 4, Study 2 (proxy-based AGC, 6 regions):
  [seagrass-agc-proxy-indonesia](https://github.com/Emha89/seagrass-agc-proxy-indonesia)

## Overview

This repository extends the proxy-based AGC framework from Chapter 4 (Study
2) to country-wide coverage across Indonesia, and compares two ways of
estimating total AGC: an **area-based** approach (persistence area x mean
density) and a **pixel-based** approach (summing per-pixel model
predictions). Uncertainty is quantified using a cluster design effect (DEFF)
/ intraclass correlation (ICC) framework, reflecting the compiled field data's
survey-cluster structure -- a different approach from Study 2, which retains
a simpler patch-size-based correction (L=30m) at regional scale.

STATUS: this repository is being assembled script by script, in the same
review process used for the Study 2 repository. This README will fill in as
scripts are reviewed -- sections below are placeholders until confirmed
against the actual code.

## Repository structure

```
seagrass-agc-countrywide-indonesia/
├── R/     -- to be filled in as scripts are reviewed
├── gee/   -- to be filled in as scripts are reviewed
└── data/  -- not included, see Data availability below
```

## Data availability

This repository contains **analysis code only**. Field survey data compiled
from third-party sources are subject to the data-sharing policies of the
original providers and are not included here (see the paper's Data
Availability statement once published).

## Open items

- [ ] Full script inventory received (R: 27 files, GEE: 24 + 7 under
      GEE_Apps/) -- content not yet reviewed. Several stages have
      multiple versioned variants to resolve as scripts come in:
  - R: `recompute_L.R` / `_v2.R` / `_v3.R` (3 versions)
  - GEE: `04_COVER_rf_modelDev_indo_v7/v8/v9`,
    `07_AGC_rf_modelDev_indo_v7/v8/v9`,
    `07_AGC_rf_calculationC_indo_v9/v10`,
    `08_AGC_rf_compare_v9/v10/v11`,
    `08_AGC_rf_modelUncer_indo{MPA,ROI}_v7/v8/v9` (MPA vs ROI may be
    different analysis domains, not just versions -- to confirm)
  - GEE_Apps/: several duplicate/backup/dated entries for the same app
- [ ] Add data availability / field-data template once column requirements
      are confirmed from the reviewed scripts

## Citation

_(to be added once the paper is accepted)_

## License

This project is licensed under the MIT License -- see the
[LICENSE](LICENSE) file for details. Field data used to produce the
results are not covered by this license -- see Data availability above.

## Contact

Muhammad Hafizt
School of the Environment, The University of Queensland, Brisbane, Australia
National Research and Innovation Agency of Indonesia (BRIN)
m.hafizt@uq.edu.au
