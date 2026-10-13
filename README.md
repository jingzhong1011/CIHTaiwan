# Uneven progress in survival across high-income Asia: A comparative analysis of age, cause, and geographic inequalities

Code and data for:

> Wang J-Z, Nomura S, Chang AY, Lin H-H. Uneven progress in survival across high-income Asia: A comparative analysis of age, cause, and geographic inequalities. *(under review)*

The analysis decomposes changes in life expectancy at birth between 2000 and 2019 in Japan, the Republic of Korea, Singapore, and Taiwan by age and cause of death (Arriaga decomposition). It also compares 2019 subnational life expectancy across Taiwan (counties), Japan (prefectures), and Korea (provinces).

## Repository structure

```
LEdecomp_national.R       National abridged life tables and Arriaga decomposition (age x cause)
LEest_subnational.R       Subnational life expectancy (Taiwan counties, Japan prefectures, Korea provinces)
lifeTable_subnational.R   Life-table functions for subnational estimates (sourced by LEest_subnational.R)
data/                     Inputs and outputs (see below)
```

All scripts are run from the repository root. Analyses were done in R 4.6.0.

## Data

### Included inputs

| File | Description | Source |
|---|---|---|
| `data/ghedata.csv` | Deaths and population by cause, age, sex and year, 2000-2021, for JPN, KOR, SGP and TWN | WHO Global Health Estimates 2021 |
| `data/tbdcdata.csv` | Taiwan deaths harmonised to GHE causes, 2000-2021; garbage codes redistributed as described in S1 Appendix, Text A | NHIRD, Taiwan |
| `data/ghecauseMap.csv` | GHE level-3 to level-2 cause mapping | WHO GHE cause list |
| `data/causeMapping.csv` | Taiwan cause-of-death categories mapped to GHE causes (S1 Appendix, Table C) | Authors |
| `data/JPNdata/DeathCount_japan2019.csv` | Deaths and population by prefecture, age and sex, 2019 | Excess and Exiguous Deaths in Japan dashboard (https://exdeaths-japan.org/en/) |
| `data/KORdata/*.csv` | Deaths and resident population by province, age and sex, 2019 | Korean Statistical Information Service (https://kosis.kr/eng/) |

### Included outputs

| File | Description |
|---|---|
| `data/LEdata.csv` | Abridged life tables by country, sex and year, 2000-2021 (`TWN(TBDC)` = Taiwan national registry, main analysis; `TWN` = WHO GHE for Taiwan, sensitivity analysis) |
| `data/LEdecomp_2000_2019.csv` | Age- and cause-specific contributions (years) to the change in life expectancy, 2000-2019 |
| `data/e0_taiwanCounty.csv` | Life expectancy at birth by Taiwan county, 2000 and 2019 |
| `data/e0_japanPrefecture2019.csv` | Life expectancy at birth by Japanese prefecture, 2019 |
| `data/e0_koreaProvince2019.csv` | Life expectancy at birth by Korean province, 2019 |

### Restricted data (not included)

- `data/restricted/TaiwanTownship/`: township-level deaths and population used by `LEest_subnational.R` for Taiwan counties.

Researchers who meet the criteria for access can request these data from the Health and Welfare Data Science Center, Ministry of Health and Welfare, Taiwan.

## How to run

1. `LEdecomp_national.R` -> `data/LEdata.csv`, `data/LEdecomp_2000_2019.csv`
2. `LEest_subnational.R` (the Taiwan section requires restricted data; the Japan and Korea sections use the included inputs) -> `data/e0_*.csv`

The data underlying each figure are in the included outputs: Fig 1 (`LEdata.csv`), Figs 2-4 and S1 Appendix Figs A-B (`LEdecomp_2000_2019.csv`), and Fig 5 (`e0_*.csv`). S1 Appendix Figs C-D compare `ghedata.csv` with `tbdcdata.csv`.

R packages: data.table, dplyr, stringr, pbmcapply.

## Licence

The code is released under the MIT License (see `LICENSE`). WHO Global Health Estimates are redistributed under WHO's terms of use, with attribution to the World Health Organization.
