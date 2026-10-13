################################################################################
### National life tables and Arriaga decomposition, 2000-2019
###
### Inputs
###   data/ghedata.csv      WHO GHE 2021 deaths and population (JPN, KOR, SGP, TWN)
###   data/tbdcdata.csv     Taiwan death registry, harmonised to GHE causes
###                         (see S1 Appendix, Text A); iso3 = "TWN(TBDC)"
###   data/ghecauseMap.csv  GHE level-3 to level-2 cause mapping
###
### Outputs
###   data/LEdata.csv              abridged life tables, 2000-2021
###   data/LEdecomp_2000_2019.csv  age- and cause-specific contributions to the
###                                change in life expectancy, 2000-2019
###
### Age groups: 0 = 0-27 days, 0.1 = 1-11 months, 1 = 1-4 years,
###             5, 10, ..., 80 = five-year groups, 85 = 85+.
### "TWN" (WHO GHE for Taiwan) is kept for the sensitivity analysis only;
### the main analysis uses "TWN(TBDC)".
################################################################################

###################
### Environment ###
###################
pkgs <- c("data.table", "dplyr", "pbmcapply")
invisible(sapply(pkgs, require, character.only = TRUE))
rm(pkgs)

# Run from the repository root (all paths are relative to it).

############
### Data ###
############
ghecauseMap <- fread("data/ghecauseMap.csv")
data <- bind_rows(fread("data/ghedata.csv"), fread("data/tbdcdata.csv"))


#################
### Functions ###
#################

### Average years lived by those dying below age 5 (ax)
### Coale-Demeny relations (Preston et al. 2001, p. 48), with a0 rescaled to the
### neonatal (0-27 days) and post-neonatal (1-11 months) intervals.
youngExtension <- function(mx, sex){
  # mx: death rates for 0-27 days, 1-11 months, 1-4 years
  if(sex == 1){
    a0   <- function(m) ifelse(m >= 0.107, 0.330, 0.045 + 2.684 * m)
    a1_4 <- ifelse(mx[3] >= 0.107, 1.352, 1.651 - 2.816 * mx[3])
  } else {
    a0   <- function(m) ifelse(m >= 0.107, 0.350, 0.053 + 2.800 * m)
    a1_4 <- ifelse(mx[3] >= 0.107, 1.361, 1.522 - 1.518 * mx[3])
  }
  c(a0(mx[1]) * 28 / 365.24,
    a0(mx[2]) * (365.24 - 28) / 365.24,
    a1_4)
}


### Kannisto extension above age 70
### A logistic (Kannisto) model is fitted to single-year qx at ages 70-84 and
### extrapolated to age 129. Returns
###   ax : average years lived in 70-74, 75-79, 80-84 and 85+ by those dying
###   e85: life expectancy at 85 from the extrapolated single-year table,
###        used to close the open interval (NA if the fit fails)
ktExtension <- function(lx70){
  fallback <- list(ax = rep(2.5, 4), e85 = NA_real_)
  if(length(lx70) != 4 || any(!is.finite(lx70)) || any(lx70 <= 0)){
    warning("ktExtension: need positive lx at ages 70, 75, 80, 85; using ax = 2.5 and the crude 85+ closure.")
    return(fallback)
  }

  # Constant hazard within 70-74, 75-79 and 80-84
  mux <- -diff(log(lx70)) / 5

  # Logit of the single-year qx at the start of each group on the group midpoint
  qx1 <- 1 - exp(-mux)
  y <- log(qx1 / (1 - qx1))
  x <- c(72.5, 77.5, 82.5)
  ok <- is.finite(y)
  if(sum(ok) < 2){
    warning("ktExtension: fewer than two valid points for the logit fit; using ax = 2.5 and the crude 85+ closure.")
    return(fallback)
  }
  fit <- lm(y[ok] ~ x[ok], weights = x[ok])

  # Single-year life table from 85 to 129
  age85 <- 85:129
  qx85 <- plogis(coef(fit)[1] + coef(fit)[2] * (age85 + 0.5))
  lx85 <- lx70[4] * cumprod(c(1, 1 - qx85[-length(qx85)]))
  dx85 <- lx85 * qx85

  # Single-year deaths at 70-84 under the constant hazard
  lx7084 <- rep(lx70[1:3], each = 5) * exp(-rep(mux, each = 5) * rep(0:4, 3))
  dx7084 <- lx7084 - c(lx7084[-1], lx70[4])

  # ax per group, assuming deaths at the middle of each single year of age
  dx  <- c(dx7084, dx85)
  grp <- rep(c(70, 75, 80, 85), c(5, 5, 5, length(age85)))
  t   <- c(rep(0:4, 3), age85 - 85) + 0.5
  ax  <- tapply(dx * t, grp, sum) / tapply(dx, grp, sum)

  list(ax = unname(ax), e85 = sum(lx85 - dx85 / 2) / lx70[4])
}


### Abridged life table for one country-sex-year
LEcalc <- function(data, year_index, sex_index, iso3_index){

  lt <- data %>%
    filter(year == year_index & sex == sex_index & iso3 == iso3_index &
             causename == "All Causes") %>%
    arrange(age) %>%
    select(iso3, sex, year, causename, age, pop, dths)
  n <- nrow(lt)

  lt <- lt %>%
    mutate(
      nx = c(28 / 365.24, (365.24 - 28) / 365.24, 4, rep(5, n - 4), NA),
      ax = c(0.2, 0.330, 1.352, rep(2.5, n - 3)),   # starting values
      Mx = dths / pop,
      qx = ifelse(row_number() == n, 1, nx * Mx / (1 + (nx - ax) * Mx)),
      qx = ifelse(qx > 1, 0.999, qx),
      px = 1 - qx,
      lx = 1e5 * cumprod(lag(px, default = 1)),
      dx = lx * qx
    )

  # Refined ax below 5 and at 70+; qx above keeps the starting ax values
  lt$ax[lt$age < 5] <- youngExtension(lt$Mx[lt$age < 5], sex_index)
  kt <- ktExtension(lt$lx[lt$age >= 70])
  lt$ax[lt$age >= 70] <- kt$ax

  # Open interval closed with the Kannisto extrapolation (fallback: lx / Mx)
  L85 <- if(is.na(kt$e85)) lt$lx[n] / lt$Mx[n] else lt$lx[n] * kt$e85

  lt %>%
    mutate(
      Lx = ifelse(row_number() < n, nx * lead(lx) + ax * dx, L85),
      Tx = rev(cumsum(rev(Lx))),
      ex = ifelse(lx == 0, 0, Tx / lx)
    )
}


### Arriaga decomposition by age
### Population 1 is the reference (year 2000); positive values are gains.
### The indirect effect of the open interval is zero.
ArriagaDecomp_age <- function(data,
                              year_index1, sex_index1, iso3_index1,
                              year_index2, sex_index2, iso3_index2){
  lt1 <- LEcalc(data, year_index1, sex_index1, iso3_index1)
  lt2 <- LEcalc(data, year_index2, sex_index2, iso3_index2)

  DE <- lt1$lx * (lt2$Lx / lt2$lx - lt1$Lx / lt1$lx) / 1e5
  IE <- (lt1$lx * lead(lt2$lx) / lt2$lx - lead(lt1$lx)) / 1e5 * lead(lt2$ex)
  IE[nrow(lt1)] <- 0

  data.frame(
    iso3 = paste0(iso3_index2, "-", iso3_index1),
    sex = paste0(sex_index2, "-", sex_index1),
    year = paste0(year_index2, "-", year_index1),
    causename = "All Causes",
    age = lt1$age,
    contri = DE + IE
  )
}


### Arriaga decomposition by age and cause
### Each age contribution is split across causes in proportion to the change in
### cause-specific death rates (proportional allocation).
ArriagaDecomp_cause <- function(data,
                                year_index1, sex_index1, iso3_index1,
                                year_index2, sex_index2, iso3_index2){

  ageDecomp <- ArriagaDecomp_age(data,
                                 year_index1, sex_index1, iso3_index1,
                                 year_index2, sex_index2, iso3_index2)
  lt1 <- LEcalc(data, year_index1, sex_index1, iso3_index1)
  lt2 <- LEcalc(data, year_index2, sex_index2, iso3_index2)
  n <- nrow(lt1)

  # Open interval (85+): use the life-table rate l85 / L85 (= 1 / e85) so that
  # the allocation follows the Kannisto-closed table, not the crude rate
  rates <- data.frame(
    age = lt1$age,
    contri = ageDecomp$contri,
    M1 = c(lt1$Mx[-n], lt1$lx[n] / lt1$Lx[n]),
    M2 = c(lt2$Mx[-n], lt2$lx[n] / lt2$Lx[n])
  )

  # Share of each cause in all deaths, by age
  causeShare <- function(year_index, sex_index, iso3_index){
    data %>%
      filter(year == year_index & sex == sex_index & iso3 == iso3_index &
               ghecause %in% ghecauseMap$ghecause_lv3) %>%
      group_by(age) %>%
      mutate(prop = dths / sum(dths)) %>%
      ungroup() %>%
      select(age, causename, prop)
  }

  causeShare(year_index1, sex_index1, iso3_index1) %>%
    inner_join(causeShare(year_index2, sex_index2, iso3_index2),
               by = c("age", "causename"), suffix = c("1", "2")) %>%
    left_join(rates, by = "age") %>%
    mutate(
      P = ifelse(prop2 == prop1, prop2, (prop2 * M2 - prop1 * M1) / (M2 - M1)),
      contri = contri * P,
      iso3 = paste0(iso3_index2, "-", iso3_index1),
      year = paste0(year_index2, "-", year_index1),
      sex = paste0(sex_index2, "-", sex_index1)
    ) %>%
    arrange(factor(causename, levels = sort(unique(causename))), age) %>%
    select(iso3, year, sex, age, causename, contri)
}

########################
### End of Functions ###
########################


###################
### Life tables ###
###################
params <- expand.grid(
  year = 2000:2021,
  sex = 1:2,
  country = c("TWN", "JPN", "KOR", "SGP", "TWN(TBDC)"),
  stringsAsFactors = FALSE
)

LEdata <- bind_rows(pbmclapply(seq_len(nrow(params)), function(i){
  LEcalc(data, params$year[i], params$sex[i], params$country[i]) %>%
    select(iso3, sex, year, age, pop, dths, nx, ax, Mx, qx, px, lx, dx, Lx, Tx, ex)
}))
rm(params)
write.csv(LEdata, "data/LEdata.csv", row.names = FALSE)


####################################################
### Contributions to the change in LE, 2000-2019 ###
####################################################
countries <- c("JPN", "TWN", "KOR", "SGP", "TWN(TBDC)")

LEdecomp_2000_2019 <- bind_rows(lapply(1:2, function(s){
  bind_rows(lapply(countries, function(cty){
    ArriagaDecomp_cause(data, 2000, s, cty, 2019, s, cty)
  }))
}))
rm(countries)
write.csv(LEdecomp_2000_2019, "data/LEdecomp_2000_2019.csv", row.names = FALSE)
