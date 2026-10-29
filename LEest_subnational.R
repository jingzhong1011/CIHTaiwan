###################
### Environment ###
###################
pkgs <- c("data.table", "dplyr", "stringr")
invisible(sapply(pkgs, require, character.only = TRUE))
rm(pkgs)
source("lifeTable_subnational.R")


############
### Data ###
############
# data format: iso3, year, sex, area, age

# Taiwan township-level deaths and population are restricted (see README) and
# are not included in this repository; place them under the path below.
path <- "data/restricted/TaiwanTownship"
path_cih <- "data"

## Taiwan's data by township
PPmean_male <- fread(paste0(path, "/smoothed/PPmean_male.csv"))
PPmean_female <- fread(paste0(path, "/smoothed/PPmean_female.csv"))

DeathCount <- fread(paste0(path, "/DeathCount_v5.csv")) %>%
  filter(! Urbanicity_AST2023 == 9) %>%
  rename(COUNT_RAW = COUNT)

DeathCount_input <- bind_rows(
  DeathCount %>%
    filter(SEX == 1) %>%
    mutate(COUNT = PPmean_male$d),
  DeathCount %>%
    filter(SEX == 2) %>%
    mutate(COUNT = PPmean_female$d)
) %>%
  mutate(
    AGP = case_when(AGP == 1 ~ 0, TRUE ~ AGP)
  ) %>%
  summarise(
    POP = sum(POP),
    COUNT = sum(COUNT),
    .by = c(YEAR, SEX, MOI, AGP)
  )

colnames(DeathCount_input) <- c("year", "sex", "area", "age", "pop", "dths", "dths_raw")

rm(DeathCount, PPmean_female, PPmean_male)

######################
## e0 calculation   ##
######################

calc_e0 <- function(DeathCount_input, year_idx, sex_idx, area_idx){
  lt <- leCompute_dt_70q0(DeathCount_input, year_idx, sex_idx, area_idx)
  e0 <- lt[agp == "0-4", ex]
  return(e0)
}

#############################
## Monte Carlo Uncertainty ##
#############################

n_sim <- 1000
seed_mc <- 1011
n_cores <- getOption("mc.cores", max(1, detectCores() - 2, na.rm = TRUE))

calc_e0_mc <- function(DeathCount_input, year_idx, sex_idx, area_idx, cell_seed){
  lt_input <- leTable_dt(DeathCount_input, year_idx, sex_idx, area_idx)
  set.seed(cell_seed)
  sapply(seq_len(n_sim), function(s){
    sim_dt <- copy(lt_input)
    sim_dt[, dths := rpois(.N, dths)]
    leCalc_dt_70q0(sim_dt)[agp == "0-4", ex]
  })
}

e0_mcDraws <- function(DeathCount_input, combo, iso3_val, seed_base){
  DeathCount_input <- as.data.table(DeathCount_input)
  draws <- pbmclapply(seq_len(nrow(combo)), function(i){
    data.table(
      year = combo$year[i], sex = combo$sex[i], area = as.character(combo$area[i]),
      sim = seq_len(n_sim),
      e0 = calc_e0_mc(DeathCount_input, combo$year[i], combo$sex[i], combo$area[i], seed_base + i)
    )
  }, mc.cores = n_cores)
  rbindlist(draws)[, iso3 := iso3_val]
}

# Attach 95% UI (2.5th and 97.5th percentiles of the draws) to point estimates
addUI <- function(e0_dt, draws){
  ui <- draws[, .(e0_lower = quantile(e0, 0.025), e0_upper = quantile(e0, 0.975)),
              by = .(year, sex, area)]
  merge(e0_dt[, area := as.character(area)], ui, by = c("year", "sex", "area"), all.x = TRUE)
}

# --- 1. Taiwan Township ---
# combo = expand.grid(
#   year = c(2000, 2019),
#   sex = 1:2,
#   area = levels(as.factor(DeathCount_input$area))
# ) %>% as.data.table()
# 
# e0_byTownship2000_2019 <- combo[, e0 := calc_e0(
#   DeathCount_input,
#   year,
#   sex,
#   area
# ), by = .(year, sex, area)]
# 
# e0_byTownship2000_2019 <- e0_byTownship2000_2019 %>% mutate(iso3 = "TWN")
# write.csv(e0_byTownship2000_2019, "data/e0_taiwanTownship.csv", row.names = FALSE)

# --- 2. Taiwan County ---
DeathCount_county <- DeathCount_input %>%
  mutate(area = str_sub(area, 1, 5)) %>%
  summarise(
    pop = sum(pop),
    dths = sum(dths_raw),
    .by = c(year, sex, area, age)
  )

combo_county = expand.grid(
  year = c(2000, 2019),
  sex = 1:2,
  area = levels(as.factor(DeathCount_county$area))
) %>% as.data.table()

e0_byCounty2000_2019 <- combo_county[, e0 := calc_e0(
  DeathCount_county,
  year,
  sex,
  area
), by = .(year, sex, area)]

mcDraws_TWN <- e0_mcDraws(DeathCount_county, combo_county, "TWN", seed_mc)
e0_byCounty2000_2019 <- addUI(e0_byCounty2000_2019, mcDraws_TWN) %>% mutate(iso3 = "TWN")
write.csv(e0_byCounty2000_2019, file = paste0(path_cih, "/e0_taiwanCounty.csv"), row.names = FALSE)

rm(combo_county, DeathCount_county, DeathCount_input)

# --- 3. Japan Prefecture (2019) ---
DeathCount_japan2019 <- fread("data/JPNdata/DeathCount_japan2019.csv")

combo_JPN = expand.grid(
  year = 2019,
  sex = 1:2,
  area = levels(as.factor(DeathCount_japan2019$area))
) %>% as.data.table()

e0_japanPrefecture2019 <- combo_JPN[, e0 := calc_e0(
  DeathCount_japan2019,
  year,
  sex,
  area
), by = .(year, sex, area)]
mcDraws_JPN <- e0_mcDraws(DeathCount_japan2019, combo_JPN, "JPN", seed_mc + 1e4)
e0_japanPrefecture2019 <- addUI(e0_japanPrefecture2019, mcDraws_JPN) %>% mutate(iso3 = "JPN")

write.csv(e0_japanPrefecture2019, paste0(path_cih, "/e0_japanPrefecture2019.csv"), row.names = FALSE)
rm(combo_JPN, DeathCount_japan2019)


# --- 4. Korea Province (2019) ---
koreaDC2019 <- fread("data/KORdata/Deaths_and_Death_rates_by_cause_104_item___sex__and_age_fiveyear_age__province_20260314165539.csv")[-(1: 2), -1]
colnames(koreaDC2019) <- c("area", "agp", "sex", "dths")

levels(as.factor(koreaDC2019$agp))
ageMapping <- data.frame(
  agp = c(
    "0 year old", "1-4 Years old", "5-9 Years old", "10-14 Years old", "15-19 Years old", 
    "20-24 Years old", "25-29 Years old", "30-34 Years old", "35-39 Years old", "40-44 Years old",
    "45-49 Years old", "50-54 Years old", "55-59 Years old", "60-64 Years old", "65-69 Years old",
    "70-74 Years old", "75-79 Years old", "80-84 Years old", "85-89 Years old", "90 Years old and over"
  ),
  age = c(0, 0, seq(5, 85, 5), 85)
)


koreaDC2019 <- koreaDC2019 %>% 
  left_join(ageMapping, by = "agp") %>% 
  mutate(
    sex = ifelse(sex == "Male", 1, 2),
    dths = as.numeric(dths)
  ) %>% 
  summarise(
    dths = sum(dths),
    .by = c(area, sex, age)
  )

koreaPop2019 <- fread("data/KORdata/Resident_Population_in_FiveYear_Age_Groups__2011_._20260314170330.csv")[-c(1, 2), ]
colnames(koreaPop2019) <- c("area", "sex", paste0("a_", seq(0, 100, 5)))

ageMapping2 <- data.frame(
  agp = c(paste0("a_", seq(0, 100, 5))),
  age = c(seq(0, 85, 5), rep(85, 3))
)


koreaPop2019 <- melt(koreaPop2019, 
                     id.vars = c("area", "sex"), 
                     measure.vars = patterns("a_"),
                     variable.name = "agp", 
                     value.name = "pop") %>% 
  left_join(ageMapping2, by = "agp") %>% 
  mutate(
    area = case_when(area == "Jeju" ~ "Jeju-do", area == "sejong" ~ "Sejong", TRUE ~ area),
    sex = ifelse(sex == "Male (Person)", 1, 2),
    pop = as.numeric(pop)
  ) %>% 
  summarise(
    pop = sum(pop),
    .by = c(area, sex, age)
  )

DeathCount_korea2019 <- koreaPop2019 %>% 
  left_join(koreaDC2019, by = c("area", "sex", "age")) %>% 
  mutate(year = 2019)

rm(koreaDC2019, koreaPop2019, ageMapping, ageMapping2)

combo_KOR = expand.grid(
  year = 2019,
  sex = 1: 2,
  area = levels(as.factor(DeathCount_korea2019$area))
) %>% as.data.table()


e0_koreaProvince2019 <- combo_KOR[, e0 := calc_e0(
  DeathCount_korea2019,
  year,
  sex,
  area
), by = .(year, sex, area)]

mcDraws_KOR <- e0_mcDraws(DeathCount_korea2019, combo_KOR, "KOR", seed_mc + 2e4)
e0_koreaProvince2019 <- addUI(e0_koreaProvince2019, mcDraws_KOR) %>% mutate(iso3 = "KOR")
write.csv(e0_koreaProvince2019, paste0(path_cih, "/e0_koreaProvince2019.csv"), row.names = FALSE)
rm(combo_KOR, DeathCount_korea2019)

# --- 5. Distribution summaries with UI (2019) ---
mcDraws_all <- rbindlist(list(mcDraws_TWN, mcDraws_JPN, mcDraws_KOR))[year == 2019]
write.csv(mcDraws_all, paste0(path_cih, "/e0_subnational_mcDraws.csv"), row.names = FALSE)

distStats <- function(dt){
  dt[, .(min = min(e0), max = max(e0), range = max(e0) - min(e0), IQR = IQR(e0)),
     by = setdiff(names(dt), c("area", "e0"))]
}

e0_subnational2019 <- rbindlist(list(
  e0_byCounty2000_2019[year == 2019], e0_japanPrefecture2019, e0_koreaProvince2019
), use.names = TRUE)

stats_point <- distStats(e0_subnational2019[, .(iso3, sex, area, e0)]) %>%
  melt(id.vars = c("iso3", "sex"), variable.name = "stat", value.name = "value")
stats_ui <- distStats(mcDraws_all[, .(iso3, sex, sim, area, e0)]) %>%
  melt(id.vars = c("iso3", "sex", "sim"), variable.name = "stat") %>%
  .[, .(lower = quantile(value, 0.025), upper = quantile(value, 0.975)), by = .(iso3, sex, stat)]

e0_subnationalSummary <- merge(stats_point, stats_ui, by = c("iso3", "sex", "stat"))
write.csv(e0_subnationalSummary, paste0(path_cih, "/e0_subnational_summary.csv"), row.names = FALSE)



