options(stringsAsFactors = FALSE)


### Load packages ==================================================================================
library(here)

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
})


### Source functions ===============================================================================



### Environment ====================================================================================
output_directory <- here("outputs", "01-sample_sheet")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0755")

project_name <- gsub("(.*)_.*", "\1", basename(here()))

data_directory <- "/disks/DATA/Projects/EpxGDM/FINNGEDI_TUBINGEN/EPIC"



### Define theme ===================================================================================
# options(ggplot2.continuous.colour = "viridis", ggplot2.continuous.fill = "viridis")
# theme_set(theme_light())


### Check Sample Sheet =============================================================================
finngedi_phenotypes <- inner_join(
  x = read_excel(here("docs", "FinnGedi_GDM_methylation_plate_design_Apr2017.xls"), sheet = 1),
  y = read_excel(here("docs", "FinnGedi_GDM_methylation_plate_design_Apr2017.xls"), sheet = 2),
  by = c("ID" = "id")
) %>% 
  select(-matches("rownum")) %>% 
  select_if(~ !all(is.na(.x))) %>%
  select(
    Sample_ID = id_str, ID, GDM, 
    AGE_m = mumage, BMI_m = prebmi, SEX = sex, gestationalweek = gestwk, status = status1,
    GDMwk, ogtt, insulin = Insulin, metformin = Metformin, 
    birthweight_c = bwt,
    Plate, ROW, COL#, rep, comment
  ) %>% 
  mutate(
    cohort = "FinnGeDi",
    SEX = ifelse(status == "Mother", 2, SEX),
    status = gsub("Child", "Offspring", status)
  ) %>% 
  arrange(status) %>%
  group_by(ID) %>% 
  fill() %>% 
  ungroup()

tubingen_phenotypes <- read_excel(here("docs", "2018-06-13_GDM_EPI_Tubingen_phenotype_data.xlsx"), sheet = 1) %>% 
  select(
    Sample_ID = NR, 
    AGE_m = AGE,
    BMI_m = BMIprepreg,
    GDM = `group (NGT0/ GDM1)`,
    gestationalweek = `SSW (Schwangerschaftswochen = Weeks of pregnancy)`
  ) %>% 
  mutate(
    cohort = "Tubingen", 
    Sample_ID = as.character(Sample_ID), 
    SEX = 2,
    status = "Mother"
  )

sample_sheet <- read_csv(file.path(data_directory, "SampleSheet_20180116.csv"), skip = 8) %>% 
  select_if(~ !all(is.na(.x))) %>% 
  rename(
    Replicate = Replicate_1,
    Sentrix_ID = SentrixBarcode_A,
    Sentrix_Position = SentrixPosition_A,
    Sample_PlateDigit = Sample_Plate,
    Sample_Plate = AMP_Plate
  ) %>% 
  mutate(
    Sample_PlateDigit = gsub(".*_plate", "", Project),
    Project = case_when(
      grepl("Finn", Project) ~ "FinnGeDi",
      grepl("Tubingen", Project) ~ "Tubingen",
      TRUE ~ "Control"
    ),
    Sample_ID = gsub("bis", "_r", Sample_ID),
    Sample_Name = Sample_ID,
    Sample_IID = gsub("_r", "", Sample_ID),
    Replicate = gsub("bis", "_r", Replicate)
  )


# Check for anomalies
left_join(
  x = select(finngedi_phenotypes, Sample_ID, Plate, ROW, COL),
  y = select(sample_sheet, Sample_ID, Sample_PlateDigit, Sample_Well),
  by = "Sample_ID"
) %>%
  mutate(COL = ifelse(Plate == 12, COL + 6, COL)) %>%
  filter(Sample_Well != paste0(ROW, sprintf("%02d", COL)))
  

write_csv(
  x = inner_join(
    x = mutate(sample_sheet, Sample_ID = gsub("(3500)_r|(3939)_r", "\\1\\2", Sample_ID)),
    y = select(bind_rows(finngedi_phenotypes, tubingen_phenotypes), -Plate, -ROW, -COL),
    by = "Sample_ID"
  ), 
  path = file.path(output_directory, "samplesheet.csv")
)
