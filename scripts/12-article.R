# Docker Container R403
### Environment ====================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", basename(here()))
output_directory <- here("outputs", "12-article")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

data_directory <- normalizePath("/disks/DATA/Projects/EpxGDM/FINNGEDI_TUBINGEN/QC")


### Load packages ==================================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(ggbeeswarm)
  library(data.table)
  library(readxl)
  library(flextable)
  library(gt)
  library(gtsummary)
  library(forcats)
  library(ragg)
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(glue)
  library(purrr)
  library(modelr)
  library(broom)
})


### Define theme ===================================================================================
options(
  ggplot2.discrete.colour = function(...) scale_colour_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.colour = function(...) scale_colour_viridis_c(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.fill = function(...) scale_fill_viridis_c(..., begin = 0.15, end = 0.85)
)
theme_set(
  theme_light(base_family = "Times", base_size = 12) + 
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.8)),
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown(),
      legend.title = element_markdown(),
      # legend.text = element_markdown(),
      strip.background = element_rect(fill = "white", colour = "grey70"),
      strip.text = element_text(colour = "black"),
      panel.grid = element_blank()
    )
)

decimal_mark <- "·"
set_gtsummary_theme(list(
  "style_number-arg:big.mark" = "",
  "style_number-arg:decimal.mark" = decimal_mark
))


### Source functions ===============================================================================
pval_trans <- function(alpha = NULL, md = FALSE, prefix = FALSE, colour = "firebrick") {
  scales::trans_new(
    name = "pval",
    domain = c(0, 1),
    transform = function(x) {x[x < .Machine$double.xmin] <- .Machine$double.xmin; -log(x, 10)},
    inverse = function(x) {10^-x},
    breaks = (function(n = 5) {
      function(x) {
        max <- floor(-log(min(c(x, alpha), na.rm = TRUE), base = 10))
        if (max == 0) 1 else sort(unique(c(10^-seq(0, max, by = floor(max / n) + 1), alpha)))
      }
    })(),
    format = (function(x) {
      if (md & nchar(system.file(package = "ggtext")) != 0) {
        prefix_text <- if (prefix) "&alpha; = " else ""
        x_fmt <- gsub(
          "^(.*)e[+]*([-]*)0*(.*)$", 
          "\\1 &times; 10<sup>\\2\\3</sup>", 
          format(x, scientific = TRUE)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- gsub("^1 &times; ", "", x_fmt)
        alpha_idx <- format(x, scientific = TRUE) == format(alpha, scientific = TRUE)
        x_fmt[alpha_idx] <- paste0("<b style='color:", colour, ";'>", prefix_text, x_fmt[alpha_idx], "</b>")
        x_fmt
      } else {
        prefix_text <- if (prefix) "alpha == " else ""
        x_fmt <- gsub(
          "^(.*)e[+]*([-]*)0*(.*)$", 
          "\\1 %*% 10^\\2\\3", 
          format(x, scientific = TRUE)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- gsub("^1 \\%\\*\\% ", "", x_fmt)
        alpha_idx <- format(x, scientific = TRUE) == format(alpha, scientific = TRUE)
        x_fmt[alpha_idx] <- paste0(prefix_text, x_fmt[alpha_idx])
        parse(text = x_fmt)
      }
    })
  )
}

fortify.manhattan <- function(data, x, y, group) {
  map_chro <- c(seq(22), "X", "Y", "X", "Y")
  names(map_chro) <- c(seq(24), "X", "Y")
  
  `:=` <- data.table::`:=`
  
  out <- data.table::as.data.table(data)
  data.table::setnames(out, c(x, y, group), c("x_pos", "y_pval", "x_chr"))
  out[, x_chr := as.character(x_chr)]
  out[, x_chr := map_chro[gsub("^chr", "", x_chr, ignore.case = TRUE)]]
  out[, x_chr := factor(x_chr, levels = intersect(c(seq(22), "X", "Y"), x_chr))]
  out[, x_pos := as.double(x_pos)]
  out[order(x_chr, x_pos)]
  out[, x_pos := scales::rescale(x = x_pos, to = c(-0.4, 0.4)), by = "x_chr"]
  out[, x_pos := x_pos + as.integer(x_chr)]
  data.table::setnames(out, c("x_pos", "y_pval", "x_chr"), c("x", "y", "group"))
  out[]
}

StatManhattan <- ggplot2::ggproto("StatManhattan", ggplot2::Stat,
  required_aes = c("x", "y", "group"),
  setup_data = function(data, params) {
    fortify.manhattan(data, "x", "y", "group")
  },
  compute_layer = function(data, scales, params) {
    data
  }
)

logit2 <- function(x) log2(x) - log2(1 - x)
ilogit2 <- function(x) 2^x / (1 + 2^x)

logit2_trans <- function(suffix = "", breaks = FALSE, drop0trailing = FALSE) {
  logit2 <- function(x) log2(x) - log2(1 - x)
  ilogit2 <- function(x) 2^x / (1 + 2^x)
  logit2_breaks <- function(n = 5) {
    function(x) {
      rng <- range(logit2(x))
      signif(ilogit2(scales::extended_breaks(n = n)(rng)), 2)
    }
  }
  
  if (breaks) {
    scales::trans_new(
      name = "logit2",
      transform = logit2,
      inverse = ilogit2,
      breaks = logit2_breaks(n = 5),
      domain = c(0, Inf), 
      format = scales::percent_format(suffix = suffix, drop0trailing = drop0trailing)
    )
  } else {
    scales::trans_new(
      name = "logit2",
      transform = logit2,
      inverse = ilogit2,
      domain = c(0, Inf), 
      format = scales::percent_format(suffix = suffix, drop0trailing = drop0trailing)
    )
  }
}


### EWAS Offspring adjusted for Mothers ============================================================
res_om <- fread(here(
  "outputs", "09-ewas_finngedi_h", "GDM", 
  "EpxGDM_EWAS_DMP_GDM_mothers_offsprings_simple_cell.csv.gz"
))[, cpg_chr := factor(cpg_chr, levels = paste0("chr", 1:22))]
res_om[UCSC_RefGene_Name == "HIST1H3E", UCSC_RefGene_Name := "H3C6"]
res_om[fdr < 0.05, label := fifelse(UCSC_RefGene_Name == "", CpG, paste0(CpG, "<br>(<i>", UCSC_RefGene_Name, "</i>)"))]


### Phenotypes Data ================================================================================
phenotype_matrix <- fread(file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes.csv"))
excel_pheno <- as.data.table(read_xlsx(
  here("docs/FinnGedi_clinical_data_2018_10_05/FinnGeDi_selection_EpxGDM_2018_10_06_TO_TOBY.xlsx")
))[
  ID %in% phenotype_matrix[, .N, by = "ID"][N == 2, ID]
][, 
  sex := factor(c("1" = "Male", "2" = "Female")[as.character(sex)], levels = c("Male", "Female"))
]

### Table 1 ========================================================================================
tbl_summary(
  data = excel_pheno, 
  by = "GDM", 
  include = c("GDM", "mumage", "prebmi", "gestwk", "sex", "bwt"), 
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} ({p})"
  ),
  label = list(
    GDM ~ "GDM",
    mumage ~ "Age of the mother (years)",
    prebmi ~ "BMI of the mother (kg/m²)",
    gestwk ~ "Gestation (weeks)",
    sex ~ "Sex of the offspring",
    bwt ~ "Birthweight (g)"
  )
) %>% 
  add_p(
    test = list(
      all_categorical() ~ "fisher.test", 
      all_continuous() ~ "t.test"
    )
  ) %>% 
  modify_header(
    stat_by = "**{c('0' = 'Control', '1' = 'GDM')[level]}** (N = {n})"
  ) %>% 
  bold_labels() %>% 
  italicize_levels() %>% 
  as_flex_table() %>% 
  save_as_docx(path = file.path(output_directory, "table1.docx"))

tbl_summary(
  data = excel_pheno, 
  by = "GDM", 
  # include = c("GDM", "mumage", "prebmi", "gestwk", "sex", "bwt"), 
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{n} ({p})"
  ),
  label = list(
    GDM ~ "GDM",
    mumage ~ "Age of the mother (years)",
    prebmi ~ "BMI of the mother (kg/m²)",
    gestwk ~ "Gestation (weeks)",
    sex ~ "Sex of the offspring",
    bwt ~ "Birthweight (g)"
  )
) %>% 
  add_p(
    test = list(
      all_categorical() ~ "fisher.test", 
      all_continuous() ~ "t.test"
    )
  ) %>% 
  modify_header(
    stat_by = "**{c('0' = 'Control', '1' = 'GDM')[level]}** (N = {n})"
  ) %>% 
  bold_labels() %>% 
  italicize_levels() %>% 
  as_flex_table() %>% 
  save_as_docx(path = file.path(output_directory, "table1_big.docx"))

### Supplementary Table 1 ==========================================================================
tbl_cross(
  data = excel_pheno[,
    .(
      "row" = factor(mumagecategorical, labels = c("< 25", "25-30", "30-35", "> 35")),
      "col" = factor(
        paste(
          factor(BMIcategorical, labels = c("≤ 24.9", "25.0-29.9", "≥ 30.0")), 
          factor(GDM, labels = c("Control", "GDM"))
        ),
        levels = unlist(lapply(c("≤ 24.9", "25.0-29.9", "≥ 30.0"), paste, c("Control", "GDM")))
      )
    )
  ], 
  row = row, 
  col = col, 
  percent = "cell", 
  label = list(row ~ "Age of the mother (years)", col ~ "BMI of the mother (kg/m²)"), 
  statistic = "{n} ({p} %)", 
  margin = NULL
) %>% 
  bold_labels() %>% 
  italicize_levels() %>% 
  as_flex_table() %>% 
  save_as_docx(path = file.path(output_directory, "supplementary_table1.docx"))


### EPIC ===========================================================================================
beta_matrix <- vroom(file.path(data_directory, "EPIC", "EPIC_QC_betavalues.csv.gz"))
samples_cols <- setdiff(colnames(beta_matrix), "cpg_id")
beta_matrix <- setDT(beta_matrix)[cpg_id %in% res_om[fdr < 0.05, CpG]]
beta_matrix[, (samples_cols) := lapply(.SD, function(x) log2(x) - log2(1 - x)), .SDcols = samples_cols]
beta_matrix <- `rownames<-`(as.matrix(beta_matrix[, ..samples_cols]), beta_matrix[["cpg_id"]])

phenotype_pcs <- phenotype_matrix %>% 
  as_tibble() %>% 
  filter(!qc_gender_discrepancy & call_rate >= 0.99) %>% 
  filter(Project == "FinnGeDi") %>% 
  filter(Sample_ID %in% colnames(beta_matrix)) %>% 
  group_by(Sample_IID) %>% 
  slice(which.max(call_rate)) %>% 
  ungroup() %>% 
  group_by(gsub("[[:alpha:]]", "", Sample_ID)) %>% 
  filter(n() == 2) %>% 
  ungroup()

default_covariates <- c(
  mothers = "AGE_m + BMI_m",
  mothers_offsprings = "SEX + log(birthweight_c) + gestationalweek"
)
default_cells <- phenotype_pcs %>% 
  select(status, starts_with("CellT")) %>% 
  group_by(status) %>% 
  summarise_all(~all(!is.na(.x))) %>% 
  ungroup() %>% 
  pivot_longer(cols = starts_with("CellT"), names_to = "CellT", values_to = "available") %>% 
  filter(available) %>% 
  group_by(status) %>% 
  summarise(CellT = paste(CellT, collapse = " + ")) %>% 
  (function(x) {
    setNames(
      object = x[["CellT"]], 
      nm = tolower(ifelse(grepl("Offspring", x[["status"]]), "mothers_offsprings", paste0(x[["status"]], "s")))
    )
  })()

om_status <- "mothers_offsprings"
trait <- "GDM"
models <- tribble(
  ~model, ~formula, ~formula_m,
  "simple_cell", as.formula(glue("~ { trait } + { default_covariates[om_status] } + { default_cells[om_status] }")), 
    as.formula(glue('~ { default_covariates["mothers"] } + { default_cells["mothers"] }'))
)

data_models <- models %>%
  mutate(
    pheno = pmap(
      .l = list(formula, model, formula_m),
      .f = function(.form, .model, .form_m) {
        phenotype_pcs %>%
          select(Sample_ID, all.vars(.form), all.vars(.form_m)) %>%
          group_by(group = gsub("[[:alpha:]]", "", Sample_ID)) %>% 
          filter(n() == 2) %>% 
          ungroup() %>% 
          select(-group)
      }
    ),
    beta_offsprings = map(pheno, ~ beta_matrix[, filter(.x, grepl("c", Sample_ID))[["Sample_ID"]]]),
    beta_mothers = map(pheno, ~ beta_matrix[, filter(.x, grepl("m", Sample_ID))[["Sample_ID"]]])
  )

data_fdr <- rbindlist(lapply(
    X = intersect(
    rownames(data_models[[1, "beta_offsprings"]][[1]]),
    rownames(data_models[[1, "beta_mothers"]][[1]])
  ),
  .form = data_models[[1, "formula"]][[1]],
  .form_m = data_models[[1, "formula_m"]][[1]],
  .beta_o = data_models[[1, "beta_offsprings"]][[1]],
  .beta_m = data_models[[1, "beta_mothers"]][[1]],
  .pheno = data_models[[1, "pheno"]][[1]],
  FUN = function(icpg, .form, .form_m, .beta_o, .beta_m, .pheno) {
    beta_icpg <- c(.beta_o[icpg, ], .beta_m[icpg, ])
    outliers <- c(
      quantile(beta_icpg, probs = 0.25) - 3 * IQR(beta_icpg),
      quantile(beta_icpg, probs = 0.75) + 3 * IQR(beta_icpg)
    )
    lm_mother_resid <- .pheno %>%
      filter(grepl("m", Sample_ID)) %>%
      mutate(Mothers = as.matrix(.beta_m)[icpg, Sample_ID]) %>%
      do(
        add_residuals(
          data = .,
          model = lm(formula = update.formula(.form_m, "Mothers ~ ."), data = .)
        )
      ) %>% 
      select(Sample_ID, Mothers = resid)
    
    lm_pheno <- .pheno %>%
      filter(grepl("c", Sample_ID)) %>%
      mutate(
        Offsprings = .beta_o[icpg, Sample_ID],
        Mothers_raw = .beta_m[icpg, gsub("c", "m", Sample_ID)]
      ) %>% 
      filter(
        between(Offsprings, outliers[1], outliers[2]) &
          between(Mothers_raw, outliers[1], outliers[2])
      ) %>% 
      mutate(MID = gsub("c", "m", Sample_ID)) %>% 
      left_join(y = lm_mother_resid, by = c("MID" = "Sample_ID")) %>% 
      select(Sample_ID, all.vars(.form), Offsprings, Mothers, Mothers_raw) %>%
      drop_na()

    lm_pheno$cpg_id <- icpg
    lm_pheno$Offspring <- ilogit2(lm_pheno$Offsprings)
    lm_pheno$Offsprings <- NULL
    lm_pheno$Mother <- ilogit2(lm_pheno$Mothers_raw)
    lm_pheno$Mothers_raw <- NULL
    lm_pheno$Mother_residual_logit2 <- lm_pheno$Mothers
    lm_pheno$Mothers <- NULL
    lm_pheno
  }
))
fwrite(x = data_fdr, file = file.path(output_directory, "data_fdr.csv.gz"))


### Table 2 ========================================================================================
dmp <- rbindlist(lapply(
    X = intersect(
    rownames(data_models[[1, "beta_offsprings"]][[1]]),
    rownames(data_models[[1, "beta_mothers"]][[1]])
  ),
  .form = data_models[[1, "formula"]][[1]],
  .form_m = data_models[[1, "formula_m"]][[1]],
  .beta_o = data_models[[1, "beta_offsprings"]][[1]],
  .beta_m = data_models[[1, "beta_mothers"]][[1]],
  .pheno = data_models[[1, "pheno"]][[1]],
  FUN = function(icpg, .form, .form_m, .beta_o, .beta_m, .pheno) {
    beta_icpg <- c(.beta_o[icpg, ], .beta_m[icpg, ])
    outliers <- c(
      quantile(beta_icpg, probs = 0.25) - 3 * IQR(beta_icpg),
      quantile(beta_icpg, probs = 0.75) + 3 * IQR(beta_icpg)
    )
    lm_mother_resid <- .pheno %>%
      filter(grepl("m", Sample_ID)) %>%
      mutate(Mothers = as.matrix(.beta_m)[icpg, Sample_ID]) %>%
      do(
        add_residuals(
          data = .,
          model = lm(formula = update.formula(.form_m, "Mothers ~ ."), data = .)
        )
      ) %>% 
      select(Sample_ID, Mothers = resid)
    
    lm_pheno <- .pheno %>%
      filter(grepl("c", Sample_ID)) %>%
      mutate(
        Offsprings = .beta_o[icpg, Sample_ID],
        Mothers_raw = .beta_m[icpg, gsub("c", "m", Sample_ID)]
      ) %>% 
      filter(
        between(Offsprings, outliers[1], outliers[2]) &
          between(Mothers_raw, outliers[1], outliers[2])
      ) %>% 
      mutate(MID = gsub("c", "m", Sample_ID)) %>% 
      left_join(y = lm_mother_resid, by = c("MID" = "Sample_ID")) %>% 
      select(Sample_ID, all.vars(.form), Offsprings, Mothers) %>%
      drop_na()

    tidy(lm(
      formula = update.formula(.form, Offsprings ~ Mothers + Mothers:GDM + .),
      data = lm_pheno
    )) %>%
      filter(grepl("GDM", term)) %>%
      mutate(
        CpG = icpg,
        avgmvalue_meth = mean(lm_pheno[["Offsprings"]]),
        avgbvalue_meth = mean(ilogit2(lm_pheno[["Offsprings"]])),
        n = nrow(lm_pheno)
      )
  }
))
merge(
  x = res_om[fdr < 0.05], 
  y = dmp[, c("CpG", "term", "n", "avgbvalue_meth")], 
  by = c("CpG", "term", "n")
)[order(pvalue),
  .(
    term, CpG, 
    estimate = paste(
      format(estimate, trim = TRUE, digits = 2, nsmall = 2, decimal.mark = decimal_mark), 
      format(se, trim = TRUE, digits = 2, nsmall = 2, decimal.mark = decimal_mark),
      sep = " ± "
    ), 
    avgbvalue_meth,
    n, pvalue, fdr, 
    Relation_to_Island = c(
      "OpenSea" = "Open Sea", "Island" = "Island", "N_Shore" = "North Shore", "S_Shelf" = "South Shelf"
    )[Relation_to_Island], 
    UCSC_RefGene_Name = fifelse(UCSC_RefGene_Name != "", paste0("*", UCSC_RefGene_Name, "*"), NA_character_), 
    UCSC_RefGene_Group = fifelse(UCSC_RefGene_Group != "", UCSC_RefGene_Group, NA_character_)
  )
][, term := gsub("Mothers:GDM", "GDM &times; Mothers", term)] %>% 
  gt() %>%
  fmt_scientific("pvalue", drop_trailing_zeros = FALSE, decimals = 2, dec_mark = decimal_mark) %>% 
  fmt_scientific("fdr", drop_trailing_zeros = FALSE, decimals = 2, dec_mark = decimal_mark) %>% 
  fmt_percent(
    columns = "avgbvalue_meth", 
    decimals = 2, drop_trailing_zeros = FALSE, incl_space = TRUE, dec_mark = decimal_mark
  ) %>% 
  fmt_markdown("UCSC_RefGene_Name") %>% 
  fmt_missing(c("UCSC_RefGene_Name", "UCSC_RefGene_Group")) %>% 
  fmt_markdown("term") %>% 
  tab_footnote(
    footnote = paste(
      'Variable term or interaction term from the linear regression,',
      'where "Mothers" is the methylation of the mothers as Mvalue.'
    ), 
    locations = cells_column_labels("term")
  ) %>% 
  tab_footnote(
    footnote = "Linear regression coefficient plus/minus standard error on the Mvalue scale.", 
    locations = cells_column_labels("estimate")
  ) %>% 
  tab_footnote(
    footnote = md(paste(
      "Number of Mother/Offspring pair after outliers exclusion,",
      "with Mvalue greater than 3 times the interquartile",
      "range below the 25<sup>th</sup> percentile and lower than 3 times",
      "the interquartile range above the 75<sup>th</sup> percentile."
    )), 
    locations = cells_column_labels("n")
  ) %>% 
  cols_label(
    term = "Term",
    CpG = "CpG ID",
    estimate = "Estimate",
    avgbvalue_meth = md("Average β<sub>value</sub>"),
    n = "Number of Mother/Offspring Pair",
    fdr = "FDR",
    Relation_to_Island = "Relation to CpG Island",
    UCSC_RefGene_Name = "Gene",
    UCSC_RefGene_Group = "Location in Gene"
  ) %>% 
  gtsave(filename = file.path(output_directory, "table2.rtf"))


### Table 2 Full ===================================================================================
dmp <- rbindlist(lapply(
    X = intersect(
    rownames(data_models[[1, "beta_offsprings"]][[1]]),
    rownames(data_models[[1, "beta_mothers"]][[1]])
  ),
  .form = data_models[[1, "formula"]][[1]],
  .form_m = data_models[[1, "formula_m"]][[1]],
  .beta_o = data_models[[1, "beta_offsprings"]][[1]],
  .beta_m = data_models[[1, "beta_mothers"]][[1]],
  .pheno = data_models[[1, "pheno"]][[1]],
  FUN = function(icpg, .form, .form_m, .beta_o, .beta_m, .pheno) {
    beta_icpg <- c(.beta_o[icpg, ], .beta_m[icpg, ])
    outliers <- c(
      quantile(beta_icpg, probs = 0.25) - 3 * IQR(beta_icpg),
      quantile(beta_icpg, probs = 0.75) + 3 * IQR(beta_icpg)
    )
    lm_mother_resid <- .pheno %>%
      filter(grepl("m", Sample_ID)) %>%
      mutate(Mothers = as.matrix(.beta_m)[icpg, Sample_ID]) %>%
      do(
        add_residuals(
          data = .,
          model = lm(formula = update.formula(.form_m, "Mothers ~ ."), data = .)
        )
      ) %>% 
      select(Sample_ID, Mothers = resid)
    
    lm_pheno <- .pheno %>%
      filter(grepl("c", Sample_ID)) %>%
      mutate(
        Offsprings = .beta_o[icpg, Sample_ID],
        Mothers_raw = .beta_m[icpg, gsub("c", "m", Sample_ID)]
      ) %>% 
      filter(
        between(Offsprings, outliers[1], outliers[2]) &
          between(Mothers_raw, outliers[1], outliers[2])
      ) %>% 
      mutate(MID = gsub("c", "m", Sample_ID)) %>% 
      left_join(y = lm_mother_resid, by = c("MID" = "Sample_ID")) %>% 
      select(Sample_ID, all.vars(.form), Offsprings, Mothers) %>%
      drop_na()

    tidy(lm(
      formula = update.formula(.form, Offsprings ~ Mothers + Mothers:GDM + .),
      data = lm_pheno
    )) %>%
      filter(grepl("GDM|mother", term, ignore.case = TRUE)) %>%
      mutate(
        CpG = icpg,
        avgmvalue_meth = mean(lm_pheno[["Offsprings"]]),
        avgbvalue_meth = mean(ilogit2(lm_pheno[["Offsprings"]])),
        n = nrow(lm_pheno)
      )
  }
))

merge(
  x = unique(
    res_om[
      fdr < 0.05, 
      c("CpG", "n", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")
    ]
  ), 
  y = dmp, 
  by = c("CpG", "n")
)[
  j = .(
    term, CpG = factor(CpG, levels = res_om[fdr < 0.05][order(pvalue), unique(CpG)]), 
    # estimate = paste(
    #   format(estimate, trim = TRUE, digits = 2, nsmall = 2, decimal.mark = decimal_mark), 
    #   format(std.error, trim = TRUE, digits = 2, nsmall = 2, decimal.mark = decimal_mark),
    #   sep = " ± "
    # ), 
    estimate,
    std.error,
    avgbvalue_meth,
    n, p.value, 
    Relation_to_Island = c(
      "OpenSea" = "Open Sea", "Island" = "Island", "N_Shore" = "North Shore", "S_Shelf" = "South Shelf"
    )[Relation_to_Island], 
    UCSC_RefGene_Name = fifelse(UCSC_RefGene_Name != "", paste0("*", UCSC_RefGene_Name, "*"), NA_character_), 
    UCSC_RefGene_Group = fifelse(UCSC_RefGene_Group != "", UCSC_RefGene_Group, NA_character_)
  )
][
  j = term := gsub("Mother:GDM", "GDM &times; Mother", gsub("Mothers", "Mother", term))
][
  j = term := factor(term, levels = c("GDM", "GDM &times; Mother", "Mother"))
][order(CpG, term)][j = term := as.character(term)] %>% 
  gt() %>%
  fmt_scientific("p.value", drop_trailing_zeros = FALSE, decimals = 2, dec_mark = decimal_mark) %>% 
  fmt_percent(
    columns = "avgbvalue_meth", 
    decimals = 2, drop_trailing_zeros = FALSE, incl_space = TRUE, dec_mark = decimal_mark
  ) %>% 
  fmt_markdown("UCSC_RefGene_Name") %>% 
  fmt_missing(c("UCSC_RefGene_Name", "UCSC_RefGene_Group")) %>% 
  fmt_markdown("term") %>% 
  fmt_number(columns = vars(estimate, std.error), n_sigfig = 3, drop_trailing_zeros = FALSE) %>% 
  cols_merge_uncert(col_val = "estimate", col_uncert = "std.error") %>% 
  tab_footnote(
    footnote = paste(
      'Variable term or interaction term from the linear regression,',
      'where "Mothers" is the methylation of the mothers as Mvalue.'
    ), 
    locations = cells_column_labels("term")
  ) %>% 
  tab_footnote(
    footnote = "Linear regression coefficient plus/minus standard error on the Mvalue scale.", 
    locations = cells_column_labels("estimate")
  ) %>% 
  tab_footnote(
    footnote = md(paste(
      "Number of Mother/Offspring pair after outliers exclusion,",
      "with Mvalue greater than 3 times the interquartile",
      "range below the 25<sup>th</sup> percentile and lower than 3 times",
      "the interquartile range above the 75<sup>th</sup> percentile."
    )), 
    locations = cells_column_labels("n")
  ) %>% 
  cols_label(
    term = "Term",
    CpG = "CpG ID",
    estimate = "Estimate",
    p.value = "P-Value",
    avgbvalue_meth = md("Average β<sub>value</sub>"),
    n = "Number of Mother/Offspring Pair",
    Relation_to_Island = "Relation to CpG Island",
    UCSC_RefGene_Name = "Gene",
    UCSC_RefGene_Group = "Location in Gene"
  ) %>% 
  cols_move(columns = "avgbvalue_meth", after = "p.value") %>% 
  cols_move(columns = "n", after = "avgbvalue_meth") %>% 
  tab_options(
    column_labels.font.weight = "bold", 
    data_row.padding = 0.05, 
    container.width = px(800), 
    table.font.size = "7pt"
  ) %>% 
  gtsave(filename = file.path(output_directory, "table2_full.html"))


### Supplementary Table 3 ==========================================================================
fwrite(
  x = fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_offsprings_simple_cell.csv.gz"
  ))[pvalue < 0.05],
  file = file.path(output_directory, "supplementary_table3.csv")
)


### Supplementary Table 4 ==========================================================================
fwrite(
  x = fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_mothers_simple_cell.csv.gz"
  ))[pvalue < 0.05],
  file = file.path(output_directory, "supplementary_table4.csv")
)


### Figure 1 =======================================================================================
cpg_order <- res_om[fdr < 0.05][order(pvalue), unique(label)]

p_a <- ggplot(
  data = res_om[order(pvalue)][,
    c("exppval", "labels") := list(
      (1:.N - 0.5) / .N,
      paste0(
        "&lambda;<sub>", gsub("Mothers:GDM", "GDM &times; Mothers", term), "</sub> = ", 
        format(median(qnorm(pvalue / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1), digits = 3, nsmall = 3)
      )
    ),
    by = "term"
  ][, labels := factor(labels, levels = sort(unique(labels), decreasing = TRUE))]
) +
  aes(x = exppval, y = pvalue, colour = labels, shape = labels) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = pval_trans(), expand = expansion(add = c(0, 0.2)), limits = c(1, NA)) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_manual(
    values = c(
      viridis_pal(begin = 0.5, end = 0.5, option = "viridis")(1),
      viridis_pal(begin = 0.5, end = 0.5, option = "plasma")(1)
    ), 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      ncol = 2, 
      override.aes = list(size = 2)
    )
  ) +
  scale_shape_discrete(
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      ncol = 2
    )
  ) +
  labs(x = "Expected P-value", y = "Observed P-value") +
  coord_cartesian(xlim = c(1, res_om[, min(pvalue)]), ylim = c(1, res_om[, min(pvalue)])) +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

p_b <- ggplot(data = res_om[term == "GDM"]) +
  aes(x = estimate, y = pvalue, colour = abs(estimate)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(data = ~ .x[fdr >= 0.05], size = 0.5) + 
  geom_point(
    mapping = aes(shape = factor(label, levels = cpg_order)),
    data = ~ .x[fdr < 0.05],
    size = 1.25,
    colour = "firebrick"
  ) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_viridis_c(trans = "sqrt", limits = c(0, NA), guide = "none", option = "viridis") +
  scale_shape_manual(
    values = 1:7,
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      byrow = TRUE, 
      nrow = 3, 
      override.aes = list(size = 2)
    )
  ) +
  labs(x = "Estimate (M<sub>value</sub>)", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

p_c <- ggplot(data = res_om[term == "Mothers:GDM"]) +
  aes(x = estimate, y = pvalue, colour = abs(estimate)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(data = ~ .x[fdr >= 0.05], size = 0.5) + 
  geom_point(
    mapping = aes(shape = factor(label, levels = cpg_order)),
    data = ~ .x[fdr < 0.05],
    size = 1.25,
    colour = "firebrick"
  ) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_viridis_c(trans = "sqrt", limits = c(0, NA), guide = "none", option = "plasma") +
  scale_shape_manual(
    values = 1:7,
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      byrow = TRUE, 
      nrow = 3, 
      override.aes = list(size = 2)
    )
  ) +
  labs(x = "Estimate (M<sub>value</sub>)", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

p_d <- ggplot(data = fortify.manhattan(res_om[term == "GDM"], "cpg_pos", "pvalue", "cpg_chr")) +
  aes(x = x, y = y, colour = group) +
  geom_point(data = ~ .x[fdr >= 0.05], size = 0.5) +
  geom_point(
    mapping = aes(shape = factor(label, levels = cpg_order)),
    data = ~ .x[fdr < 0.05],
    size = 1.25,
    colour = "firebrick"
  ) +
  scale_x_continuous(breaks = 1:24, labels = c(1:22, "X", "Y"), expand = expansion(add = 0.25)) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_manual(values = rep(viridis_pal(begin = 1/4, end = 3/4, option = "viridis")(2), 12), guide = "none") +
  scale_shape_manual(
    values = 1:7,
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      byrow = TRUE, 
      nrow = 3, 
      override.aes = list(size = 2)
    )
  ) +
  labs(x = "Chromosome", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

p_e <- ggplot(data = fortify.manhattan(res_om[term == "Mothers:GDM"], "cpg_pos", "pvalue", "cpg_chr")) +
  aes(x = x, y = y, colour = group) +
  geom_point(data = ~ .x[fdr >= 0.05], size = 0.5) +
  geom_point(
    mapping = aes(shape = factor(label, levels = cpg_order)),
    data = ~ .x[fdr < 0.05],
    size = 1.25,
    colour = "firebrick"
  ) +
  scale_x_continuous(breaks = 1:24, labels = c(1:22, "X", "Y"), expand = expansion(add = 0.25)) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_manual(values = rep(viridis_pal(begin = 1/4, end = 3/4, option = "plasma")(2), 12), guide = "none") +
  scale_shape_manual(
    values = 1:7,
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      byrow = TRUE, 
      nrow = 3, 
      override.aes = list(size = 2)
    )
  ) +
  labs(x = "Chromosome", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

agg_tiff(
  filename = file.path(output_directory, "figure1.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
  print(
    wrap_plots(
      p_b, p_d,
      p_c, p_e,
      p_a, guide_area(), 
      widths = c(1/3, 2/3),
      nrow = 3, 
      ncol = 2,
      guides = "collect"
    ) + 
      plot_annotation(
        tag_levels = "A", 
        theme = theme(
          legend.spacing.y = unit(0, unit = "lines"),
          legend.box.margin = margin(12, 0, 0, 0, "pt")
        )
      )
  )
invisible(dev.off())


### Figure 2 =======================================================================================
fdr_meth_value <- rbindlist(lapply(
    X = intersect(
    rownames(data_models[[1, "beta_offsprings"]][[1]]),
    rownames(data_models[[1, "beta_mothers"]][[1]])
  ),
  .form = data_models[[1, "formula"]][[1]],
  .form_m = data_models[[1, "formula_m"]][[1]],
  .beta_o = data_models[[1, "beta_offsprings"]][[1]],
  .beta_m = data_models[[1, "beta_mothers"]][[1]],
  .pheno = data_models[[1, "pheno"]][[1]],
  FUN = function(icpg, .form, .form_m, .beta_o, .beta_m, .pheno) {
    beta_icpg <- c(.beta_o[icpg, ], .beta_m[icpg, ])
    outliers <- c(
      quantile(beta_icpg, probs = 0.25) - 3 * IQR(beta_icpg),
      quantile(beta_icpg, probs = 0.75) + 3 * IQR(beta_icpg)
    )
    lm_mother_resid <- .pheno %>%
      filter(grepl("m", Sample_ID)) %>%
      mutate(Mothers = as.matrix(.beta_m)[icpg, Sample_ID]) %>%
      do(
        add_residuals(
          data = .,
          model = lm(formula = update.formula(.form_m, "Mothers ~ ."), data = .)
        )
      ) %>% 
      select(Sample_ID, Mothers = resid)
    
    lm_pheno <- .pheno %>%
      filter(grepl("c", Sample_ID)) %>%
      mutate(
        Offsprings = .beta_o[icpg, Sample_ID],
        Mothers_raw = .beta_m[icpg, gsub("c", "m", Sample_ID)]
      ) %>% 
      filter(
        between(Offsprings, outliers[1], outliers[2]) &
          between(Mothers_raw, outliers[1], outliers[2])
      ) %>% 
      mutate(MID = gsub("c", "m", Sample_ID)) %>% 
      left_join(y = lm_mother_resid, by = c("MID" = "Sample_ID")) %>% 
      select(Sample_ID, all.vars(.form), Offsprings, Mothers, Mothers_raw) %>%
      drop_na()

    lm_pheno$cpg_id <- icpg
    lm_pheno$Offspring <- ilogit2(lm_pheno$Offsprings)
    lm_pheno$Mother <- ilogit2(lm_pheno$Mothers_raw)
    lm_pheno
  }
))[,
  x_axis := paste0(c("Control", "GDM")[GDM + 1], "<br>(<i style='font-size:8pt'>N = ", .N, "</i>)"),
  by = c("GDM", "cpg_id")
]

fdr_meth_value <- merge(
  x = fdr_meth_value, 
  y = unique(res_om[fdr < 0.05 & term == "GDM", .(cpg_id = CpG, label)]),
  by = "cpg_id"
)[, label := gsub("(.*)<br>\\((.*)\\)", "\\1 - \\2", label)]

top_cpg_plot <- fdr_meth_value[cpg_id %in% res_om[fdr < 0.05 & term == "GDM", CpG],
  list(ggplot = (function(data, label) {
    m_limits <- range(logit2(data[, c("Mother", "Offspring")]), na.rm = TRUE)
    m_breaks <- signif(ilogit2(scales::extended_breaks(n = 5)(m_limits)), 2)
    m_breaks <- c(0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 10, 15) / 100

    p1 <- ggplot(data = data) +
      aes(x = x_axis, y = Offspring, colour = x_axis) +
      geom_boxplot(outlier.shape = NA) +
      geom_quasirandom(shape = 16, size = 0.50, alpha = 0.50) +
      stat_summary(fun = mean, geom = "point", shape = 18, colour = "firebrick", size = 1.5) +
      scale_x_discrete(
        labels = function(x) {
          fifelse(
            grepl("Control", x),
            paste0("<b style='color:#3B528B;'>", x, "</b>"),
            paste0("<b style='color:#5DC863;'>", x, "</b>")
          )
        }
      ) +
      scale_y_continuous(
        trans = logit2_trans(drop0trailing = TRUE), 
        expand = expansion(mult = c(0.05, 0.20)), 
        breaks = m_breaks
      ) +
      scale_colour_viridis_d(begin = 0.25, end = 0.75) +
      labs(
        x = "Offspring", 
        y = paste(unique(label), "<br>Offspring (<i style='font-size:8pt'>&beta;<sub>value</sub> in %</i>)")
      ) +
      theme(legend.position = "none") +
      theme(
        axis.title = element_markdown(size = rel(0.8)),
        axis.text = element_markdown(size = rel(0.6))
      ) +
      coord_cartesian(ylim = range(m_breaks))

    p3 <- ggplot(data = data) +
      aes(x = Mother, y = Offspring, colour = factor(GDM, labels = c("Control", "GDM"))) +
      geom_point(
        shape = 16, 
        size = 0.50, 
        alpha = 0.50, 
        show.legend = FALSE
      ) +
      geom_smooth(
        method = "lm", 
        formula = y ~ x, 
        se = FALSE, 
        size = 0.75
      ) +
      scale_x_continuous(
        trans = logit2_trans(drop0trailing = TRUE), 
        expand = expansion(mult = c(0.05, 0.05)), 
        breaks = m_breaks
      ) +
      scale_y_continuous(
        trans = logit2_trans(drop0trailing = TRUE), 
        expand = expansion(mult = c(0.05, 0.20)), 
        breaks = m_breaks
      ) +
      scale_colour_viridis_d(begin = 0.25, end = 0.75) +
      labs(
        x = paste(unique(label), "<br>Mother (<i style='font-size:8pt'>&beta;<sub>value</sub> in %</i>)"), 
        y = paste(unique(label), "<br>Offspring (<i style='font-size:8pt'>&beta;<sub>value</sub> in %</i>)"),
        colour = NULL
      ) + 
      guides(colour = guide_legend(direction = "horizontal")) +
      theme(
        legend.position = c(0.5, 0.99),
        legend.justification = c("centre", "top"),
        legend.box.just = "centre",
        legend.margin = margin(1.5, 1.5, 1.5, 1.5),
        legend.spacing.x = unit(2.5, "pt"),
        legend.spacing.y = unit(0, "pt"),
        legend.background = element_rect(fill = NA, colour = NA)
      ) +
      theme(
        axis.title = element_markdown(size = rel(0.8)),
        axis.text = element_markdown(size = rel(0.6))
      ) +
      coord_cartesian(xlim = range(m_breaks), ylim = range(m_breaks))
    
    list(p1, p3)
  })(.SD, label)),
  by = "cpg_id"
][, ggplot]

agg_tiff(
  filename = file.path(output_directory, "figure2.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  wrap_plots(top_cpg_plot, nrcol = 2) + 
    plot_annotation(tag_levels = "A")
)
invisible(dev.off())


### Figure 3 =======================================================================================
fdr_meth_value <- rbindlist(lapply(
    X = intersect(
    rownames(data_models[[1, "beta_offsprings"]][[1]]),
    rownames(data_models[[1, "beta_mothers"]][[1]])
  ),
  .form = data_models[[1, "formula"]][[1]],
  .form_m = data_models[[1, "formula_m"]][[1]],
  .beta_o = data_models[[1, "beta_offsprings"]][[1]],
  .beta_m = data_models[[1, "beta_mothers"]][[1]],
  .pheno = data_models[[1, "pheno"]][[1]],
  FUN = function(icpg, .form, .form_m, .beta_o, .beta_m, .pheno) {
    beta_icpg <- c(.beta_o[icpg, ], .beta_m[icpg, ])
    outliers <- c(
      quantile(beta_icpg, probs = 0.25) - 3 * IQR(beta_icpg),
      quantile(beta_icpg, probs = 0.75) + 3 * IQR(beta_icpg)
    )
    lm_mother_resid <- .pheno %>%
      filter(grepl("m", Sample_ID)) %>%
      mutate(Mothers = as.matrix(.beta_m)[icpg, Sample_ID]) %>%
      do(
        add_residuals(
          data = .,
          model = lm(formula = update.formula(.form_m, "Mothers ~ ."), data = .)
        )
      ) %>% 
      select(Sample_ID, Mothers = resid)
    
    lm_pheno <- .pheno %>%
      filter(grepl("c", Sample_ID)) %>%
      mutate(
        Offsprings = .beta_o[icpg, Sample_ID],
        Mothers_raw = .beta_m[icpg, gsub("c", "m", Sample_ID)]
      ) %>% 
      filter(
        between(Offsprings, outliers[1], outliers[2]) &
          between(Mothers_raw, outliers[1], outliers[2])
      ) %>% 
      mutate(MID = gsub("c", "m", Sample_ID)) %>% 
      left_join(y = lm_mother_resid, by = c("MID" = "Sample_ID")) %>% 
      select(Sample_ID, all.vars(.form), Offsprings, Mothers, Mothers_raw) %>%
      drop_na()

    lm_pheno$cpg_id <- icpg
    lm_pheno$Offspring <- ilogit2(lm_pheno$Offsprings)
    lm_pheno$Mother <- ilogit2(lm_pheno$Mothers_raw)
    lm_pheno
  }
))[,
  x_axis := paste0(c("Control", "GDM")[GDM + 1], "<br>(<i style='font-size:8pt'>N = ", .N, "</i>)"),
  by = c("GDM", "cpg_id")
]

fdr_meth_value <- merge(
  x = fdr_meth_value, 
  y = unique(res_om[fdr < 0.05 & term == "Mothers:GDM", .(cpg_id = CpG, label)]),
  by = "cpg_id"
)[, label := gsub("(.*)<br>\\((.*)\\)", "\\1 - \\2", label)]

fdr_meth_value <- fdr_meth_value[!cpg_id %in% res_om[fdr < 0.05 & term == "GDM"]]

agg_tiff(
  filename = file.path(output_directory, "figure3.tiff"), 
  width = 107, 
  height = 107, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  wrap_plots(
    c(
      fdr_meth_value[cpg_id %in% res_om[fdr < 0.05 & term == "Mothers:GDM"][order(pvalue), CpG]][
        order(factor(cpg_id, levels = res_om[fdr < 0.05 & term == "Mothers:GDM"][order(pvalue), CpG])),
        list(ggplot = (function(data, label) {
          p <- ggplot(data = data) +
            aes(x = Mother, y = Offspring, colour = factor(GDM, labels = c("Control", "GDM"))) +
             geom_point(
              shape = 16, 
              size = 0.50, 
              alpha = 0.50, 
              show.legend = FALSE
            ) +
            geom_smooth(
              method = "lm", 
              formula = y ~ x, 
              se = FALSE, 
              size = 0.5
            ) +
            scale_x_continuous(trans = logit2_trans(breaks = TRUE)) +
            scale_y_continuous(trans = logit2_trans(breaks = TRUE)) +
            scale_colour_viridis_d(begin = 0.25, end = 0.75) +
            labs(
              x = paste(unique(label), "<br>Mother (<i style='font-size:8pt'>&beta;<sub>value</sub> in %</i>)"), 
              y = paste(unique(label), "<br>Offspring (<i style='font-size:8pt'>&beta;<sub>value</sub> in %</i>)"),
              colour = NULL
            ) +
            theme(
              axis.title = element_markdown(size = rel(0.8)),
              axis.text = element_markdown(size = rel(0.6))
            )
          list(p)
        })(.SD, label)),
        by = "cpg_id"
      ][, ggplot]
    ),
    ncol = 2,
    byrow = TRUE,
    guides = "collect"
  ) + 
    plot_annotation(tag_levels = "A")
)
invisible(dev.off())


### Supplementary Figure 1 =========================================================================
res_o <- rbind(
  fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_offsprings_simple.csv.gz"
  ))[, type := "RAW"],
  fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_offsprings_simple_cell.csv.gz"
  ))[, type := "CELL"]
)[, cpg_chr := factor(cpg_chr, levels = paste0("chr", 1:22))]

agg_tiff(
  filename = file.path(output_directory, "supplementary_figure1.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(
    data = res_o[order(pvalue)][,
      c("exppval", "labels") := list(
        (1:.N - 0.5) / .N,
        paste0(
          "&lambda;<sub>", type, "</sub> = ", 
          format(median(qnorm(pvalue / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1), digits = 3, nsmall = 3)
        )
      ),
      by = "type"
    ][order(factor(type, levels = c("RAW", "CELL")))][, labels := factor(labels, levels = unique(labels))][order(labels)]
  ) +
    aes(x = exppval, y = pvalue, colour = labels, shape = labels) +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
    geom_point(size = 0.5) +
    scale_x_continuous(trans = pval_trans(), expand = expansion(add = c(0, 0.2)), limits = c(1, NA)) +
    scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
    scale_colour_viridis_d(
      begin = 0.25, 
      end  = 0.75,
      guide = guide_legend(
        direction = "horizontal", 
        title = NULL, 
        ncol = 2, 
        override.aes = list(size = 2)
      )
    ) +
    scale_shape_discrete(
      drop = FALSE, 
      guide = guide_legend(
        direction = "horizontal", 
        title = NULL, 
        ncol = 2
      )
    ) +
    labs(x = "Expected P-value", y = "Observed P-value") +
    coord_cartesian(xlim = c(1, res_o[, min(pvalue)]), ylim = c(1, res_o[, min(pvalue)])) +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      legend.text = element_markdown(),
      legend.position = c(0.5, 1),
      legend.justification = c(0.5, 1.05)
    )
)
invisible(dev.off())


### Supplementary Figure 2 =========================================================================
res_m <- rbind(
  fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_mothers_simple.csv.gz"
  ))[, type := "RAW"],
  fread(here(
    "outputs", "08-ewas_finngedi", "GDM", 
    "EpxGDM_EWAS_DMP_GDM_mothers_simple_cell.csv.gz"
  ))[, type := "CELL"]
)[, cpg_chr := factor(cpg_chr, levels = paste0("chr", 1:22))]

agg_tiff(
  filename = file.path(output_directory, "supplementary_figure2.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(
    data = res_m[order(pvalue)][,
      c("exppval", "labels") := list(
        (1:.N - 0.5) / .N,
        paste0(
          "&lambda;<sub>", type, "</sub> = ", 
          format(median(qnorm(pvalue / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1), digits = 3, nsmall = 3)
        )
      ),
      by = "type"
    ][order(factor(type, levels = c("RAW", "CELL")))][, labels := factor(labels, levels = unique(labels))][order(labels)]
  ) +
    aes(x = exppval, y = pvalue, colour = labels, shape = labels) +
    geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
    geom_point(size = 0.5) +
    scale_x_continuous(trans = pval_trans(), expand = expansion(add = c(0, 0.2)), limits = c(1, NA)) +
    scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
    scale_colour_viridis_d(
      begin = 0.25, 
      end  = 0.75,
      guide = guide_legend(
        direction = "horizontal", 
        title = NULL, 
        ncol = 2, 
        override.aes = list(size = 2)
      )
    ) +
    scale_shape_discrete(
      drop = FALSE, 
      guide = guide_legend(
        direction = "horizontal", 
        title = NULL, 
        ncol = 2
      )
    ) +
    labs(x = "Expected P-value", y = "Observed P-value") +
    coord_cartesian(xlim = c(1, res_m[, min(pvalue)]), ylim = c(1, res_m[, min(pvalue)])) +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      legend.text = element_markdown(),
      legend.position = c(0.5, 1),
      legend.justification = c(0.5, 1.05)
    )
)
invisible(dev.off())


### Supplementary Figure 3 =========================================================================
res_lmm <- fread(here(
  "outputs", "10-ewas_finngedi_lmm", "GDM", 
  "EpxGDM_EWAS_DMP_GDM_lmm_simple_cell.csv.gz"
))[, cpg_chr := factor(cpg_chr, levels = paste0("chr", 1:22))][term == "GDM"]

p_a <- ggplot(
  data = res_lmm[order(pvalue)][,
    c("exppval", "labels") := list(
      (1:.N - 0.5) / .N,
      paste0(
        "&lambda;<sub>gc</sub> = ", 
        format(median(qnorm(pvalue / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1), digits = 3, nsmall = 3)
      )
    ),
    by = "term"
  ]
) +
  aes(x = exppval, y = pvalue, colour = labels, shape = labels) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
  geom_point(size = 0.5) +
  scale_x_continuous(trans = pval_trans(), expand = expansion(add = c(0, 0.2)), limits = c(1, NA)) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_viridis_d(
    begin = 0.50, 
    end  = 0.50,
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      ncol = 2, 
      override.aes = list(size = 2)
    )
  ) +
  scale_shape_discrete(
    drop = FALSE, 
    guide = guide_legend(
      direction = "horizontal", 
      title = NULL, 
      ncol = 2
    )
  ) +
  labs(x = "Expected P-value", y = "Observed P-value") +
  coord_cartesian(xlim = c(1, res_lmm[, min(pvalue)]), ylim = c(1, res_lmm[, min(pvalue)])) +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown(),
    legend.position = c(0.5, 1),
    legend.justification = c(0.5, 1.05)
  )

p_b <- ggplot(data = res_lmm) +
  aes(x = estimate, y = pvalue, colour = abs(estimate)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 0.5) + 
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_viridis_c(trans = "sqrt", limits = c(0, NA), guide = "none", option = "viridis") +
  labs(x = "Estimate (M<sub>value</sub>)", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

p_c <- ggplot(data = fortify.manhattan(res_lmm[, -"group"], "cpg_pos", "pvalue", "cpg_chr")) +
  aes(x = x, y = y, colour = group) +
  geom_point(size = 0.5) +
  scale_x_continuous(breaks = 1:24, labels = c(1:22, "X", "Y"), expand = expansion(add = 0.25)) +
  scale_y_continuous(trans = pval_trans(), expand = expansion(mult = c(0, 0.2)), limits = c(1, NA)) +
  scale_colour_manual(values = rep(viridis_pal(begin = 1/4, end = 3/4, option = "viridis")(2), 12), guide = "none") +
  labs(x = "Chromosome", y = "P-value") +
  theme(
    axis.title = element_markdown(size = rel(0.8)),
    axis.text = element_markdown(size = rel(0.6)),
    legend.text = element_markdown()
  )

agg_tiff(
  filename = file.path(output_directory, "supplementary_figure3.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  wrap_plots(
    p_a, p_b, p_c,
    design = "AB\nCC"
  ) + 
    plot_annotation(
      tag_levels = "A", 
      theme = theme(
        legend.spacing.y = unit(0, unit = "lines"),
        legend.box.margin = margin(12, 0, 0, 0, "pt")
      )
    )
)
invisible(dev.off())


### Supplementary Figure 4 =========================================================================
# gtex_tsv <- "https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip"
# download.file(
#   url = gtex_tsv, 
#   destfile = file.path(tempdir(), basename(gtex_tsv))
# )
# unzip(zipfile = file.path(tempdir(), basename(gtex_tsv)), exdir = tempdir())
# gtex_data <- fread(file.path(tempdir(), gsub(".zip$", "", basename(gtex_tsv))))

gtex_tsv <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
download.file(
  url = gtex_tsv, 
  destfile = file.path(tempdir(), basename(gtex_tsv))
)
gtex_data <- melt(
  data = fread(file.path(tempdir(), basename(gtex_tsv))), 
  id.vars = c("Name", "Description"), 
  variable.name = "Tissue", 
  value.name = "TPM"
)[, group := gsub("(.*) - .*", "\\1", Tissue)]

agg_tiff(
  filename = file.path(output_directory, "supplementary_figure4.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(data = gtex_data[like(Description, "^TFCP2$")]) +
    aes(x = fct_reorder(Tissue, TPM), y = TPM, fill = group) +
    geom_col() +
    scale_x_discrete(labels = toupper, guide = guide_axis(angle = 45)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_viridis_d(begin = 0, end = 1) +
    labs(x = NULL, y = "TPM (Median)") +
    guides(fill = "none") +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      axis.text.x = element_markdown(size = rel(0.75)), 
      plot.margin = margin(5.5, 5.5, 5.5, 5.5 * 1, "pt")
    )
)
invisible(dev.off())


### Supplementary Figure 5 =========================================================================
agg_tiff(
  filename = file.path(output_directory, "supplementary_figure5.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(data = gtex_data[like(Description, "^DLGAP2$")]) +
    aes(x = fct_reorder(Tissue, TPM), y = TPM, fill = group) +
    geom_col() +
    scale_x_discrete(labels = toupper, guide = guide_axis(angle = 45)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_viridis_d(begin = 0, end = 1) +
    labs(x = NULL, y = "TPM (Median)") +
    guides(fill = "none") +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      axis.text.x = element_markdown(size = rel(0.75)), 
      plot.margin = margin(5.5, 5.5, 5.5, 5.5 * 6, "pt")
    )
)
invisible(dev.off())


### Supplementary Figure 6 =========================================================================
agg_tiff(
  filename = file.path(output_directory, "supplementary_figure6.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(data = gtex_data[like(Description, "^FAM13A$")]) +
    aes(x = fct_reorder(Tissue, TPM), y = TPM, fill = group) +
    geom_col() +
    scale_x_discrete(labels = toupper, guide = guide_axis(angle = 45)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_viridis_d(begin = 0, end = 1) +
    labs(x = NULL, y = "TPM (Median)") +
    guides(fill = "none") +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      axis.text.x = element_markdown(size = rel(0.75)), 
      plot.margin = margin(5.5, 5.5, 5.5, 5.5 * 10, "pt")
    )
)
invisible(dev.off())


### Supplementary Figure 7 =========================================================================
agg_tiff(
  filename = file.path(output_directory, "supplementary_figure7.tiff"), 
  width = 107, 
  height = 80.25, 
  units = "mm",
  res = 300,
  pointsize = 12,
  scaling = 0.66
)
print(
  ggplot(data = gtex_data[like(Description, "^HIST1H3E$")]) +
    aes(x = fct_reorder(Tissue, TPM), y = TPM, fill = group) +
    geom_col() +
    scale_x_discrete(labels = toupper, guide = guide_axis(angle = 45)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_viridis_d(begin = 0, end = 1) +
    labs(x = NULL, y = "TPM (Median)") +
    guides(fill = "none") +
    theme(
      axis.title = element_markdown(size = rel(0.8)),
      axis.text = element_markdown(size = rel(0.6)),
      axis.text.x = element_markdown(size = rel(0.75)), 
      plot.margin = margin(5.5, 5.5, 5.5, 5.5 * 12, "pt")
    )
)
invisible(dev.off())


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
