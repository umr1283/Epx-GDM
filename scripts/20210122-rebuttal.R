### Project Setup ==================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", basename(here()))
output_directory <- here("outputs", "20210122-rebuttal")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

data_directory <- normalizePath("/disks/DATA/Projects/EpxGDM/FINNGEDI_TUBINGEN/QC")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(ragg)
  library(stats)
  library(ggdendro)
  library(ggplot2)
  library(ggtext)
  library(scales)
  library(patchwork)
  library(Hmisc)
  library(data.table)
  library(flashpcaR)
  library(readxl)
  library(gtsummary)
  library(gt)
})


### Tables and Figures Theme =======================================================================
options(
  ggplot2.discrete.colour = function(...) scale_colour_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.colour = function(...) scale_colour_viridis_c(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.fill = function(...) scale_fill_viridis_c(..., begin = 0.15, end = 0.85)
)
theme_set(
  theme_minimal(base_family = "Times") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.8)),
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown()
    )
)


### Functions ======================================================================================


### Get Data =======================================================================================
sample_sheet <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes_SNPs.csv")
)[
  (!qc_gender_discrepancy) & call_rate >= 0.99 & Project %in% "FinnGeDi"
][
  j = .SD[call_rate == max(call_rate)],
  by = "Sample_IID"
]

beta_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues_SNPs.csv.gz"), 
  header = TRUE
)

common_samples <- intersect(sample_sheet[["Sample_ID"]], colnames(beta_matrix))
sample_sheet <- sample_sheet[Sample_ID %in% common_samples]
beta_matrix <- beta_matrix[j = .SD, .SDcols = c("cpg_id", common_samples)]

manifest_hg19 <- get(data("EPIC.manifest.pop.hg19", package = "ChAMPdata"))
which_population <- which(manifest_hg19[, paste("MASK_general", "FIN", sep = "_")])
beta_matrix_fin <- beta_matrix[cpg_id %in% rownames(manifest_hg19)[which_population]]
setnames(beta_matrix_fin, old = names(beta_matrix_fin), new = gsub("_r", "", names(beta_matrix_fin)))


### Pair Matching ==================================================================================
pairs <- table(gsub("[mc]$", "", colnames(beta_matrix_fin)[-1]))
good_pairs <- names(pairs)[pairs == 2]

best_diff <- optimise(
  f = function(x) {
    fin_snps <- which(apply(apply(beta_matrix_fin[, -1], 1, range), 2, diff) > x)
    samples_corr <- rcorr(as.matrix(beta_matrix_fin[fin_snps, -1]), type = "pearson")
    
    pairs_corr <- melt(
      data = as.data.table(samples_corr$r, keep.rownames = "x"), 
      id.vars = "x",
      measure.vars = names(samples_corr$r), 
      variable.name = "y",
      value.name = "pearson"
    )[
      gsub("[cm]$", "", x) %in% good_pairs & gsub("[cm]$", "", y) %in% good_pairs
    ][
      grepl("m", x) & grepl("c", y)
    ][
      j = rank := paste0(rank(1 - pearson), "/", max(rank(1 - pearson))), 
      by = "x"
    ]
    pairs_matched <- pairs_corr[gsub("[cm]$", "", x) == gsub("[cm]$", "", y)][order(-rank, as.numeric(gsub("[cm]$", "", x)))]
    
    length(unique(pairs_matched$rank)) - 1 + x
  }, 
  interval = c(0.65, 0.85)
)$minimum

fin_snps <- which(apply(apply(beta_matrix_fin[, -1], 1, range), 2, diff) > best_diff)
samples_corr <- rcorr(as.matrix(beta_matrix_fin[fin_snps, -1]), type = "pearson")

pairs_corr <- melt(
  data = as.data.table(samples_corr$r, keep.rownames = "x"), 
  id.vars = "x",
  measure.vars = names(samples_corr$r), 
  variable.name = "y",
  value.name = "pearson"
)[
  gsub("[cm]$", "", x) %in% good_pairs & gsub("[cm]$", "", y) %in% good_pairs
][
  grepl("m", x) & grepl("c", y)
][
  j = rank := paste0(rank(1 - pearson), "/", max(rank(1 - pearson))), 
  by = "x"
]
pairs_matched <- pairs_corr[gsub("[cm]$", "", x) == gsub("[cm]$", "", y)][order(-rank, as.numeric(gsub("[cm]$", "", x)))]

fwrite(
  x = pairs_matched, 
  file = file.path(output_directory, "pair_matching.csv")
)

agg_png(filename = file.path(output_directory, "pair_matching.png"), width = 16, height = 12, units = "cm", res = 300)
  print(
    ggplot(data = pairs_matched) +
      aes(x = pearson) +
      geom_density() +
      geom_rug() +
      labs(
        x = "Pearson Correlation Between Mother/Offspring", 
        y = "Density",
        title = "Pair Matching Between Mother And Offspring",
        subtitle = paste(
          "Using a total of", format(length(fin_snps), big.mark = ","), 
          "likely CpG-SNPs<sup>&dagger;</sup> for the Finnish population"
        ),
        caption = paste0("<sup>&dagger;</sup>CpG-SNPs defined with an observed change (Max-Min) of ", round(best_diff, 5), ".")
      ) +
      coord_cartesian(xlim = c(0, 1)) +
      scale_x_continuous(labels = percent_format(suffix = " %"))
  )
invisible(dev.off())


### PCA ============================================================================================
pairs <- table(gsub("[mc]$", "", colnames(beta_matrix)[-1]))
good_pairs <- unlist(lapply(c("m", "c"), function(x) paste0(names(pairs)[pairs == 2], x)), use.names = FALSE)


## All CpGs ----------------------------------------------------------------------------------------
pca_methylation <- `rownames<-`(
  as.matrix(beta_matrix[j = .SD, .SDcols = c(good_pairs)]), 
  beta_matrix[["cpg_id"]]
)

pca_all_cpgs <- wrap_plots(
  lapply(X = c("m", "c"), FUN = function(x) {
    pca_res <- flashpca(X = t(pca_methylation[, grepl(x, colnames(pca_methylation))]), stand = "sd", ndim = 10)
    
    pca_dfxy <- merge(
      x = sample_sheet,
      y = as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_IID"),
      by = "Sample_IID"
    )
    pca_dfxy[j = GDM := factor(GDM, labels = c("Control", "GDM"))]
    pc_names <- sprintf("PC%02d (%s %%)", seq_along(pca_res[["pve"]]), round(pca_res[["pve"]] * 100, 2))
    setnames(
      x = pca_dfxy, 
      old = paste0("V", seq_along(pca_res[["pve"]])), 
      new = pc_names
    )
  
    ivar <- "GDM"
    p <- apply(
      X = combn(pc_names[1:3], 2),
      MARGIN = 2,
      FUN = function(x) {
        ggplot(data = pca_dfxy[, .SD, .SDcols = c(ivar, x)]) +
          aes(x = .data[[x[1]]], y = .data[[x[2]]], colour = .data[[ivar]]) +
          geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
          geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
          geom_point(na.rm = TRUE) +
          {
            if (is.numeric(pca_dfxy[[ivar]])) {
              scale_colour_viridis_c(
                name = NULL,
                begin = 0,
                end = 0.75
              )
            } else {
              list(
                stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE),
                scale_colour_viridis_d(
                  name = NULL,
                  begin = if (pca_dfxy[, uniqueN(.SD), .SDcols = ivar] == 2) 0.25 else 0,
                  end = 0.75, 
                  guide = guide_legend(override.aes = list(size = 4))
                ),
                if (length(unique(pca_dfxy[[ivar]])) > 10) {
                  theme(legend.position = "none")
                } else {
                  NULL
                }
              )
            }
          }
      }
    )
    wrap_plots(p, ncol = 3, tag_level = "new", guides = "collect")
  }), 
  nrow = 2, 
  guides = "collect"
) + 
  plot_annotation(
    title = "Structure Detection Using Principal Component Analysis: All CpGs",
    subtitle = "With **A)** Mothers and **B)** Offspring",
    tag_levels = c("A", 1)
  )

agg_png(
  filename = file.path(output_directory, "pca_all_cpgs.png"),
  width = 16 / 2.54, height = 12 / 2.54, units = "in", res = 300, scaling = 0.75
)
  print(pca_all_cpgs)
invisible(dev.off())

## FIN CpGs ----------------------------------------------------------------------------------------
fin_snps <- which(apply(apply(beta_matrix_fin[j = .SD, .SDcols = c(good_pairs)], 1, range), 2, diff) > 0.75)

pca_methylation <- `rownames<-`(
  as.matrix(beta_matrix_fin[j = .SD, .SDcols = c(good_pairs)]), 
  beta_matrix_fin[["cpg_id"]]
)

pca_fin_cpgs <- wrap_plots(
  lapply(X = c("m", "c"), FUN = function(x) {
    pca_res <- flashpca(X = t(pca_methylation[, grepl(x, colnames(pca_methylation))]), stand = "sd", ndim = 10)
    
    pca_dfxy <- merge(
      x = sample_sheet,
      y = as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_IID"),
      by = "Sample_IID"
    )
    pca_dfxy[j = GDM := factor(GDM, labels = c("Control", "GDM"))]
    pc_names <- sprintf("PC%02d (%s %%)", seq_along(pca_res[["pve"]]), round(pca_res[["pve"]] * 100, 2))
    setnames(
      x = pca_dfxy, 
      old = paste0("V", seq_along(pca_res[["pve"]])), 
      new = pc_names
    )
  
    ivar <- "GDM"
    p <- apply(
      X = combn(pc_names[1:3], 2),
      MARGIN = 2,
      FUN = function(x) {
        ggplot(data = pca_dfxy[, .SD, .SDcols = c(ivar, x)]) +
          aes(x = .data[[x[1]]], y = .data[[x[2]]], colour = .data[[ivar]]) +
          geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
          geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
          geom_point(na.rm = TRUE) +
          {
            if (is.numeric(pca_dfxy[[ivar]])) {
              scale_colour_viridis_c(
                name = NULL,
                begin = 0,
                end = 0.75
              )
            } else {
              list(
                stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE),
                scale_colour_viridis_d(
                  name = NULL,
                  begin = if (pca_dfxy[, uniqueN(.SD), .SDcols = ivar] == 2) 0.25 else 0,
                  end = 0.75, 
                  guide = guide_legend(override.aes = list(size = 4))
                ),
                if (length(unique(pca_dfxy[[ivar]])) > 10) {
                  theme(legend.position = "none")
                } else {
                  NULL
                }
              )
            }
          }
      }
    )
    wrap_plots(p, ncol = 3, tag_level = "new", guides = "collect")
  }), 
  nrow = 2, 
  guides = "collect"
) + 
  plot_annotation(
    title = "Structure Detection Using Principal Component Analysis: Finnish CpG-SNPs",
    subtitle = "With **A)** Mothers and **B)** Offspring",
    tag_levels = c("A", 1)
  )

agg_png(
  filename = file.path(output_directory, "pca_finnish_cpgs.png"),
  width = 16 / 2.54, height = 12 / 2.54, units = "in", res = 300, scaling = 0.75
)
  print(pca_fin_cpgs)
invisible(dev.off())

## Check Variables ---------------------------------------------------------------------------------
pheno <- setDT(
  read_excel(here("docs", "FinnGedi_clinical_data_2018_10_05", "FinnGeDi_selection_EpxGDM_2018_10_06_TO_TOBY.xlsx"))
)
tbl_summary(
  data = pheno, 
  by = "GDM", 
  include = setdiff(names(which(sapply(pheno, uniqueN) >= 2)), c("ID", "DNAtrio", "DNAMother1", "DNAMother2", "DNAchild"))
) %>% 
  add_p() %>% 
  as_gt() %>% 
  tab_options(
    column_labels.font.weight = "bold", 
    data_row.padding = 0.025,
    # container.width = px(1200), 
    table.font.size = "8pt"
  ) %>% 
  gtsave(filename = file.path(output_directory, "pheno_asso.html"))

all_vars <- c(
  "GDM", "mumage", "mumagecategorical", "prebmi", "BMIcategorical", 
  "GDMwk", "ogtt", "Insulin", "GDMpharmacologytreated", "Metformin", 
  "primipara", "previousdeliveries", "previousdeliverycategorical", 
  "gestwk", "sex", "ModeOfDelivery", "BabyRespSupport", "BabyGluc1", 
  "BabyGluc2", "BabyGlucLowesth", "BabyGlucLowestmin", "BabyGluc2ndLowesth", 
  "BabyGluc2ndLowestmin", "BabySugar", "bwt", "blt", "bhc", "BWSDSfi13", 
  "BLTSDSfi13", "BHCSDSfi13", "lga", "sga", "apgar1", "apgar5", 
  "BabyRespiratorySupport", "mumwtgain", "mumwtgainnormal", "glucabnn", 
  "glucearlyabn", "gluc0hfirstabn", "gluc1hfirstabn", "gluc2hfirstabn", 
  "gluc0hCtrl24wk", "gluc1hCtrl24wk", "gluc2hCtrl24wk", "PE", "CH", 
  "GH", "mgdm", "mgdm2", "mfatherdm", "mmotherdm", "dgkotiseuranta", 
  "hospital", "ses_4lk", "education", "GDMearlierpreg", "BPsHighestEarly", 
  "BPdHighestEarly", "BPsHighestLate", "BPdHighestLate"
)

ggplot_association <- lapply(X = c("m", "c"), FUN = function(x) {
  pca_res <- flashpca(X = t(pca_methylation[, grepl(x, colnames(pca_methylation))]), stand = "sd", ndim = 10)
  
  pca_dfxy <- merge(
    x = merge(
      x = sample_sheet[j = c("Sample_IID", "ID", "Sample_Plate", "Sentrix_ID", "Sentrix_Position")],
      y = as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_IID"),
      by = "Sample_IID"
    ),
    y = pheno[j = .SD, .SDcols = c("ID", all_vars)],
    by = "ID"
  )
  pc_names <- sprintf("PC%02d (%s %%)", seq_along(pca_res[["pve"]]), round(pca_res[["pve"]] * 100, 2))
  setnames(
    x = pca_dfxy, 
    old = paste0("V", seq_along(pca_res[["pve"]])), 
    new = pc_names
  )
  
  asso_dt <- melt(
    data = pca_dfxy,
    measure.vars = pc_names,
    variable.name = "pc",
    value.name = "values"
  )[pc %in% pc_names[1:5]][,
    {
      m <- model.matrix(
        object = as.formula(
          object = paste0("values ~ ", paste(all_vars, collapse = " + "))
        ),
        data = .SD
      )
  
      if (qr(m)$rank == ncol(m)) {
        out <- as.data.table(
          anova(
            lm(
              formula = as.formula(
                object = paste0("values ~ ", paste(all_vars, collapse = " + "))
              ),
              data = .SD
            )
          ),
          keep.rownames = "term"
        )[term != "Residuals"]
      } else {
        out <- rbindlist(
          lapply(X = all_vars, .data = .SD, FUN = function(.x, .data) {
            as.data.table(
              anova(
                lm(
                  formula = as.formula(paste0("values ~ ", .x)),
                  data = .SD
                )
              ),
              keep.rownames = "term"
            )[term != "Residuals"]
          })
        )
      }
      out[, full_rank := qr(m)$rank == ncol(m)]
    },
    by = "pc"
  ]
  
  y_levels <- setorderv(
    x = dcast(
      data = asso_dt[j = list(pc, term, `Pr(>F)` = fifelse(`Pr(>F)` <= 0.1, `Pr(>F)`, NA_real_))], 
      formula = term ~ pc, 
      value.var = "Pr(>F)"
    ), 
    cols = levels(asso_dt[["pc"]])[1:5], 
    order = -1
  )[["term"]]
  
  asso_dt[
    j = term := factor(x = term, levels = y_levels)
  ]
  
  p_association <- ggplot(data = asso_dt) +
    aes(
      x = factor(.data[["pc"]]), 
      y = term,
      fill = .data[["Pr(>F)"]]
    ) +
    geom_tile(colour = "white", na.rm = TRUE) +
    geom_richtext(
      mapping = aes(
        label = gsub(
          pattern = "(.*)e([-+]*)0*(.*)",
          replacement = "\\1 &times; 10<sup>\\2\\3</sup>",
          x = format(.data[["Pr(>F)"]], digits = 2, nsmall = 2, scientific = TRUE)
        )
      ),
      colour = "white",
      fill = NA,
      label.colour = NA,
      size = 1.25,
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(na.value = "grey85", end = 0.75, limits = c(0, 0.1)) +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(expand = c(0, 0), labels = function(x) gsub(" \\(", "<br>(", x)) +
    scale_y_discrete(labels = toupper, expand = c(0, 0)) +
    labs(
      x = "Principal Component",
      y = "Variable",
      # title = paste(
      #   "Association Tests Between Variables And", c("m" = "Mothers", "c" = "Offspring")[x], "Principal Components"
      # ),
      # caption = "<sup>&dagger;</sup>Each variable is tested against principal components using ANOVA.",
      fill = "P-Value<i><sup>&dagger;</sup></i>"
    ) +
    theme(legend.title = element_markdown())
  
  # agg_png(
  #   filename = file.path(output_directory, paste0("pca_asso_", c("m" = "mothers", "c" = "offspring")[x], ".png")),
  #   width = 16,
  #   height = 24.7,
  #   units = "cm",
  #   res = 300, 
  #   scaling = 0.75
  # )
  #   print(p_association)
  # invisible(dev.off())
  
  p_association
})

p_association <- wrap_plots(
  lapply(X = ggplot_association, FUN = function(.gg) {
    .gg$data$term <- factor(.gg$data$term, levels = sort(as.character(unique(.gg$data$term)), decreasing = TRUE))
    .gg
  }), 
  ncol = 2, 
  guides = "collect"
) + 
  plot_annotation(
    title = "Association Tests Between Variables And Principal Components",
    subtitle = "With **A)** Mothers and **B)** Offspring",
    caption = "<sup>&dagger;</sup>Each variable is tested against principal components using ANOVA.",
    tag_levels = "A"
  )

agg_png(
  filename = file.path(output_directory, "pca_asso.png"),
  width = 24.7, height = 16, units = "cm", res = 300, scaling = 0.75
)
  print(p_association)
invisible(dev.off())


### Adjustment for mumwtgain =======================================================================
results_files <- c(
  here("outputs/09-ewas_finngedi_h/GDM/EpxGDM_EWAS_DMP_GDM_mothers_offsprings_simple_cell.csv.gz"),
  list.files(path = here("outputs/16-ewas_finngedi_h_rebuttal/GDM"), pattern = ".csv.gz$", full.names = TRUE)
)

results <- rbindlist(lapply(X = results_files, FUN = function(.file) {
  model <- gsub(".*_cell_*(.*).csv.gz", "\\1", basename(.file))
  if (model == "") model <- "default"
  fread(.file)[j = rebuttal_model := model][]
}))

top_cpgs <- results[
  i = CpG %in% results[fdr <= 0.05, CpG],
  j = list(
    CpG,
    Gene = fifelse(UCSC_RefGene_Name != "", UCSC_RefGene_Name, NA_character_), 
    term, estimate, se, pvalue, fdr, 
    rebuttal_model = gsub("_", "", rebuttal_model)
  )
][j = min_p := min(pvalue), by = "CpG"][order(min_p)][
  j = `:=`(CpG = factor(CpG, levels = unique(CpG)), min_p = NULL)
][
  j = term := gsub("Mother:GDM", "GDM × Mother", gsub("Mothers", "Mother", term))
][
  j = term := factor(term, levels = c("GDM", "GDM × Mother", "Mother"))
][
  order(CpG, term, rebuttal_model)
]

dt_summary <- dcast(
  data = top_cpgs, 
  formula = ... ~ rebuttal_model, 
  value.var = names(which(sapply(top_cpgs, is.numeric)))
)

fmt_spanner <- function(.data, label) {
  .data %>% 
    fmt_number(
      columns = paste0("fdr_", label),
      rows = dt_summary[[paste0("fdr_", label)]] >= 0.01, 
      n_sigfig = 3
    ) %>% 
    fmt_scientific(
      columns = paste0("fdr_", label),
      rows = dt_summary[[paste0("fdr_", label)]] < 0.01,
      decimals = 2
    ) %>% 
    fmt_number(
      columns = paste0("pvalue_", label),
      rows = dt_summary[[paste0("pvalue_", label)]] >= 0.01, 
      n_sigfig = 3
    ) %>% 
    fmt_scientific(
      columns = paste0("pvalue_", label),
      rows = dt_summary[[paste0("pvalue_", label)]] < 0.01,
      decimals = 2
    ) %>% 
    cols_merge_uncert(col_val = paste0("estimate_", label), col_uncert = paste0("se_", label)) %>% 
    tab_spanner(
      label = c(
        "default" = "Default", 
        "mumwtgainprimipara" = "Default + Mother Weight Gain + Primipara", 
        "mumwtgain" = "Default + Mother Weight Gain", 
        "primipara" = "Default + Primipara"
      )[label], 
      columns = matches(paste0("_", label, "$")), 
      gather = TRUE
    )
}

gt(data = dt_summary, auto_align = "center") %>% 
  fmt_number(columns = matches("estimate|se"), n_sigfig = 3) %>% 
  fmt_missing(columns = "Gene") %>% 
  fmt(columns = "Gene", fns = function(x) fifelse(is.na(x), NA_character_, paste0("<i>", x, "</i>"))) %>% 
  (function(.data, .label) {
    for (x in .label) {.data <- fmt_spanner(.data, x) } 
    .data 
  })(unique(top_cpgs[["rebuttal_model"]])) %>% 
  cols_label(
    .list = `names<-`(
      gsub(
        paste(paste0("_", unique(top_cpgs[["rebuttal_model"]]), "$"), collapse = "|"), 
        "", 
        names(dt_summary)
      ),
      names(dt_summary)
    )
  ) %>% 
  opt_all_caps() %>% 
  opt_row_striping() %>% 
  tab_options(
    column_labels.font.weight = "bold", 
    data_row.padding = 0.025,
    # container.width = px(1200), 
    table.font.size = "8pt"
  ) %>% 
  gtsave(filename = file.path(output_directory, "new_ewas.html"))


### Check Probes Quality ===========================================================================
if (!file.exists(file.path(output_directory, "detection_p.csv.gz"))) {
  rgSet <- readRDS(here("outputs/02-qc_idats/EPIC_idats.rds"))
  data_detP <- minfi::detectionP(rgSet$rgset)
  data_detP[is.na(data_detP)] <- 1
  data_detP[data_detP == 0] <- min(data_detP[data_detP != 0])
  fwrite(
    x = as.data.table(x = data_detP, keep.rownames = "cpg_id"), 
    file = file.path(output_directory, "detection_p.csv.gz")
  )
} else {
  detection_pvalue <- fread(file.path(output_directory, "detection_p.csv.gz"))
}

agg_png(
  filename = file.path(output_directory, "probes_quality.png"),
  width = 16, height = 12, units = "cm", res = 300, scaling = 0.75
)
  print(wrap_plots(
    melt(
      data = detection_pvalue[
        cpg_id %in% c("cg22790973", "cg03456133", "cg23355087", "cg19107264", "cg17065901", "cg24440941", "cg20002843")
      ], 
      id.vars = "cpg_id"
    )[
      j = list(list(
        ggplot(data = .SD) +
          aes(x = -log10(value)) +
          geom_density() +
          geom_rug() +
          labs(x = "Detection P-value (-log<sub>10</sub>)", y = "Density", title = cpg_id) +
          theme(legend.position = "none", panel.grid.minor.y = element_blank())
    
      )),
      by = "cpg_id"  
    ][["V1"]], 
    ncol = 3
  ))
invisible(dev.off())


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
