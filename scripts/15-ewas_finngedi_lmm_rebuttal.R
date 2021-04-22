# Docker Container R403
### Load packages ==================================================================================
library(here)

suppressPackageStartupMessages({
  library(DMRcate)
  library(tibble)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tidyr)
  library(dplyr)
  library(parallel)
  library(glue)
  library(purrr)
  library(modelr)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(patchwork)
  library(ggtext)
  library(forcats)
  library(gt)
  library(readxl)
})


### Source functions ===============================================================================
#' pval_trans
#'
#' @return trans
#' @export
pval_trans <- function(alpha = NULL, md = FALSE, prefix = FALSE) {
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
      g <- function(x) {
        gsub(
          "1 %*% ", "",
          gsub(
            "(.*)e([-+]*)0*(.*)", "\\1 %*% 10^\\2\\3",
            gsub(
             "1e+00", "1",
              scales::scientific(x),
             fixed = TRUE
            )),
          fixed = TRUE
        )
      }
      highlight_alpha <- function(x, md = FALSE, prefix = FALSE) {
        if (md & nchar(system.file(package = "ggtext")) != 0) {
          prefix_text <- if (prefix) "&alpha; = " else ""
          out <- paste0(
            "<b style = 'color:firebrick2;'>", 
            gsub("(.*) [%][*][%] .*\\^(.*)", paste(prefix_text, "\\1 &times; 10<sup>\\2</sup>"), g(x)), 
            "</b>"
          )
        } else {
          prefix_text <- if (prefix) "alpha ==" else ""
          out <- parse(text = paste(prefix_text, g(x)))
        }
        
        out
      }
      
      if (!is.null(alpha)) {
        ifelse(
          test = scales::scientific(x) == scales::scientific(alpha), 
          yes = highlight_alpha(x, md, prefix), 
          no = g(x)
        )
      } else {
        parse(text = g(x))
      }
    })
  )
}

#' stat_manhattan
#'
#' @param mapping []
#' @param data []
#' @param geom []
#' @param position []
#' @param na.rm []
#' @param show.legend []
#' @param inherit.aes []
#' @param shape []
#' @param fill []
#' @param ... []
#'
#' @return theme
#' @export
stat_manhattan <- function(
  mapping = NULL,
  data = NULL,
  geom = "point",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  shape = 16,
  size = 0.5,
  fill = NA,
  ...
) {
  ggplot2::layer(
    stat = StatManhattan,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      shape = shape,
      fill = fill,
      size = size,
      ...
    )
  )
}


#' @rdname stat_manhattan
#' @export
geom_manhattan <- function(
  mapping = NULL,
  data = NULL,
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  shape = 16,
  fill = NA,
  size = 0.5,
  ...
) {
  list(
    ggplot2::layer(
      stat = StatManhattan,
      data = data,
      mapping = mapping,
      geom = "point",
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        shape = shape,
        fill = fill,
        size = size,
        ...
      )
    ),
    ggplot2::scale_x_continuous(
      breaks = 1:24,
      labels = c(1:22, "X", "Y"),
      expand = ggplot2::expansion(add = 0.25)
    ),
    ggplot2::scale_y_continuous(
      trans = "pval",
      expand = ggplot2::expansion(mult = c(0, 0.10)),
      limits = c(1, NA)
    ),
    ggplot2::scale_colour_manual(values = rep(scales::viridis_pal(begin = 1/4, end = 3/4)(2), 12)),
    ggplot2::theme(
      panel.grid.minor.x = ggplot2::element_blank()
    ),
    ggplot2::labs(colour = "Chromosome", x = "Chromosome", y = "P-Value")
  )
}

#' @rdname stat_manhattan
#'
#' @format NULL
#' @usage NULL
#'
#' @export
StatManhattan <- ggplot2::ggproto("StatManhattan", ggplot2::Stat,
  required_aes = c("x", "y", "colour"),
  default_aes = ggplot2::aes(group = ggplot2::stat(colour)),
  setup_data = function(data, params) {
    map_chro <- c(seq(22), "X", "Y", "X", "Y") 
    names(map_chro) <- c(seq(24), "X", "Y")
      
    data %>%
      dplyr::mutate(
        x_chr = map_chro[colour],
        x_chr = factor(x_chr, levels = unique(map_chro)),
        colour = x_chr,
        x_pos = as.integer(x),
        y_pval = as.numeric(y)
      ) %>%
      dplyr::group_by(x_chr) %>%
      dplyr::arrange(x_pos) %>%
      dplyr::mutate(x_pos = scales::rescale(x = x_pos, to = c(-0.4, 0.4))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(x_pos = x_pos + as.integer(x_chr)) %>%
      dplyr::select(-x, -y, -group) %>%
      dplyr::rename(x = x_pos, y = y_pval, group = x_chr)
  },
  compute_panel = function(data, scales, params) {
    data
  }
)

fortify.manhattan <- function(data, x, y, group) {
  map_chro <- c(seq(22), "X", "Y", "X", "Y") 
  names(map_chro) <- c(seq(24), "X", "Y")
    
  data %>%
    dplyr::mutate(
      x_chr = map_chro[.data[[group]]],
      x_chr = factor(x_chr, levels = unique(map_chro)),
      colour = x_chr,
      x_pos = as.integer(.data[[x]]),
      y_pval = as.numeric(.data[[y]])
    ) %>%
    dplyr::group_by(x_chr) %>%
    dplyr::arrange(x_pos) %>%
    dplyr::mutate(x_pos = scales::rescale(x = x_pos, to = c(-0.4, 0.4))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x_pos = x_pos + as.integer(x_chr)) %>%
    dplyr::select(-x, -y, -group) %>%
    dplyr::rename(x = x_pos, y = y_pval, group = x_chr)
}


### Environment ====================================================================================
project_name <- gsub("(.*)_.*", "\\1", basename(here()))
output_directory <- here("outputs", "15-ewas_finngedi_lmm_rebuttal")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0755")

data_directory <- normalizePath("/disks/DATA/Projects/EpxGDM/FINNGEDI_TUBINGEN/QC")

n_cores <- 40


### Define theme ===================================================================================
# options(ggplot2.continuous.colour = "viridis", ggplot2.continuous.fill = "viridis")
theme_set(theme_light())


### Traits to be analysed ==========================================================================
traits <- c("Gestational Diabetes Mellitus" = "GDM")


### Data ===========================================================================================
beta_matrix <- fread(file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues.csv.gz")) %>%
  column_to_rownames("cpg_id") %>%
  as.matrix() %>%
  (function(x) log2(x) - log2(1 - x))() # logit2

# ethnicity_pcs <- read_csv(
#   file = file.path(data_directory, "Omni2.5", glue("{project_name}_ethnicity.csv")),
#   col_types = cols(
#     .default = col_double(),
#     sample = col_character(),
#     cohort = col_character(),
#     pop = col_character(),
#     super_pop = col_character()
#   )
# ) %>%
#   separate(col = "sample", into = c("family_id", "sample_id"), convert = FALSE)

phenotype_matrix <- full_join(
  x = read_csv(file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes.csv")),
  y = here(
    "docs", "FinnGedi_clinical_data_2018_10_05", "FinnGeDi_selection_EpxGDM_2018_10_06_TO_TOBY.xlsx"
  ) %>% 
    read_excel() %>% 
    select(ID, primipara, mumwtgain),
  by = "ID"
)


# phenotype_pcs <- left_join(x = phenotype_matrix, y = ethnicity_pcs, by = c("Sample_ID" = "sample_id"))
phenotype_pcs <- phenotype_matrix %>% 
  filter(!qc_gender_discrepancy & call_rate >= 0.99) %>% 
  filter(Project == "FinnGeDi") %>% 
  filter(Sample_ID %in% colnames(beta_matrix)) %>% 
  group_by(Sample_IID) %>% 
  slice(which.max(call_rate)) %>% 
  ungroup() %>% 
  group_by(gsub("[[:alpha:]]", "", Sample_ID)) %>% 
  filter(n() == 2) %>% 
  ungroup()


### Analysis =======================================================================================
for (.model_type in c("mumwtgain", "primipara", "mumwtgain + primipara")) {
  default_covariates <- c(
    mothers = paste0("AGE_m + BMI_m + ", .model_type),
    offsprings = "log(birthweight_c) + gestationalweek"
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
    (function(x) {setNames(x[["CellT"]], tolower(paste0(x[["status"]], "s")))})()
  
  walk(traits, function(trait) {
    om_status <- "lmm"
    models <- list(
      model = c("simple", "simple_cell")[2],
      formula_o = c(
        as.formula(glue('~ { default_covariates["offsprings"] }')),
        as.formula(glue('~ { default_covariates["offsprings"] } + { default_cells["offsprings"] }'))
      ),
      formula_m = c(
        as.formula(glue('~ { default_covariates["mothers"] }')),
        as.formula(glue('~ { default_covariates["mothers"] } + { default_cells["mothers"] }'))
      )
    )
    
    if (!file.exists(file.path(output_directory, trait, glue("ewas_{trait}_{om_status}.rds")))) {
      dmp <- pmap(models, function(model, formula_o, formula_m) {
        pheno_ids <- phenotype_pcs %>%
          select(
            Sample_ID,
            !!trait,
            unique(c(all.vars(formula_o), all.vars(formula_m))),
            -starts_with("CellT")
          ) %>%
          drop_na() %>%
          group_by(group = gsub("[[:alpha:]]", "", Sample_ID)) %>%
          filter(n() == 2) %>%
          ungroup()
    
        out <- do.call("rbind", mclapply(
          X = split(
            x = as.data.frame(beta_matrix),
            f = (1:nrow(beta_matrix) - 1) %/% ( (nrow(beta_matrix) - 1) %/% n_cores )
          ),
          .form_o = formula_o,
          .form_m = formula_m,
          .pheno = phenotype_pcs %>%
            select(
              Sample_ID, !!trait, status, ID, SEX,
              unique(c(all.vars(formula_o), all.vars(formula_m)))
            ) %>%
            filter(Sample_ID %in% pheno_ids[["Sample_ID"]]),
          mc.preschedule = FALSE,
          mc.cores = n_cores,
          FUN = function(beta_chunk, .form_o, .form_m, .pheno) {
            map_df(.x = rownames(beta_chunk), .f = function(icpg) {
              lm_pheno <- bind_rows(
                .pheno %>%
                  filter(grepl("c", Sample_ID)) %>%
                  mutate(!!icpg := as.matrix(beta_chunk)[icpg, Sample_ID]) %>%
                  do(
                    add_residuals(
                      data = .,
                      model = lm(formula = update.formula(.form_o, paste0(icpg, " ~ .")), data = .)
                    )
                  ),
                .pheno %>%
                  filter(grepl("m", Sample_ID)) %>%
                  mutate(!!icpg := as.matrix(beta_chunk)[icpg, Sample_ID]) %>%
                  do(
                    add_residuals(
                      data = .,
                      model = lm(formula = update.formula(.form_m, paste0(icpg, " ~ .")), data = .)
                    )
                  )
              )
              suppressMessages({
                res <- lmer(formula = resid ~ GDM * status + SEX + (1 | ID), data = lm_pheno)
              })
    
              tidy(res) %>%
                filter(grepl("GDM|status", term)) %>%
                mutate(CpG = icpg, n = nrow(lm_pheno), is_singular = isSingular(res))
            })
          }
        ))
        out %>%
          rename(pvalue = p.value, se = std.error) %>%
          group_by(term) %>%
          mutate(fdr = p.adjust(pvalue, method = "BH")) %>%
          ungroup() %>%
          inner_join(
            y = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Locations %>%
              as.data.frame() %>%
              rownames_to_column("CpG"),
            by = "CpG"
          ) %>%
          rename(cpg_chr = chr, cpg_pos = pos, cpg_strand = strand) %>%
          inner_join(
            y = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC %>%
              as.data.frame() %>%
              rownames_to_column("CpG"),
            by = "CpG"
          ) %>%
          inner_join(
            y = IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other %>%
              as.data.frame() %>%
              select(UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group) %>%
              rownames_to_column("CpG"),
            by = "CpG"
          )
      })
    
      data_models <- tibble(
        dmp = dmp,
        model = models$model,
        formula = list(as.formula("~ GDM * status + SEX + (1 | ID)")),
        formula_o = models$formula_o,
        formula_m = models$formula_m
      )
    
      dir.create(file.path(output_directory, trait), recursive = TRUE, showWarnings = FALSE, mode = "0755")
      write_rds(
        x = data_models,
        file = file.path(output_directory, trait, glue("ewas_{trait}_{om_status}_{gsub(' + ', '_', .model_type, fixed = TRUE)}.rds"))
      )
    } else {
      data_models <- read_rds(file.path(output_directory, trait, glue("ewas_{trait}_{om_status}_{gsub(' + ', '_', .model_type, fixed = TRUE)}.rds")))
    }
    
    data_models_dmp <- data_models %>%
      select(model, formula, formula_o, formula_m, dmp) %>% 
      unnest(dmp) %>% 
      nest(
        dmp = c(
          effect, group, df,
          estimate, se, statistic, pvalue, CpG, 
          n, fdr, cpg_chr, cpg_pos, cpg_strand, Islands_Name, Relation_to_Island, 
          UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group,
          is_singular
        )
      )
    
    
    ### Output =========================================================================================
    ## DMPs Tables -------------------------------------------------------------------------------------
    pwalk(
      .l = list(data_models$model, data_models$dmp),
      .f = function(.model, .dmp) {
        write_csv(
          x = .dmp,
          file = file.path(
            output_directory, trait, 
            glue("{project_name}_EWAS_DMP_{trait}_{om_status}_{.model}_{gsub(' + ', '_', .model_type, fixed = TRUE)}.csv.gz")
          )
        )
      }
    )
    
    ## DMPs Plots --------------------------------------------------------------------------------------
    pwalk(
      .l = list(
        data_models_dmp$model, 
        data_models_dmp$formula,
        data_models_dmp$term, 
        data_models_dmp$dmp,
        data_models_dmp$formula_o,
        data_models_dmp$formula_m
      ),
      .f = function(.model, .form, .term, .dmp, .form_o, .form_m) {
        ## QQ plot
        p_qq <- ggplot(
          data = .dmp %>%
            select(CpG, pvalue) %>%
            arrange(pvalue) %>%
            mutate(
              exppval = (1:n() - 0.5) / n(),
              gc = median(qnorm(pvalue / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1)
            ) %>%
            mutate(
              labels = paste0("lambda[", .term, "]=='", number_format(accuracy = 0.001)(gc), "'"),
              labels = factor(labels, levels = unique(labels))
            ),
          mapping = aes(x = exppval, y = pvalue, colour = labels, shape = labels)
        ) +
          coord_cartesian(xlim = rev(range(c(1, .dmp[["pvalue"]]))), ylim = rev(range(c(1, .dmp[["pvalue"]])))) +
          geom_abline(intercept = 0, slope = 1, colour = "black") +
          geom_point(size = 0.75) +
          scale_x_continuous(trans = "pval", expand = expansion(c(0, 0.1))) +
          scale_y_continuous(trans = "pval", expand = expansion(c(0, 0.1))) +
          scale_colour_viridis_d(labels = parse_format(), begin = 0.5, end = 0.5) +
          scale_shape_discrete(solid = TRUE, labels = parse_format()) +
          labs(x = "Expected P-value", y = "Observed P-value", colour = NULL, shape = NULL) +
          theme(legend.position = c(0, 1), legend.justification = c(-0.05, 1.05)) +
          guides(shape = guide_legend(override.aes = list(size = 3)))
    
        ## Volcano Plot
        p_volcano <- ggplot(
          data = .dmp,
          mapping = aes(x = estimate, y = pvalue, colour = abs(estimate))
        )  +
          annotate(
            geom = "rect",
            xmin = -Inf, xmax = Inf, ymin = 1, ymax = 0.05,
            fill = "firebrick2", alpha = 0.2, colour = NA
          ) +
          geom_vline(xintercept = 0, linetype = 2) +
          geom_point(size = 0.75) +
          scale_colour_viridis_c(trans = "sqrt", limits = c(0, NA)) +
          scale_y_continuous(trans = "pval", expand = expansion(mult = c(0, 0.2))) +
          labs(x = "Estimate (M-value)", y = "P-value") +
          theme(legend.position = "none")
    
        ## Manhattan Plot
        p_manhattan <- ggplot(
          data = .dmp,
          mapping = aes(x = cpg_pos, y = pvalue, colour = gsub("chr", "", cpg_chr))
        ) +
          annotate(
            geom = "rect",
            xmin = -Inf, xmax = Inf, ymin = 1, ymax = 0.05,
            fill = "firebrick2", alpha = 0.2, colour = NA
          ) +
          geom_manhattan(size = 0.75) +
          scale_colour_manual(values = rep(viridis_pal(begin = 0.25, end = 0.75)(2), 12)) +
          theme(legend.position = "none")
    
        p_all <- ((p_qq + p_volcano) / p_manhattan) +
          plot_annotation(
            title = glue("[{project_name}] EWAS {toupper(om_status)}"),
            # subtitle = gsub(
            #   pattern = paste0(" ", .term, " +"), 
            #   replacement = paste0(' <b style="color:#21908C;">', .term, "</b> +"), 
            #   x = "Model: M<sub>residual</sub> = GDM + GDM &times; STATUS + STATUS + SEX + (1 | FID)",
            #   fixed = TRUE
            # ),
            subtitle = case_when(
              grepl("GDM", .term) & !grepl("status", .term) ~ 
                'Model: M<sub>residual</sub> = <b style="color:#21908C;">GDM</b> + GDM &times; STATUS + STATUS + SEX + (1 | FID)',
              grepl("GDM", .term) & grepl("status", .term) ~ 
                'Model: M<sub>residual</sub> = GDM + <b style="color:#21908C;">GDM &times; STATUS</b> + STATUS + SEX + (1 | FID)',
              !grepl("GDM", .term) & grepl("status", .term) ~ 
                'Model: M<sub>residual</sub> = GDM + GDM &times; STATUS + <b style="color:#21908C;">STATUS</b> + SEX + (1 | FID)'
            ),
            caption = glue(.sep = "<br>",
              "With M<sub>residual</sub> the residuals (&epsilon;) from:",
              '&bull; **Offsprings**: M<sub>value</sub> = { glue_collapse(gsub("CellT_", "", all.vars(.form_o)), sep = " + ", last = " + ") } + &epsilon;',
              '&bull; **Mothers**: M<sub>value</sub> = { glue_collapse(gsub("CellT_", "", all.vars(.form_m)), sep = " + ", last = " + ") } + &epsilon;'
            ),
            tag_levels = "A",
            theme = theme(
              plot.subtitle = element_markdown(), 
              plot.caption = element_markdown()
            )
          )
        ggsave(
          filename = file.path(
            output_directory, trait, 
            glue('{project_name}_EWAS_DMP_{trait}_{om_status}_{.model}_{gsub(":|Offspring", "", .term)}_{gsub(" + ", "_", .model_type, fixed = TRUE)}.png')
          ),
          plot = p_all,
          width = 29.7 - 5, height = 21 - 5, units = "cm",
          dpi = 120
        )
      }
    )
  })
}


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
