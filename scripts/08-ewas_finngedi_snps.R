# Docker Container R403
### Load packages ==================================================================================
library(here)

suppressPackageStartupMessages({
  library(Homo.sapiens)
  # Use prefix instead of attaching these big bad written packages
  # library(limma)
  # library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  # library(DMRcate)
  # library(DMRcatedata)
  # library(ggbio)

  library(readr)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(glue)
  library(forcats)
  library(patchwork)
  library(gt)
  library(utils)
  library(data.table)
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


### Functions ======================================================================================
draw_dmr <- function(.dmr, .dmp, base_colour = viridis_pal(begin = 0.5, end = 0.5)(1)) {
  pval_trans_inf <- function() {
    neglog10_breaks <- function(n = 5) {
      function(x) {
        rng <- -log(range(x, na.rm = TRUE), base = 10)
        min <- 0
        max <- floor(rng[1])
        if (max == min) {
          10^-min
        } else {
          by <- floor((max - min) / n) + 1
          10^-seq(min, max, by = by)
        }
      }
    }
    trans_new(
      name = "pval",
      transform = function(x) {
        x[x == Inf] <- 0
        x[x == -Inf] <- 1
        -log(x, 10)
      },
      inverse = function(x) 10^-x,
      breaks = neglog10_breaks(),
      format = function(x) {
        parse(
          text = gsub("e", " %*% 10^", gsub("1e+00", "1", scientific_format()(x), fixed = TRUE))
        )
      },
      domain = c(0, 1)
    )
  }

  dmr_range <- as.data.frame(.dmr@ranges)
  .dmp <- .dmp[.dmp[["cpg_pos"]] >= dmr_range[["start"]] & .dmp[["cpg_pos"]] <= dmr_range[["end"]], ]

  if (!is.na(.dmr$overlapping.genes)) {
    p_gene <- ggbio::autoplot(
      object = Homo.sapiens::Homo.sapiens,
      which = .dmr,
      colour = base_colour,
      fill = base_colour
    )@ggplot +
      scale_y_continuous(expand = expansion(0.15), position = "right") +
      facet_grid(rows = vars("Gene"), switch = "y") +
      theme(plot.margin = unit(x = c(5.5, 5.5, 0, 5.5), units = "pt"))
  }

  p_estimate <- ggplot(
    data = .dmp,
    mapping = aes(x = cpg_pos, y = estimate, ymin = estimate - se, ymax = estimate + se)
  ) +
    geom_hline(yintercept = 0, linetype = 2, size = 0.2)  +
    geom_pointrange(size = 0.2, show.legend = FALSE, colour = base_colour) +
    geom_line(colour = base_colour) +
    scale_y_continuous(position = "right") +
    facet_grid(rows = vars("Estimate (M-value)"), switch = "y") +
    if (!is.na(.dmr$overlapping.genes)) {
      theme(plot.margin = unit(x = c(0, 5.5, 0, 5.5), units = "pt"))
    } else {
      theme(plot.margin = unit(x = c(5.5, 5.5, 0, 5.5), units = "pt"))
    }

  p_pval <- ggplot(
    data = distinct(.dmp, CpG, cpg_chr, cpg_pos, pvalue), #cpg_strand
    mapping = aes(x = cpg_pos, y = pvalue)
  ) +
    geom_point(colour = base_colour) +
    annotate(
      geom = "rect",
      xmin = -Inf, xmax = Inf, ymin = 1, ymax = 0.05,
      fill = "firebrick2", alpha = 0.2, colour = NA
    ) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "firebrick2") +
    scale_y_continuous(
      position = "right",
      trans = pval_trans_inf(),
      expand = expansion(mult = c(0.002, 0.05))
    ) +
    facet_grid(rows = vars("P-value"), switch = "y") +
    theme(plot.margin = unit(x = c(0, 5.5, 0, 5.5), units = "pt"))

  p_density <- ggplot(
    data = distinct(.dmp, CpG, cpg_chr, cpg_pos), #cpg_strand,
    mapping = aes(x = cpg_pos)
  ) +
    geom_density(stat = "density", colour = base_colour, fill = base_colour, alpha = 0.6) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      position = "right",
      labels = function(x) {
        parse(
          text = gsub("e", " %*% 10^", gsub("1e+00", "1", x, fixed = TRUE))
        )
      }
    ) +
    facet_grid(rows = vars("Density of CpG"), switch = "y") +
    theme(plot.margin = unit(x = c(0, 5.5, 5.5, 5.5), units = "pt"))

  if (!is.na(.dmr$overlapping.genes)) {
    lp <- list(p_gene, p_estimate, p_pval, p_density)
  } else {
    lp <- list(p_estimate, p_pval, p_density)
  }

  lp <- map2(
    .x = lp,
    .y = as.list(c(rep(FALSE, length(lp) - 1), TRUE)),
    .range = range(unlist(map(lp, ~ layer_scales(.x)$x$range$range))),
    .f = function(p, x_axis, .range) {
      out <- p +
        scale_x_continuous(
          labels = function(l) paste0(prettyNum(l / 1e6, digits = 5), " Mb"),
          limits = c(.range)
        ) +
        theme_light() +
        theme(
          axis.title = element_blank(),
          legend.title = element_blank(),
          strip.background = element_rect(fill = NA),
          strip.text = element_text(colour = "grey30")
        )
        if (!x_axis) {
          out <- out +
            theme(
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank()
            )
        }
      out
    }
  )

  lp <- map(lp, `+`,
    geom_rect(
      data = dmr_range,
      mapping = aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "dodgerblue",
      colour = NA,
      alpha = 0.1
    )
  )
  wrap_plots(lp, ncol = 1)
}


### Environment ====================================================================================
data_directory <- normalizePath("/disks/DATA/Projects/EpxGDM/FINNGEDI_TUBINGEN/QC")
output_directory <- here("outputs", "08-ewas_finngedi_snps")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0755")

project_name <- gsub("(.*)_.*", "\\1", basename(here()))

theme_set(theme_light())


### Traits to be analysed ==========================================================================
traits <- c("Gestational Diabetes Mellitus" = "GDM")


### Data ===========================================================================================
beta_matrix <- fread(
  file = file.path(data_directory, "EPIC", "EPIC_QC_betavalues_SNPs.csv.gz")
) %>%
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

phenotype_matrix <- read_csv(
  file = file.path(data_directory, "EPIC", "EPIC_QC_phenotypes_SNPs.csv")
)

# phenotype_pcs <- left_join(x = phenotype_matrix, y = ethnicity_pcs, by = c("Sample_ID" = "sample_id"))
phenotype_pcs <- phenotype_matrix %>% 
  filter(!qc_gender_discrepancy & call_rate >= 0.99) %>% 
  filter(Project == "FinnGeDi") %>% 
  filter(Sample_ID %in% colnames(beta_matrix)) %>% 
  group_by(Sample_IID) %>% 
  slice(which.max(call_rate)) %>% 
  ungroup()
beta_matrix <- beta_matrix[, as.character(phenotype_pcs[["Sample_ID"]])]


### Analyses =======================================================================================
default_covariates <- c(
  mothers = "AGE_m + BMI_m",
  offsprings = "SEX + log(birthweight_c) + gestationalweek"
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
  walk(c("mothers", "offsprings"), function(om_status) {
    models <- tribble(
      ~model, ~formula,
      "simple", as.formula(glue("~ { trait } + { default_covariates[om_status] }")),
      "simple_cell", as.formula(glue("~ { trait } + { default_covariates[om_status] } + { default_cells[om_status] }"))
    )
  
    data_models <- models %>%
      mutate(
        pheno = pmap(
          .l = list(formula, model),
          .f = function(.form, .model) {
            phenotype_pcs %>%
              filter(tolower(paste0(status, "s")) == !!om_status) %>% 
              select(Sample_ID, all.vars(.form)) %>%
              drop_na()
          }
        ),
        beta = map(pheno, ~ beta_matrix[, as.character(.x[["Sample_ID"]])])
      ) %>%
      mutate(
        dmp = pmap(
          .l = list(formula, beta, pheno),
          .f = function(.form, .beta, .pheno) {
            limma_fit1 <- limma::lmFit(.beta, model.matrix(object = .form, data = .pheno))
            limma_fit2 <- limma::eBayes(limma_fit1)
            limma::topTable(
              fit = limma_fit2,
              coef = trait,
              number = nrow(.beta),
              adjust.method = "BH",
              p.value = 1,
              sort.by = "none"
            ) %>%
              rownames_to_column("CpG") %>%
              rename(
                estimate = logFC,
                avgmvalue_meth = AveExpr,
                t_statistic = t,
                pvalue = P.Value
              ) %>%
              mutate(
                se = sqrt(limma_fit2$s2.post) * limma_fit2$stdev.unscaled[, trait],
                fdr = p.adjust(pvalue, method = "BH"),
                n = nrow(.pheno)
              ) %>%
              select(-adj.P.Val, -B) %>%
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
          }
        )
      )
  
    data_models <- data_models %>%
      mutate(
        dmr = pmap(
          .l = list(formula, beta, pheno),
          .f = function(.form, .beta, .pheno) {
            suppressPackageStartupMessages({
              .dmp_for_dmr <- DMRcate::cpg.annotate(
                datatype = "array",
                object = .beta,
                what = "M",
                arraytype = "EPIC",
                analysis.type = "differential",
                design = model.matrix(object = .form, data = .pheno),
                coef = 2,
                fdr = 1
              )
            })
  
            DMRcate::extractRanges(
              DMRcate::dmrcate(
                object = .dmp_for_dmr,
                lambda = 1000, # default
                C = NULL, # default
                min.cpgs = 2 # default
              ),
              genome = "hg19"
            )
          }
        )
      )
  
    dir.create(file.path(output_directory, trait), recursive = TRUE, showWarnings = FALSE, mode = "0755")
    write_rds(
      x = select(data_models, -pheno, -beta), 
      file = file.path(output_directory, trait, glue("ewas_snps_{trait}_{om_status}.rds"))
    )
  
  
    ### Output =====================================================================================
    ## DMPs Tables ---------------------------------------------------------------------------------
    pwalk(
      .l = list(data_models$model, data_models$dmp),
      .f = function(.model, .dmp) {
        write_csv(
          x = .dmp,
          file = file.path(
            output_directory, trait, 
            glue("{project_name}_EWAS_DMP_SNPs_{trait}_{om_status}_{.model}.csv.gz")
          )
        )
      }
    )
  
    ## DMPs Plots ----------------------------------------------------------------------------------
    pwalk(
      .l = list(data_models$model, data_models$formula, data_models$dmp),
      .f = function(.model, .form, .dmp) {
        ## QQ plot
        p_qq <- ggplot(
          data = .dmp %>%
            select(CpG, pvalue) %>%
            pivot_longer(cols = "pvalue", names_to = "type", values_to = "value") %>%
            group_by(type) %>%
            arrange(value) %>%
            mutate(
              exppval = (1:n() - 0.5) / n(),
              gc = median(qnorm(value / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1)
            ) %>%
            ungroup() %>%
            arrange(type) %>%
            mutate(
              labels = paste0("lambda[", type, "]=='", number_format(accuracy = 0.001)(gc), "'"),
              labels = factor(labels, levels = unique(labels))
            ),
          mapping = aes(x = exppval, y = value, colour = labels, shape = labels)
        ) +
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
            title = glue("[{project_name}] EWAS SNPs {toupper(om_status)}"),
            subtitle = glue('Model: M-value = { glue_collapse(gsub("CellT_", "", all.vars(.form)), sep = " + ", last = " + ") }'),
            tag_levels = "A"
          )
  
        ggsave(
          filename = file.path(
            output_directory, trait, 
            glue("{project_name}_EWAS_DMP_SNPs_{trait}_{om_status}_{.model}.png")
          ),
          plot = p_all,
          width = 29.7 - 5, height = 21 - 5, units = "cm",
          dpi = 120
        )
      }
    )
  
    ## DMPs QQ Plot --------------------------------------------------------------------------------
    p_qq <- ggplot(
      data = map(data_models$dmp, ~ select(.x, CpG, pvalue)) %>%
        reduce(.f = full_join, by = "CpG") %>%
        set_names(nm = c("CpG", data_models$model)) %>%
        pivot_longer(cols = starts_with("simple"), names_to = "type", values_to = "value") %>%
        mutate(
          type = factor(
            x = type,
            levels = c("simple", "simple_cell"),
            labels = c("RAW", "Cells")
          ),
          type_group = fct_collapse(
            .f = type,
            "Without Cellular Composition" = "RAW",
            "With Cellular Composition" = "Cells"
          ),
          type_lambda = fct_collapse(
            .f = type,
            "RAW" = c("RAW", "Cells")
          )
        ) %>%
        group_by(type) %>%
        arrange(value) %>%
        mutate(
          exppval = (1:n() - 0.5) / n(),
          gc = median(qnorm(value / 2)^2, na.rm = TRUE) / qchisq(0.5, df = 1)
        ) %>%
        ungroup() %>%
        arrange(type) %>%
        mutate(
          labels = paste0("lambda['", type, "']=='", number_format(accuracy = 0.001)(gc), "'"),
          labels = factor(labels, levels = unique(labels))
        ),
      mapping = aes(x = exppval, y = value, colour = labels, shape = labels)
    ) +
      geom_abline(intercept = 0, slope = 1, colour = "black") +
      geom_point(size = 0.75) +
      scale_x_continuous(trans = "pval", expand = expansion(c(0, 0.1))) +
      scale_y_continuous(trans = "pval", expand = expansion(c(0, 0.1))) +
      scale_colour_viridis_d(labels = parse_format(), begin = 0, end = 0.95) +
      scale_shape_discrete(solid = TRUE, labels = parse_format()) +
      labs(x = "Expected P-value", y = "Observed P-value", colour = NULL, shape = NULL) +
      theme(legend.position = "top") +
      guides(shape = guide_legend(override.aes = list(size = 3))) +
      facet_grid(rows = vars(type_lambda), cols = vars(type_group))
  
    ggsave(
      filename = file.path(
        output_directory, trait, 
        glue("{project_name}_EWAS_DMP_SNPs_{trait}_{om_status}_ALL.png")
      ),
      plot = p_qq,
      width = 16, height = 12, units = "cm",
      dpi = 120
    )
  
    ## DMRs Tables ---------------------------------------------------------------------------------
    pwalk(
      .l = list(data_models$model, data_models$dmr),
      .f = function(.model, .dmr) {
        write_csv(
          x = as.data.frame(.dmr),
          file = file.path(
            output_directory, trait, 
            glue("{project_name}_EWAS_DMR_SNPs_{trait}_{om_status}_{.model}.csv.gz")
          )
        )
      }
    )
  
    ## DMRs Plots ----------------------------------------------------------------------------------
    # draw_dmr(.dmr = data_models$dmr[[1]][2, ], .dmp = data_models$dmp[[1]])
  
    # summary(data.frame(data_models$dmr[[1]]))
  })
})


### Summary ========================================================================================
walk(traits, function(trait) {
  data_models <- map_df(c("mothers", "offsprings"), function(om_status) {
    data_models <- read_rds(file.path(output_directory, trait, glue("ewas_snps_{trait}_{om_status}.rds")))
    mutate(data_models, status = om_status)
  }) %>%
    mutate(
      trait = trait,
      dmp_dmr = map2(
        .x = dmp, .y = dmr,
        .f = function(x, y) {
          full_join(
            x = map_df(c(0.10, 0.05, 0.01), function(alpha) {
              as.data.table(x)[, 
                .(
                  alpha = alpha, 
                  pvalue = sum(pvalue <= alpha), 
                  fdr = sum(fdr <= alpha), 
                  n_dmp = .N
                )
              ]
            }),
            y = map_df(c(0.10, 0.05, 0.01), function(alpha) {
              as.data.table(y)[, 
                .(
                  alpha = alpha, 
                  mean_n_cpgs = mean(no.cpgs),
                  sd_n_cpgs = sd(no.cpgs),
                  min_smoothed_fdr = sum(min_smoothed_fdr <= alpha),
                  HMFDR = sum(HMFDR <= alpha),
                  Fisher = sum(Fisher <= alpha),
                  n_dmr = .N
                )
              ]
            }),
            by = "alpha"
          ) %>%
            arrange(desc(alpha)) %>%
            mutate(alpha = paste("α =", alpha))
        }
      ),
      dmp = NULL,
      dmr = NULL,
      model = NULL
    ) %>%
    unnest("dmp_dmr") %>%
    mutate(
      formula = map(formula, function(x) {
        glue('Model: M-value = { glue_collapse(gsub("CellT_", "", all.vars(x)), sep = " + ", last = " + ") }')
      })
    ) %>%
    gt(rowname_col = "alpha", groupname_col = c("trait", "status", "formula")) %>%
    fmt_number(
      columns = vars(pvalue, fdr, n_dmp, mean_n_cpgs, sd_n_cpgs, min_smoothed_fdr, HMFDR, Fisher, n_dmr),
      decimals = 0
    ) %>%
    tab_spanner(
      label = "DMP",
      columns = vars(pvalue, fdr, n_dmp)
    ) %>%
    tab_spanner(
      label = "DMR",
      columns = vars(mean_n_cpgs, min_smoothed_fdr, HMFDR, Fisher, n_dmr)
    ) %>%
    cols_merge_uncert(
      col_val = vars(mean_n_cpgs, sd_n_cpgs),
      col_uncert = vars(sd_n_cpgs)
    ) %>%
    fmt_missing(columns = everything()) %>%
    cols_label(
      pvalue = "P-value",
      fdr = "FDR",
      n_dmp = "N",
      mean_n_cpgs = "# CpGs",
      min_smoothed_fdr = "min(FDR)",
      HMFDR = "HM-FDR",
      Fisher = "Fisher",
      n_dmr = "N"
    ) %>%
    tab_header(
      title = "Summary Results",
      subtitle = md(glue("**{ project_name }** Project"))
    ) %>%
    tab_footnote("FDR: false discovery rate.", cells_column_labels("fdr")) %>%
    tab_footnote("Harmonic mean of the individual CpG FDRs.", cells_column_labels("HMFDR")) %>%
    tab_footnote("Fisher combined probability transform of the individual CpG FDRs.", cells_column_labels("Fisher")) %>%
    tab_footnote("Minimum FDR of the smoothed estimate.", cells_column_labels("min_smoothed_fdr")) %>%
    tab_footnote("Number of CpGs in the DMR as mean (± standard deviation).", cells_column_labels("mean_n_cpgs")) %>%
    tab_footnote("DMP: differentially methylated position.", cells_column_spanners("DMP")) %>%
    tab_footnote("DMR: differentially methylated region.", cells_column_spanners("DMR")) %>%
    opt_row_striping() %>%
    tab_style(
      style = list(
        cell_fill(color = "firebrick2"),
        cell_text(color = "white")
      ),
      locations = cells_body(columns = vars(fdr), rows = fdr == 0)
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = "firebrick2"),
        cell_text(color = "white")
      ),
      locations = cells_body(columns = vars(HMFDR), rows = HMFDR == 0)
    ) %>%
    tab_style(
      style = list(
        cell_fill(color = "firebrick2"),
        cell_text(color = "white")
      ),
      locations = cells_body(columns = vars(Fisher), rows = Fisher == 0)
    ) %>%
    gtsave(file.path(output_directory, trait, glue("{ project_name }_EWAS_SNPs_summary.html")))
})


### Archive ========================================================================================
owd <- getwd()
walk(traits, function(trait) {
  setwd(file.path(output_directory, trait))
  archive_name <- file.path(
    output_directory, 
    glue('{gsub("-", "", Sys.Date())}_{project_name}_EWAS_SNPs_{trait}.zip')
  )
  zip(
    zipfile = archive_name,
    files = list.files(pattern = "gz$|png$|html$", full.names = TRUE)
  )
  setwd(owd)
  system(glue("fexsend {archive_name} ."))
  unlink(archive_name)
})
