#!/usr/bin/env Rscript
# Regenerate Figure 3 (single-pass multi-population LD decay on chromosome 22)
# as VECTOR PDF + high-resolution PNG.
#
# Why regenerate: the copy embedded in the manuscript is a 833 x 614 px raster,
# which is only ~240 dpi at the 3.475 in width it is printed at. Wiley/iMeta
# line art is normally required at 300-600 dpi, so that copy is likely to fail
# production. The underlying data is fully reproducible -- it is the output of
# one PopLDdecay2 -SubPop run -- so a vector version can be produced instead.
#
# The styling deliberately follows the original figure (title, axis labels,
# legend position, population colours) so it can be dropped in as a replacement.
# It is NOT byte-identical to the authors' original: that image was made by a
# script that is not in the repository.
#
#   Rscript plot_lddecay.r [outprefix]
#
# Reads mp.<POP>.stat.gz from this directory (written by
# step1a_runPopLDdecay2Multi.sh) and writes
# Fig3_lddecay.png / .pdf here, so simply run: Rscript plot_lddecay.r
suppressMessages(library(ggplot2))

.f <- commandArgs(trailingOnly = FALSE)
.fi <- grep("^--file=", .f, value = TRUE)
.HERE <- if (length(.fi)) dirname(normalizePath(sub("^--file=", "", .fi[1]))) else getwd()

args <- commandArgs(trailingOnly = TRUE)
W <- if (length(args) >= 1) args[1] else .HERE
OUT <- if (length(args) >= 2) args[2] else file.path(.HERE, "Fig3_lddecay")

POPS <- c(AFR = "#E31A1C", AMR = "#2CA25F", EAS = "#1F78B4",
          EUR = "#6A3D9A", SAS = "#FF7F00")

read_pop <- function(p) {
  f <- file.path(W, sprintf("mp.%s.stat.gz", p))
  if (!file.exists(f)) stop("missing ", f)
  con <- gzfile(f, "rt")
  hdr <- strsplit(sub("^#", "", readLines(con, n = 1)), "\t")[[1]]
  close(con)
  # comment.char="#" skips the header line; header=FALSE stops R from
  # promoting the first data row to column names; check.names=FALSE keeps
  # "Mean_r^2" intact instead of mangling it to "Mean_r.2".
  d <- read.delim(gzfile(f), comment.char = "#", header = FALSE,
                  col.names = hdr, check.names = FALSE)
  stopifnot("Dist" %in% names(d), "Mean_r^2" %in% names(d))
  data.frame(dist_kb = d$Dist / 1000, r2 = d[["Mean_r^2"]], pop = p)
}

df <- do.call(rbind, lapply(names(POPS), read_pop))
df$pop <- factor(df$pop, levels = names(POPS))
cat(sprintf("rows: %s across %d populations\n", format(nrow(df), big.mark = ","), nlevels(df$pop)))

p <- ggplot(df, aes(x = dist_kb, y = r2, colour = pop)) +
  geom_line(linewidth = 0.55) +
  scale_colour_manual(values = POPS, name = NULL) +
  scale_x_continuous(limits = c(0, 300), breaks = seq(0, 300, 50),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Single-pass multi-populations LD Decay (Chr22)",
       x = "Distance(kb)", y = expression(r^2)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.position = "right",
    legend.key.width = unit(1.1, "lines"),
    legend.text = element_text(size = 11),
    plot.margin = margin(8, 8, 6, 8))

# PNG for the manuscript, PDF for editing
ggsave(paste0(OUT, ".png"), p, width = 5.2, height = 3.9, dpi = 600)
ggsave(paste0(OUT, ".pdf"), p, width = 5.2, height = 3.9, device = cairo_pdf)
cat("wrote", paste0(OUT, ".png"), " (600 dpi)\n")
cat("wrote", paste0(OUT, ".pdf"), " (vector)\n")
