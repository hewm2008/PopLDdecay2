#!/usr/bin/env Rscript
# Figure 2: (A) real time, (B) peak RSS, (C) strong-scaling efficiency.
#
# Every number comes from measurements; nothing is typed in below. Two sources:
#
#   (default)        the .time files that step1_runThreads.sh wrote in this
#                    directory; the aggregates are mirrored to perf_data.tsv
#   PERF_DATA=<tsv>  read an aggregate table instead of the .time files
#                    (same layout as the perf_data.tsv this script writes);
#                    nothing is written back
#
# Conventions, matching the paper caption:
#   * error bars are +/- 1 SD of the measured quantity (sample SD, divisor n-1)
#   * the speed-up label is a ratio of MEANS: mean(PopLDdecay)/mean(PopLDdecay2)
#   * no error bar on the speed-up: it is a ratio of means, and we do not define
#     a spread for a derived ratio
#   * efficiency E = T1_wall / (T * wall_at_T)
#
# Panel A keeps the broken y axis: the baseline is ~18x the single-threaded bar,
# so a linear axis would flatten all the others.
#
#   Rscript plot_performance.r [out.png]

suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

.f <- commandArgs(trailingOnly = FALSE)
.fi <- grep("^--file=", .f, value = TRUE)
.HERE <- if (length(.fi)) dirname(normalizePath(sub("^--file=", "", .fi[1]))) else getwd()

THREADS <- c(1, 4, 8, 16, 32)
REPEATS <- 1:5

# Two data sources: the per-run .time records written by step1_runThreads.sh,
# or, with PERF_DATA set, an aggregate table such as the manuscript's.
PERF <- Sys.getenv("PERF_DATA", "")

if (nzchar(PERF)) {
  # aggregate table: "#" comments, key/value pairs, then a 7-column table
  kv <- list(); tab <- list()
  for (ln in readLines(PERF, warn = FALSE)) {
    if (!nzchar(ln) || startsWith(ln, "#")) next
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(f) == 2) kv[[f[1]]] <- suppressWarnings(as.numeric(f[2]))
    else if (f[1] != "threads") tab[[length(tab) + 1]] <- as.numeric(f)
  }
  if (length(tab) != length(THREADS))
    stop("expected ", length(THREADS), " rows in ", PERF, ", found ", length(tab))
  col <- function(i) vapply(tab, `[`, numeric(1), i)
  OLD_WALL <- kv[["old_wall"]]; OLD_WALL_SD <- kv[["old_wall_sd"]]
  OLD_RSS  <- kv[["old_rss"]];  OLD_RSS_SD  <- kv[["old_rss_sd"]]
  NEW_WALL <- col(2); NEW_WALL_SD <- col(3)
  NEW_RSS  <- col(4); NEW_RSS_SD  <- col(5)
  SPEEDUP  <- col(6); EFFIC       <- col(7)
  cat(sprintf("baseline : %.2f +/- %.2f s, %.0f MB (n=%d, from %s)\n",
              OLD_WALL, OLD_WALL_SD, OLD_RSS, kv[["old_n"]], PERF))
} else {
  # "    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:12.38"  -> seconds
  # "    Maximum resident set size (kbytes): 135168"            -> MB
  parse_time <- function(path) {
    x <- readLines(path, warn = FALSE)
    el <- sub(".*: ", "", grep("Elapsed", x, value = TRUE)[1])
    p <- as.numeric(strsplit(el, ":")[[1]])
    wall <- if (length(p) == 3) p[1] * 3600 + p[2] * 60 + p[3] else p[1] * 60 + p[2]
    rss <- as.numeric(sub(".*: ", "", grep("Maximum resident", x, value = TRUE)[1])) / 1024
    c(wall, rss)
  }

  # mean and sample SD over the repeats; the 5th value is the repeat count
  summ <- function(paths) {
    if (!length(paths) || !all(file.exists(paths))) return(rep(NA_real_, 5))
    m <- vapply(paths, parse_time, numeric(2))
    c(mean(m[1, ]), sd(m[1, ]), mean(m[2, ]), sd(m[2, ]), length(paths))
  }

  new <- t(vapply(THREADS,
                  function(T) summ(sprintf("%s/t%d.r%d.time", .HERE, T, REPEATS)),
                  numeric(5)))
  old <- summ(sprintf("%s/old.r%d.time", .HERE, REPEATS))
  if (anyNA(new[, 1]) || is.na(old[1]))
    stop("missing .time files in ", .HERE, " -- run: sh step1_runThreads.sh")

  OLD_WALL <- old[1]; OLD_WALL_SD <- old[2]
  OLD_RSS  <- old[3]; OLD_RSS_SD  <- old[4]
  NEW_WALL <- new[, 1]; NEW_WALL_SD <- new[, 2]
  NEW_RSS  <- new[, 3]; NEW_RSS_SD  <- new[, 4]
  SPEEDUP  <- OLD_WALL / NEW_WALL
  EFFIC    <- 100 * NEW_WALL[1] / (THREADS * NEW_WALL)

  cat(sprintf("baseline : %.2f +/- %.2f s, %.0f MB (n=%d)\n",
              OLD_WALL, OLD_WALL_SD, OLD_RSS, old[5]))

  # mirror the numbers so the figure can be checked against the table
  KV <- c(sprintf("old_wall\t%.2f", OLD_WALL),
          sprintf("old_wall_sd\t%.2f", OLD_WALL_SD),
          sprintf("old_rss\t%.0f", OLD_RSS),
          sprintf("old_rss_sd\t%.2f", OLD_RSS_SD),
          sprintf("old_n\t%d", old[5]),
          "threads\twall\twall_sd\trss\trss_sd\tspeedup\tefficiency")
  TAB <- sprintf("%d\t%.2f\t%.2f\t%.0f\t%.2f\t%.1f\t%.1f",
                 THREADS, NEW_WALL, NEW_WALL_SD, NEW_RSS, NEW_RSS_SD, SPEEDUP, EFFIC)
  writeLines(c("# written by plot_performance.r from the *.time files -- do not edit",
               KV, TAB), file.path(.HERE, "perf_data.tsv"))
}

# ---------------------------------------------------------------- style
BLUE <- "#87CEEB"; GREEN <- "#7CB95D"
GREY <- "grey30"; ACCENT <- "#2E8B57"; ACCENT2 <- "#3FA9F5"
DODGE_W <- 0.72
# Round the axis break so the tick labels stay tidy.
GAP_LO <- ceiling(max(NEW_WALL) / 50) * 50 + 50
GAP_HI <- OLD_WALL * 0.90
GAP_F  <- 0.10
squish <- function(y) ifelse(y >= GAP_HI, GAP_LO + (y - GAP_HI) * GAP_F, y)
theme_p <- theme_bw(base_size = 11) +
  theme(panel.border = element_rect(colour = "grey25"),
        axis.text.x = element_text(colour = "grey25"),
        axis.text.y = element_text(colour = "grey25"))

# ---------------------------------------------------------------- panel A
a_df <- rbind(
  data.frame(x = factor(THREADS, levels = THREADS), tool = "PopLDdecay2",
             wall = NEW_WALL, sd = NEW_WALL_SD),
  data.frame(x = factor(1, levels = THREADS), tool = "PopLDdecay",
             wall = OLD_WALL, sd = OLD_WALL_SD))
a_df$y <- squish(a_df$wall)
# An SD of a few seconds shrinks to well under a pixel once the axis is
# squashed, so it is drawn but not allowed to dominate the plot.
a_df$sd_y <- pmin(squish(a_df$wall + a_df$sd) - squish(a_df$wall),
                  squish(a_df$wall) - squish(a_df$wall - a_df$sd))
# dodge width is DODGE_W, so the PopLDdecay2 bar sits half a dodge to the right
# of the category centre; the labels have to follow it or they land on the gap.
a_lbl <- data.frame(x = factor(THREADS, levels = THREADS),
                    y = squish(NEW_WALL) + 14,
                    txt = sprintf("%.1fx", SPEEDUP))
a_old <- data.frame(x = factor(1, levels = THREADS),
                    y = squish(OLD_WALL) + 14,
                    txt = sprintf("%.0f s", OLD_WALL))

pA <- ggplot() +
  geom_col(data = a_df, aes(x = x, y = y, fill = tool), width = DODGE_W,
           position = position_dodge(width = DODGE_W)) +
  geom_errorbar(data = a_df, aes(x = x, y = y, ymin = y - sd_y, ymax = y + sd_y),
                width = 0.12, colour = GREY, linewidth = 0.3, position = position_dodge(width = DODGE_W)) +
  geom_text(data = a_lbl, aes(x = as.numeric(x) + DODGE_W / 2, y = y, label = txt),
            size = 3.4, colour = ACCENT, fontface = "bold", vjust = 0) +
  geom_text(data = a_old, aes(x = as.numeric(x) - DODGE_W / 2, y = y, label = txt),
            size = 3.4, colour = GREY, fontface = "bold", vjust = 0) +
  annotate("text", x = 0.62, y = GAP_LO, label = "//", size = 4.0,
           colour = "firebrick", fontface = "bold") +
  scale_fill_manual(values = c("PopLDdecay" = BLUE, "PopLDdecay2" = GREEN),
                    name = NULL) +
  scale_y_continuous(breaks = c(seq(0, GAP_LO, by = 50)),
                     limits = c(0, GAP_LO + 60), expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0.95)) +
  labs(x = "Threads", y = "real time (s)") +
  theme_p + theme(legend.position = "top",
                  legend.text = element_text(size = 8, colour = "grey25"),
                  legend.key.size = unit(0.7, "lines"),
                  plot.margin = margin(6, 4, 4, 4))

# ---------------------------------------------------------------- panel B
b_df <- rbind(
  data.frame(x = factor(THREADS, levels = THREADS), tool = "PopLDdecay2",
             rss = NEW_RSS, sd = NEW_RSS_SD),
  data.frame(x = factor(1, levels = THREADS), tool = "PopLDdecay",
             rss = OLD_RSS, sd = OLD_RSS_SD))
save_pct <- 100 * (OLD_RSS - NEW_RSS) / OLD_RSS
b_lbl <- data.frame(x = factor(THREADS, levels = THREADS),
                    y = NEW_RSS + 90, txt = sprintf("%.1f%%", save_pct))

pB <- ggplot() +
  geom_col(data = b_df, aes(x = x, y = rss, fill = tool), width = DODGE_W,
           position = position_dodge(width = DODGE_W)) +
  geom_errorbar(data = b_df, aes(x = x, ymin = rss - sd, ymax = rss + sd),
                width = 0.12, colour = GREY, linewidth = 0.3,
                position = position_dodge(width = DODGE_W)) +
  geom_text(data = b_lbl, aes(x = as.numeric(x) + DODGE_W / 2, y = y, label = txt),
            size = 3.4, colour = ACCENT2, fontface = "bold") +
  scale_fill_manual(values = c("PopLDdecay" = BLUE, "PopLDdecay2" = GREEN),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  scale_x_discrete(expand = expansion(add = 0.95)) +
  labs(x = "Threads", y = "Peak RSS (Mb)") +
  theme_p

# ---------------------------------------------------------------- panel C
c_df <- data.frame(x = factor(THREADS, levels = THREADS), eff = EFFIC)
c_lbl <- data.frame(x = factor(THREADS, levels = THREADS),
                    y = EFFIC + 7, txt = sprintf("%.1f%%", EFFIC))

pC <- ggplot() +
  geom_col(data = c_df, aes(x = x, y = eff, fill = "PopLDdecay2"),
           width = DODGE_W) +
  geom_hline(yintercept = 100, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_text(data = c_lbl, aes(x = x, y = y, label = txt), size = 3.4,
            colour = GREY, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("PopLDdecay2" = GREEN), guide = "none") +
  scale_x_discrete(expand = expansion(add = 1.6)) +
  scale_y_continuous(breaks = seq(0, 100, by = 25),
                     limits = c(0, 118), expand = expansion(mult = c(0, 0))) +
  annotate("text", x = 5.35, y = 100, label = "ideal 100%", size = 2.8,
           colour = "grey45", vjust = -0.5) +
  labs(x = "Threads", y = "Strong-scaling efficiency (%)") +
  theme_p

# ---------------------------------------------- layout: A|B on top, C below
top <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1.15, 1))
fig <- plot_grid(top, pC, ncol = 1, rel_heights = c(1.55, 1))

# Write both a high-resolution PNG and a vector PDF. The PDF is the one to edit:
# every panel, axis, label and error bar stays vector.
ARGS <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(ARGS) >= 1) ARGS[1] else file.path(.HERE, "Fig2_performance.png")
OUT <- sub("\\.(png|pdf)$", "", OUT, ignore.case = TRUE)
ggsave(paste0(OUT, ".png"), fig, width = 11, height = 7.4, dpi = 300)
ggsave(paste0(OUT, ".pdf"), fig, width = 11, height = 7.4, device = cairo_pdf)
cat("wrote", paste0(OUT, ".png"), "and", paste0(OUT, ".pdf"), "\n")

cat(sprintf("speedup  : %s\n", paste(sprintf("T%d=%.1fx", THREADS, SPEEDUP), collapse = "  ")))
cat(sprintf("effic.   : %s\n", paste(sprintf("T%d=%.1f%%", THREADS, EFFIC), collapse = "  ")))
cat(sprintf("RSS saved: %s\n", paste(sprintf("T%d=%.1f%%", THREADS, save_pct), collapse = "  ")))
