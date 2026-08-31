#!/usr/bin/env Rscript
# Generate paper performance figures from /usr/bin/time -v output files.
# Usage: Rscript fig_perf.R <workdir> <outprefix>
# Reads <workdir>/<name>.time files (bench_chr22.sh / bench_subpop.sh output).
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
work <- if (length(args) >= 1) args[1] else "Test//paper_chr22"
outp <- if (length(args) >= 2) args[2] else "figures/fig_perf"

wall_secs <- function(s) {
  vapply(s, function(x) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    if (length(parts) == 3) sum(as.numeric(parts) * c(3600, 60, 1))
    else sum(as.numeric(parts) * c(60, 1))
  }, numeric(1))
}

files <- list.files(work, pattern = "\\.time$", full.names = TRUE)
rows <- list()
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  name <- sub("\\.time$", "", basename(f))
  wall_line <- grep("Elapsed", txt, value = TRUE)
  rss_line <- grep("Maximum resident", txt, value = TRUE)
  if (length(wall_line) == 0 || length(rss_line) == 0) {
    cat("skip (incomplete):", name, "\n")
    next
  }
  wall <- wall_secs(sub(".*: ", "", wall_line) %>% trimws())
  rss_kb <- as.numeric(gsub("[^0-9]", "", sub(".*: ", "", rss_line) %>% trimws()))
  rows[[name]] <- data.frame(name = name, wall_s = wall, rss_mb = rss_kb / 1024)
}
df <- bind_rows(rows)

df$T <- as.numeric(ifelse(grepl("^new.t", df$name), sub("new.t", "", df$name), NA))
df$tool <- ifelse(grepl("^old", df$name), "PopLDdecay 3.45", "PopLDdecay2")

# --- Figure: thread scaling (PopLDdecay2 only) ---
newdf <- df %>% filter(tool == "PopLDdecay2", !is.na(T)) %>% arrange(T)
newdf$speedup <- newdf$wall_s[newdf$T == 1] / newdf$wall_s

p1 <- ggplot(newdf, aes(x = T, y = wall_s)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_x_continuous(breaks = newdf$T) +
  labs(x = "Threads (-T)", y = "Wall time (s)", title = "PopLDdecay2 thread scaling") +
  theme_bw()
ggsave(paste0(outp, "_scaling.pdf"), p1, width = 6, height = 4)
ggsave(paste0(outp, "_scaling.png"), p1, width = 6, height = 4, dpi = 300)

p2 <- ggplot(newdf, aes(x = T, y = speedup)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 2.5) +
  scale_x_continuous(breaks = newdf$T) +
  labs(x = "Threads (-T)", y = "Speedup vs -T1", title = "PopLDdecay2 speedup") +
  theme_bw()
ggsave(paste0(outp, "_speedup.pdf"), p2, width = 6, height = 4)
ggsave(paste0(outp, "_speedup.png"), p2, width = 6, height = 4, dpi = 300)

# --- Figure: memory comparison (3.45 vs PopLDdecay2 -T1) ---
cmp <- df %>% filter((grepl("^old", name) & grepl("^old.345", name)) | name == "new.t1")
if (nrow(cmp) >= 2) {
  cmp$tool <- factor(cmp$tool, levels = c("PopLDdecay 3.45", "PopLDdecay2"))
  p3 <- ggplot(cmp, aes(x = tool, y = rss_mb, fill = tool)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.0f MB", rss_mb)), vjust = -0.4) +
    labs(x = "", y = "Peak RSS (MB)", title = "Peak memory: 3.45 vs PopLDdecay2 (-T1)") +
    scale_fill_manual(values = c("gray60", "steelblue")) +
    theme_bw() + theme(legend.position = "none")
  ggsave(paste0(outp, "_mem.pdf"), p3, width = 5, height = 4)
  ggsave(paste0(outp, "_mem.png"), p3, width = 5, height = 4, dpi = 300)
}

# --- Summary table to stdout ---
cat("== perf summary ==")
print(df)

# --- Table 1 data row (TSV) ---
cat("\n== Table 1 rows (name\twall_s\trss_mb) ==")
write.table(df[c("name", "wall_s", "rss_mb")], row.names = FALSE, sep = "\t", quote = FALSE)
