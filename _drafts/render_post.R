# Render an R Markdown draft into a Jekyll post.
#
# Usage (from RStudio's Knit button, via the `knit:` field in the Rmd YAML, or
# from the command line):
#   Rscript -e 'source("_drafts/render_post.R"); render_post("_drafts/2026-09-07-few-clusters.Rmd")'
#
# What it does beyond rmarkdown::render():
#   * writes the .md into _posts/ (the Rmd file name supplies the date prefix),
#   * writes figures to images/<slug>_files/figure-gfm/ and points the post at
#     /images/... so they resolve on the live site,
#   * strips the render-only YAML keys (everything after the
#     "# --- render config" marker) so Jekyll sees a clean front matter.
# publish = FALSE writes the .md next to the .Rmd in _drafts/ instead of _posts/
# (Jekyll only builds _drafts/ with `jekyll serve --drafts`), which is handy while
# a post is still being written. Flip to TRUE (or drop the argument) to publish.
render_post <- function(input, encoding = "UTF-8", publish = TRUE) {
  input <- normalizePath(input)
  drafts_dir <- dirname(input)
  site_dir <- normalizePath(file.path(drafts_dir, ".."))
  slug <- sub("^\\d{4}-\\d{2}-\\d{2}-", "", tools::file_path_sans_ext(basename(input)))
  fig_dir <- file.path(site_dir, "images", paste0(slug, "_files"), "figure-gfm")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  old_wd <- setwd(drafts_dir); on.exit(setwd(old_wd), add = TRUE)
  out <- rmarkdown::render(
    input, encoding = encoding, output_dir = if (publish) file.path(site_dir, "_posts") else drafts_dir,
    # gfm with $...$ / $$...$$ math (pandoc >= 3 would otherwise emit $`x`$ and
    # ```math fences, which Jekyll/MathJax do not understand); no line wrapping so
    # inline math never straddles a line break.
    output_options = list(variant = "gfm-tex_math_gfm+tex_math_dollars", preserve_yaml = TRUE,
                          pandoc_args = c("--wrap=none")),
    knit_root_dir = drafts_dir, quiet = TRUE,
    envir = new.env(parent = globalenv())
  )

  # knitr's default fig.path is "<basename>_files/figure-gfm/"; move to images/.
  local_figs <- file.path(drafts_dir, paste0(tools::file_path_sans_ext(basename(input)), "_files"), "figure-gfm")
  if (dir.exists(local_figs)) {
    for (f in list.files(local_figs, full.names = TRUE)) file.copy(f, fig_dir, overwrite = TRUE)
    unlink(dirname(local_figs), recursive = TRUE)
  }
  posts_figs <- file.path(if (publish) file.path(site_dir, "_posts") else drafts_dir,
                          paste0(tools::file_path_sans_ext(basename(input)), "_files"))
  if (dir.exists(posts_figs)) {
    for (f in list.files(file.path(posts_figs, "figure-gfm"), full.names = TRUE)) file.copy(f, fig_dir, overwrite = TRUE)
    unlink(posts_figs, recursive = TRUE)
  }

  md <- readLines(out, warn = FALSE)
  # Strip render-only YAML keys.
  fm_end <- which(md == "---")[2]
  marker <- grep("^# --- render config", md[seq_len(fm_end)])
  if (length(marker)) md <- md[-(marker:(fm_end - 1))]
  # Point figure links at /images/<slug>_files/figure-gfm/.
  # rmarkdown may have made the figure paths absolute (when output_dir != input dir),
  # so replace everything up to and including "<basename>_files/figure-gfm/".
  md <- gsub(paste0("[^() ]*", tools::file_path_sans_ext(basename(input)), "_files/figure-gfm/"),
             paste0("/images/", slug, "_files/figure-gfm/"), md)
  writeLines(md, out)
  message("Wrote ", out)
  invisible(out)
}
