# ============================================================
#  GSE245301 — PUBLIC DATA DOWNLOAD
# ============================================================
# Downloads the publicly available per-sample RNA-seq count files for
# GSE245301 from GEO and extracts them into GSE245301/raw/, where
# GSE245301_pipeline.R expects them.
#
# Source: GSE245301 (mouse liver, ND / HFD-4w / HFD-12w)
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245301
#
# Running this script is optional: GSE245301_pipeline.R downloads the
# same archive automatically if GSE245301/raw/ is missing or empty.
# ============================================================

options(timeout = 900)

tar_path <- "GSE245301_RAW.tar"
raw_dir  <- file.path("GSE245301", "raw")

url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/",
  "GSE245301/suppl/GSE245301_RAW.tar"
)

if (!file.exists(tar_path)) {
  cat("Downloading", url, "\n")
  ok <- tryCatch(
    {
      utils::download.file(url, tar_path, mode = "wb", quiet = TRUE)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!isTRUE(ok)) {
    system(paste("curl -4 -fL -C - -o", shQuote(tar_path), shQuote(url)))
  }
} else {
  cat("Archive already present:", tar_path, "\n")
}

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
utils::untar(tar_path, exdir = raw_dir)

cat("Raw count files extracted to", raw_dir, ":\n")
print(list.files(raw_dir))
