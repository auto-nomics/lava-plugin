# LAVA scan plugin node: the official R program of the legacy
# lava_scan_container wrapper with every parameter moved from Rust
# string-building into LAVA_SCAN_* environment variables. Conventions (same
# as the mvmr/ldsc plugin scripts):
#   - optional params render as empty strings; nzchar() is the [ -n ] test
#   - booleans render as "true"/"false" and convert with as.logical()
#   - string arrays render space-joined and rebuild with strsplit()
# The official multiple-locus workflow stays script-side, exactly like the
# legacy wrapper: process.input runs once, the script iterates the
# user-provided loci table, and run.univ.bivar executes per processable
# locus with failures logged and skipped.
#
# Static input port contract (the v0 plugin DSL has no dynamic ports; the
# legacy wrapper varied port 2 and the sumstats ports with spec fields):
#   AUTONOMICS_INPUT0  input.info table (phenotype, cases, controls)
#   AUTONOMICS_INPUT1  loci table
#   AUTONOMICS_INPUT2  sample.overlap table (always connected)
#   AUTONOMICS_INPUT3  sumstats for phenotypes[1]
#   AUTONOMICS_INPUT4  sumstats for phenotypes[2]

options(width = 200)
sink(Sys.getenv("AUTONOMICS_OUTPUT3"), split = TRUE)
on.exit(sink(), add = TRUE)

expected_inputs <- 5L
if (as.integer(Sys.getenv("AUTONOMICS_INPUT_COUNT")) != expected_inputs) {
  stop("LAVA scan input port count does not match its specification")
}

# Space-joined env lists -> R character vectors (empty text -> length 0).
split_list <- function(text)
  if (nzchar(text)) strsplit(text, " ", fixed = TRUE)[[1]] else character(0)

phenos <- split_list(Sys.getenv("LAVA_SCAN_PHENOS"))
target <- split_list(Sys.getenv("LAVA_SCAN_TARGET"))
if (length(target) == 0) target <- NULL
requested_locus_ids <- split_list(Sys.getenv("LAVA_SCAN_LOCUS_IDS"))
if (length(requested_locus_ids) == 0) requested_locus_ids <- NULL
selected_chr <- suppressWarnings(as.integer(Sys.getenv("LAVA_SCAN_CHR")))
if (is.na(selected_chr)) selected_chr <- NULL
univ_threshold <- as.numeric(Sys.getenv("LAVA_SCAN_UNIV_THRESHOLD"))
adap_thresh <- as.numeric(split_list(Sys.getenv("LAVA_SCAN_ADAP_THRESH")))
if (length(adap_thresh) == 0) adap_thresh <- NULL
p_values <- as.logical(Sys.getenv("LAVA_SCAN_P_VALUES"))
cis <- as.logical(Sys.getenv("LAVA_SCAN_CIS"))

# Legacy validate(): phenotype list shape.
if (length(phenos) == 0 || any(!nzchar(phenos))) {
  stop("phenotypes cannot be empty and cannot contain empty IDs")
}
if (anyDuplicated(phenos)) {
  stop("phenotypes cannot contain duplicate IDs")
}
if (length(phenos) != 2L) {
  stop("the lava_scan port contract provides exactly two sumstats ports; phenotypes must have two entries")
}
# Legacy validate(): target shape.
if (!is.null(target) && (length(target) != 1 || !target %in% phenos)) {
  stop("target must contain exactly one configured phenotype")
}
# Legacy validate(): locus subsetting.
if (!is.null(requested_locus_ids)) {
  if (length(requested_locus_ids) == 0 || any(!nzchar(requested_locus_ids))) {
    stop("locus_ids cannot be empty when provided")
  }
  if (anyDuplicated(requested_locus_ids)) {
    stop("locus_ids cannot contain duplicates")
  }
}
if (!is.null(selected_chr) && !(selected_chr >= 1L && selected_chr <= 23L)) {
  stop("chr must be between 1 and 23 (23 is chromosome X)")
}
if (!is.finite(univ_threshold) || univ_threshold <= 0 || univ_threshold > 1) {
  stop("univ_threshold must be finite and in (0, 1]")
}
if (!is.null(adap_thresh) &&
  (any(!is.finite(adap_thresh)) || any(adap_thresh <= 0))) {
  stop("adap_thresh must be finite and greater than zero")
}

input_info <- read.table(
  Sys.getenv("AUTONOMICS_INPUT0"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("phenotype", "cases", "controls") %in% names(input_info))) {
  stop("input.info must contain phenotype, cases, and controls columns")
}
if (anyDuplicated(input_info$phenotype)) {
  stop("input.info phenotype IDs must be unique")
}
if (!all(phenos %in% input_info$phenotype)) {
  stop("one or more configured phenotypes are missing from input.info")
}
input_info <- input_info[match(phenos, input_info$phenotype), , drop = FALSE]
sumstats_paths <- unname(Sys.getenv(sprintf("AUTONOMICS_INPUT%d", seq_along(phenos) + 2L)))
input_info$filename <- sumstats_paths
normalized_input_info <- file.path(
  Sys.getenv("AUTONOMICS_WORKDIR"),
  ".autonomics",
  "lava_scan",
  "input.info.txt"
)
dir.create(dirname(normalized_input_info), recursive = TRUE, showWarnings = FALSE)
write.table(
  input_info,
  normalized_input_info,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

input <- LAVA::process.input(
  input.info.file = normalized_input_info,
  sample.overlap.file = Sys.getenv("AUTONOMICS_INPUT2"),
  ref.prefix = "/panels/lava_ref/lava-ukb-v1.1",
  phenos = phenos
)
loci <- LAVA::read.loci(Sys.getenv("AUTONOMICS_INPUT1"))
if (anyDuplicated(loci$LOC)) {
  stop("loci table contains duplicate LOC IDs")
}

if (!is.null(requested_locus_ids)) {
  matches <- match(requested_locus_ids, as.character(loci$LOC))
  if (any(is.na(matches))) {
    stop("requested loci are absent: ", paste(requested_locus_ids[is.na(matches)], collapse = ", "))
  }
  loci <- loci[matches, , drop = FALSE]
}
if (!is.null(selected_chr)) {
  loci <- loci[loci$CHR == selected_chr, , drop = FALSE]
}
if (nrow(loci) == 0) {
  stop("loci selection is empty")
}

empty_result <- function() {
  data.frame(
    locus = character(),
    chr = integer(),
    start = integer(),
    stop = integer(),
    n.snps = integer(),
    n.pcs = integer(),
    stringsAsFactors = FALSE
  )
}

univ_results <- list()
bivar_results <- list()
raw_results <- list()
failed <- 0L
skipped <- 0L
cat("Starting official LAVA scan for", nrow(loci), "loci\n")
for (row_index in seq_len(nrow(loci))) {
  locus_row <- loci[row_index, , drop = FALSE]
  locus_id <- as.character(locus_row$LOC)
  cat("Processing locus", locus_id, "\n")
  locus <- tryCatch(
    LAVA::process.locus(
      locus_row,
      input,
      phenos = phenos,
      min.K = 2,
      prune.thresh = 99,
      max.prop.K = 0.75,
      drop.failed = TRUE,
      max.block.size = 3000,
      cap.estimates = TRUE
    ),
    error = function(error) {
      failed <<- failed + 1L
      message <- paste0("Error processing locus ", locus_id, ": ",
                        conditionMessage(error))
      cat(message, "\n")
      NULL
    }
  )
  if (is.null(locus)) {
    skipped <<- skipped + 1L
    next
  }

  result <- tryCatch(
    LAVA::run.univ.bivar(
      locus,
      phenos = phenos,
      target = target,
      univ.thresh = univ_threshold,
      adap.thresh = adap_thresh,
      p.values = p_values,
      CIs = cis,
      param.lim = 1.25,
      cap.estimates = TRUE
    ),
    error = function(error) {
      failed <<- failed + 1L
      message <- paste0("Error analyzing locus ", locus_id, ": ",
                        conditionMessage(error))
      cat(message, "\n")
      NULL
    }
  )
  if (is.null(result)) {
    skipped <<- skipped + 1L
    next
  }

  raw_results[[locus_id]] <- result
  locus_info <- data.frame(
    locus = locus$id,
    chr = locus$chr,
    start = locus$start,
    stop = locus$stop,
    n.snps = locus$n.snps,
    n.pcs = locus$K,
    stringsAsFactors = FALSE
  )
  if (!is.null(result$univ)) {
    univ_results[[length(univ_results) + 1L]] <- cbind(locus_info, result$univ)
  }
  if (!is.null(result$bivar)) {
    bivar_results[[length(bivar_results) + 1L]] <- cbind(locus_info, result$bivar)
  }
}

univ <- if (length(univ_results)) do.call(rbind, univ_results) else empty_result()
bivar <- if (length(bivar_results)) do.call(rbind, bivar_results) else empty_result()
cat(
  "Finished official LAVA scan:",
  nrow(loci), "requested;",
  nrow(univ), "univariate rows;",
  nrow(bivar), "bivariate rows;",
  skipped, "skipped;",
  failed, "failed\n"
)
write.table(univ, Sys.getenv("AUTONOMICS_OUTPUT0"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bivar, Sys.getenv("AUTONOMICS_OUTPUT1"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(loci = loci, univ = univ, bivar = bivar, raw_results = raw_results),
  Sys.getenv("AUTONOMICS_OUTPUT2")
)
