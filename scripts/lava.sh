# LAVA plugin node: the official R program of the legacy lava_container
# wrapper with every parameter moved from Rust string-building into LAVA_*
# environment variables. Conventions (same as the mvmr/ldsc plugin scripts):
#   - optional params render as empty strings; nzchar() is the [ -n ] test
#   - booleans render as "true"/"false" and convert with as.logical()
#   - string arrays render space-joined and rebuild with strsplit()
# The validation block reproduces the message set of the legacy Rust
# validate(); it runs inside the container instead of at DAG build time.
# The LAVA API tokens (process.input -> read.loci -> process.locus -> run.*)
# are identical to the legacy script.
#
# Static input port contract (the v0 plugin DSL has no dynamic ports; the
# legacy wrapper varied port 2 and the sumstats ports with spec fields):
#   AUTONOMICS_INPUT0  input.info table (phenotype, cases, controls)
#   AUTONOMICS_INPUT1  loci table
#   AUTONOMICS_INPUT2  sample.overlap table (always connected)
#   AUTONOMICS_INPUT3  sumstats for phenotypes[1]
#   AUTONOMICS_INPUT4  sumstats for phenotypes[2]

options(width = 200)
sink(Sys.getenv("AUTONOMICS_OUTPUT2"), split = TRUE)
on.exit(sink(), add = TRUE)

expected_inputs <- 5L
if (as.integer(Sys.getenv("AUTONOMICS_INPUT_COUNT")) != expected_inputs) {
  stop("LAVA input port count does not match its specification")
}

# Space-joined env lists -> R character vectors (empty text -> length 0).
split_list <- function(text)
  if (nzchar(text)) strsplit(text, " ", fixed = TRUE)[[1]] else character(0)

phenos <- split_list(Sys.getenv("LAVA_PHENOS"))
analysis <- Sys.getenv("LAVA_ANALYSIS")
target <- split_list(Sys.getenv("LAVA_TARGET"))
if (length(target) == 0) target <- NULL

# Legacy validate(): analysis membership (a serde enum rejection in Rust).
if (!analysis %in% c("univ", "bivar", "pcor", "multireg")) {
  stop("unknown LAVA analysis; known analyses: univ, bivar, pcor, multireg")
}
# Legacy validate(): phenotype list shape.
if (length(phenos) == 0 || any(!nzchar(phenos))) {
  stop("phenotypes cannot be empty and cannot contain empty IDs")
}
if (anyDuplicated(phenos)) {
  stop("phenotypes cannot contain duplicate IDs")
}
if (length(phenos) != 2L) {
  stop("the lava port contract provides exactly two sumstats ports; phenotypes must have two entries")
}
# Legacy validate(): target subset plus per-analysis shape.
if (!is.null(target) && !all(target %in% phenos)) {
  stop("target must reference one of phenotypes")
}
if (analysis == "univ" && !is.null(target)) {
  stop("target is not used by the univ analysis")
}
if (analysis == "bivar" && (!is.null(target) && length(target) != 1)) {
  stop("bivar target must contain one phenotype or be null")
}
if (analysis == "pcor" && (is.null(target) || length(target) != 2 || target[1] == target[2])) {
  stop("pcor target must contain exactly two distinct phenotypes")
}
if (analysis == "multireg" && (is.null(target) || length(target) != 1)) {
  stop("multireg target must contain one outcome phenotype")
}
# Legacy validate(): adaptive thresholds and max.r2.
adap_thresh <- as.numeric(split_list(Sys.getenv("LAVA_ADAP_THRESH")))
if (length(adap_thresh) == 0) adap_thresh <- NULL
if (!is.null(adap_thresh) &&
  (any(!is.finite(adap_thresh)) || any(adap_thresh <= 0))) {
  stop("adap_thresh must be finite and greater than zero")
}
max_r2 <- as.numeric(Sys.getenv("LAVA_MAX_R2"))
if (!is.finite(max_r2) || max_r2 <= 0) {
  stop("max_r2 must be finite and greater than zero")
}
p_values <- as.logical(Sys.getenv("LAVA_P_VALUES"))
cis <- as.logical(Sys.getenv("LAVA_CIS"))
variances <- as.logical(Sys.getenv("LAVA_VARIANCES"))
only_full_model <- as.logical(Sys.getenv("LAVA_ONLY_FULL_MODEL"))
# Legacy container_spec(): run.multireg relaxes param.lim to 1.5; every
# other analysis uses 1.25.
param_lim <- if (analysis == "multireg") 1.5 else 1.25

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
normalized_input_info <- file.path(Sys.getenv("AUTONOMICS_WORKDIR"), ".autonomics", "lava", "input.info.txt")
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

# Locus selection: an explicit LOC id wins over the one-based row index.
# The legacy container_spec() rendered either match("<id>",
# as.character(loci$LOC)) or the bare index and shared the bounds check.
requested_locus_id <- Sys.getenv("LAVA_LOCUS_ID")
if (nzchar(requested_locus_id)) {
  locus_selector <- match(requested_locus_id, as.character(loci$LOC))
} else {
  locus_selector <- as.integer(Sys.getenv("LAVA_LOCUS_INDEX"))
  if (is.na(locus_selector) || locus_selector < 1L) {
    stop("locus_index must be one-based")
  }
}
if (is.na(locus_selector) || locus_selector > nrow(loci)) {
  stop("selected locus is outside the loci table")
}
locus_row <- loci[locus_selector, , drop = FALSE]
locus <- LAVA::process.locus(
  locus_row,
  input,
  phenos = phenos,
  min.K = 2,
  prune.thresh = 99,
  max.prop.K = 0.75,
  drop.failed = TRUE,
  max.block.size = 3000,
  cap.estimates = TRUE
)
if (is.null(locus)) { stop("official LAVA could not process the selected locus") }

result <- switch(
  analysis,
  univ = LAVA::run.univ(locus, phenos = phenos, var = variances, cap.estimates = TRUE),
  bivar = LAVA::run.bivar(locus, phenos = phenos, target = target, adap.thresh = adap_thresh, p.values = p_values, CIs = cis, param.lim = param_lim, cap.estimates = TRUE),
  pcor = LAVA::run.pcor(locus, phenos = phenos, target = target, adap.thresh = adap_thresh, p.values = p_values, CIs = cis, max.r2 = max_r2, param.lim = param_lim),
  multireg = LAVA::run.multireg(locus, phenos = phenos, target = target, adap.thresh = adap_thresh, only.full.model = only_full_model, p.values = p_values, CIs = cis, param.lim = param_lim, suppress.message = FALSE)
)
if (is.null(result)) { stop("official LAVA returned no result") }

flatten_result <- function(value, model = character()) {
  if (is.data.frame(value)) {
    output <- value
    if (length(model) > 0) output$model <- paste(model, collapse = "/")
    return(output)
  }
  if (is.null(names(value)) || identical(names(value), character(0))) {
    child_results <- lapply(value, flatten_result, model = model)
  } else {
    child_results <- Map(function(name, child) flatten_result(child, c(model, name)), names(value), value)
  }
  results <- Filter(Negate(is.null), child_results)
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}
flat_result <- flatten_result(result)
print(locus_row)
print(result)
write.table(
  flat_result,
  Sys.getenv("AUTONOMICS_OUTPUT0"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(list(locus = locus_row, analysis = analysis, result = result), Sys.getenv("AUTONOMICS_OUTPUT1"))
