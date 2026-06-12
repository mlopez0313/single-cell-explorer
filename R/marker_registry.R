# ============================================================================
# Marker registry
# ----------------------------------------------------------------------------
# Typed, versioned, state-level resource. Lives on `state$marker_registry`
# and is queried by annotation engines (`marker_score`, future SingleR /
# CellTypist) and the annotation UI. No module hardcodes marker definitions.
#
# Schema (one entry per cell type):
#   cell_type            character(1)   display label, e.g. "T cell"
#   ontology_id          character(1)   Cell Ontology id, e.g. "CL:0000084"
#                                       (NA is allowed; the slot must exist
#                                       because retrofitting it later means
#                                       touching every entry.)
#   parent_cell_type     character(1)   hierarchy parent (NA for roots)
#   species              character(1)   "human" / "mouse" / NA
#   tissue               character()    one or more tissues, or "various"
#   markers              list(marker_gene)  list of per-gene entries:
#       gene             character(1)
#       role             "positive" | "negative" | "specific"
#       weight           numeric(1)    contribution to the cell-type score
#       evidence         "protein_validated" | "scRNA_only" | "inferred"
#       aliases          named character: species -> gene symbol (Cd3d/CD3D)
#   source               character(1)  "builtin" | "PanglaoDB" | "CellMarker"
#                                      | "Azimuth" | "user"
#   source_version       character(1)
#   confidence_threshold numeric(1)    source-recommended cutoff, NA if none
#   notes                character(1)
#
# The full registry object additionally carries:
#   version              character(1)  used by `marker_registry_version`
#                                      stamping on annotation results
#   entries              list(marker_entry)
#   created_at           POSIXct
#
# TODO (future ingest paths, not in this file):
#   - PanglaoDB              tsv ingest -> marker_entries
#   - CellMarker 2.0         tsv ingest
#   - Azimuth reference      programmatic extract
#   - User-defined overlay   highest precedence layered on top
# ============================================================================

MARKER_REGISTRY_SCHEMA_VERSION <- "marker_registry_v1"

#' Build a single marker-gene record.
marker_gene <- function(gene, role = c("positive", "negative", "specific"),
                        weight = 1.0,
                        evidence = c("protein_validated", "scRNA_only", "inferred"),
                        aliases = list()) {
  role     <- match.arg(role)
  evidence <- match.arg(evidence)
  list(
    gene     = as.character(gene),
    role     = role,
    weight   = as.numeric(weight),
    evidence = evidence,
    aliases  = as.list(aliases)
  )
}

#' Build a single cell-type entry for the registry.
#'
#' `ontology_id` defaults to `NA_character_` -- the slot is mandatory so that
#' future ontology-aware features (cross-dataset merging, hierarchical
#' coloring) do not require touching every entry.
marker_entry <- function(cell_type,
                         ontology_id          = NA_character_,
                         parent_cell_type     = NA_character_,
                         species              = NA_character_,
                         tissue               = character(),
                         markers              = list(),
                         source               = "builtin",
                         source_version       = NA_character_,
                         confidence_threshold = NA_real_,
                         notes                = "") {
  stopifnot(is.character(cell_type), length(cell_type) == 1L, nzchar(cell_type))
  stopifnot(is.list(markers))
  list(
    cell_type            = cell_type,
    ontology_id          = as.character(ontology_id),
    parent_cell_type     = as.character(parent_cell_type),
    species              = as.character(species),
    tissue               = as.character(tissue),
    markers              = markers,
    source               = as.character(source),
    source_version       = as.character(source_version),
    confidence_threshold = as.numeric(confidence_threshold),
    notes                = as.character(notes)
  )
}

#' Built-in starter registry. Covers the cell types referenced by the mock
#' dataset's example genes (CD3D, MS4A1, LST1, EPCAM, COL1A1, NKG7).
#'
#' Intentionally small. Real deployments should ingest from PanglaoDB /
#' CellMarker / Azimuth and merge into this structure.
default_marker_registry <- function() {
  entries <- list(
    marker_entry(
      cell_type = "T cell", ontology_id = "CL:0000084",
      parent_cell_type = "Lymphocyte",
      species = "human", tissue = c("PBMC", "lymph_node"),
      markers = list(
        marker_gene("CD3D",  "positive", 1.0, "protein_validated",
                    aliases = list(human = "CD3D", mouse = "Cd3d")),
        marker_gene("CD3E",  "positive", 1.0, "protein_validated"),
        marker_gene("MS4A1", "negative", 0.8, "protein_validated"),
        marker_gene("NKG7",  "negative", 0.4, "scRNA_only")
      ),
      source = "builtin"
    ),
    marker_entry(
      cell_type = "NK cell", ontology_id = "CL:0000623",
      parent_cell_type = "Lymphocyte",
      species = "human", tissue = "PBMC",
      markers = list(
        marker_gene("NKG7", "positive", 1.0, "protein_validated"),
        marker_gene("CD3D", "negative", 0.6, "protein_validated"),
        marker_gene("MS4A1", "negative", 0.5, "protein_validated")
      ),
      source = "builtin"
    ),
    marker_entry(
      cell_type = "B cell", ontology_id = "CL:0000236",
      parent_cell_type = "Lymphocyte",
      species = "human", tissue = c("PBMC", "lymph_node"),
      markers = list(
        marker_gene("MS4A1", "positive", 1.0, "protein_validated"),
        marker_gene("CD3D",  "negative", 0.8, "protein_validated")
      ),
      source = "builtin"
    ),
    marker_entry(
      cell_type = "Myeloid cell", ontology_id = "CL:0000763",
      parent_cell_type = NA_character_,
      species = "human", tissue = "PBMC",
      markers = list(
        marker_gene("LST1", "positive", 1.0, "scRNA_only"),
        marker_gene("CD3D", "negative", 0.5, "protein_validated")
      ),
      source = "builtin"
    ),
    marker_entry(
      cell_type = "Epithelial cell", ontology_id = "CL:0000066",
      parent_cell_type = NA_character_,
      species = "human", tissue = c("various"),
      markers = list(
        marker_gene("EPCAM", "positive", 1.0, "protein_validated")
      ),
      source = "builtin"
    ),
    marker_entry(
      cell_type = "Fibroblast", ontology_id = "CL:0000057",
      parent_cell_type = NA_character_,
      species = "human", tissue = c("stroma", "various"),
      markers = list(
        marker_gene("COL1A1", "positive", 1.0, "protein_validated")
      ),
      source = "builtin"
    )
  )
  list(
    schema_version = MARKER_REGISTRY_SCHEMA_VERSION,
    version        = "builtin_v0.1.0",
    source         = "builtin",
    created_at     = Sys.time(),
    entries        = entries
  )
}

# ---- Query helpers --------------------------------------------------------

#' All cell-type labels in the registry (display strings).
marker_registry_cell_types <- function(registry) {
  if (is.null(registry) || is.null(registry$entries)) return(character())
  vapply(registry$entries, `[[`, character(1), "cell_type")
}

#' Filter the registry. NULL arguments are "no filter on this axis".
#' Returns the matching list of marker_entry records (not a wrapped registry).
marker_registry_filter <- function(registry,
                                   species = NULL, tissue = NULL,
                                   cell_type = NULL, gene = NULL) {
  if (is.null(registry) || is.null(registry$entries)) return(list())
  out <- registry$entries
  if (!is.null(species)) {
    out <- Filter(function(e) is.na(e$species) || identical(e$species, species), out)
  }
  if (!is.null(tissue)) {
    out <- Filter(function(e) {
      length(intersect(e$tissue, tissue)) > 0L || "various" %in% e$tissue
    }, out)
  }
  if (!is.null(cell_type)) {
    out <- Filter(function(e) identical(e$cell_type, cell_type), out)
  }
  if (!is.null(gene)) {
    out <- Filter(function(e) {
      any(vapply(e$markers, function(m) identical(m$gene, gene), logical(1)))
    }, out)
  }
  out
}

#' Lookup a single entry by cell_type. Returns NULL if not found.
marker_registry_get <- function(registry, cell_type) {
  hits <- marker_registry_filter(registry, cell_type = cell_type)
  if (length(hits) == 0L) NULL else hits[[1]]
}

#' Union of marker gene symbols in the registry (or a filtered subset).
marker_registry_genes <- function(registry, cell_type = NULL) {
  entries <- if (is.null(cell_type)) (registry$entries %||% list())
             else marker_registry_filter(registry, cell_type = cell_type)
  if (!length(entries)) return(character())
  unique(unlist(lapply(entries, function(e) {
    vapply(e$markers, `[[`, character(1), "gene")
  })))
}

#' Ontology-id lookup for a label (or NA if not in registry).
marker_registry_ontology_id <- function(registry, cell_type) {
  e <- marker_registry_get(registry, cell_type)
  if (is.null(e)) NA_character_ else (e$ontology_id %||% NA_character_)
}

#' Build the ontology map for a vector of labels: label -> ontology_id.
#' Labels not in the registry get NA_character_. Useful for stamping onto
#' annotation_result_v1$ontology_map.
build_ontology_map <- function(registry, labels) {
  if (!length(labels)) return(setNames(character(), character()))
  vapply(unique(stats::na.omit(labels)), function(lab) {
    marker_registry_ontology_id(registry, lab)
  }, FUN.VALUE = character(1)) -> ids
  names(ids) <- unique(stats::na.omit(labels))
  ids
}
