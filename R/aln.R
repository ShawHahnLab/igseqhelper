# make a data frame of factors, one for each list combination from the given
# annotations, naming the regions within the alignment based off of the selected
# sequences.
prep_seq_annotations <- function(alignment, annots) {
  annots_grid <- list()
  for (seq_name in names(annots)) {
    if (seq_name %in% names(alignment)) {
      labelled_seq <- alignment[[seq_name]]
      gaps <- cumsum(strsplit(labelled_seq, "")[[1]] == "-")
      for (item in names(annots[[seq_name]])) {
        annots_here <- annots[[seq_name]][[item]]
        # Use gaps to shift positions listed in the annotations, and then set
        # up the factors.  careful; need to shift iteratively I think.
        for (idx in 1:length(annots_here)) {
          shift_1 <- gaps[annots_here[[idx]][1]]
          shift_2 <- gaps[annots_here[[idx]][2]]
          annots_here[[idx]][1] <- annots_here[[idx]][1] + shift_1
          annots_here[[idx]][2] <- annots_here[[idx]][2] + shift_2
          if (idx < length(annots_here)) {
            for (idx2 in (idx + 1):length(annots_here)) {
              annots_here[[idx2]][1] <- annots_here[[idx2]][1] + shift_2
              annots_here[[idx2]][2] <- annots_here[[idx2]][2] + shift_2
            }
          }
        }
        # Now, set up factors with these shifted positions.
        # Note that the first matching position wins, if there are overlapping
        # ranges.
        matches <- lapply(1:length(annots_here), function(annot_idx) {
          annots_here[[annot_idx]][1] <= seq(1, nchar(labelled_seq)) &
            annots_here[[annot_idx]][2] >= seq(1, nchar(labelled_seq))
        })
        idx <- apply(do.call(cbind, matches), 1, function(vec) match(TRUE, vec))
        vec <- factor(names(annots_here)[idx], levels = names(annots_here))
        annots_grid[[paste(seq_name, item)]] <- vec
      }
    }
  }
  annots_grid <- do.call(cbind.data.frame, annots_grid)
  annots_grid
}

#' Make matrix version of alignment
#'
#' Take a sequence alignment in character vector form and return a matrix
#' equivalent.  If the string lengths aren't identical, shorter sequences will
#' be padded with NA in the matrix.
#'
#' @param aln_vec character vector of aligned sequences
#' @return matrix of individual characters
make_aln_mat <- function(aln_vec) {
  if (length(unique(nchar(aln_vec))) > 1) {
    warning("length mismatch")
  }
  if (! is.null(names(aln_vec)) &&
      length(unique(names(aln_vec))) != length(aln_vec)) {
    warning("duplicate names")
  }
  do.call(rbind, lapply(strsplit(aln_vec, ""), `[`, 1:max(nchar(aln_vec))))
}

#' Make data frame version of alignment
#'
#' Take character vector of already-aligned sequences and create a long-format
#' data frame with a row for each sequence and position, noting characteristics
#' relative to a reference sequence within the alignment (by default the first).
#'
#' @param aln_vec named character vector of already-aligned sequences
#' @param index of reference within aln_vec
#' @param squeeze vector of one or more characters (.e.g, "-").  If any column
#'   contains only a single one of these characters that column will be removed.
#' @importFrom magrittr "%>%"
make_aln_df <- function(aln_vec, ref = 1, squeeze = NULL) {
  # Convert vector of strings into matrix and corresponding matrix of logicals
  # for match/mismatch to alignment
  aln_mat <- make_aln_mat(aln_vec)
  aln_mat_mask <- t(apply(aln_mat, 1, `==`, aln_mat[ref, ]))
  if (! is.null(squeeze)) {
    idxl_keep <- ! apply(aln_mat, 2, function(vec) {
      length(unique(vec)) == 1 && unique(vec) %in% squeeze
    })
    aln_mat <- aln_mat[, idxl_keep]
    aln_mat_mask <- aln_mat_mask[, idxl_keep]
  }
  # Create long-format data frames for each of these
  # Handle names carefully to avoid any assumptions
  aln_mat_names <- rownames(aln_mat)
  rownames(aln_mat) <- NULL
  rownames(aln_mat_mask) <- NULL
  reshape2::melt(
    aln_mat,
    varnames = c("SeqNum", "Position"),
    value.name = "Base") %>%
    as.data.frame() %>%
    dplyr::mutate(SeqNumPos = paste(SeqNum, Position)) ->
    grid_al
  reshape2::melt(
    aln_mat_mask,
    varnames = c("SeqNum", "Position"),
    value.name = "MatchesReference") %>%
    as.data.frame() %>%
    dplyr::mutate(SeqNumPos = paste(SeqNum, Position)) %>%
    dplyr::select(-c(SeqNum, Position)) ->
    grid_mask
  # Combine into one full data frame
  aln_df <- dplyr::full_join(grid_al, grid_mask, by = "SeqNumPos")
  aln_df$SeqNumPos <- NULL
  if (is.null(aln_mat_names)) {
    aln_df$Name <- NA
  } else {
    aln_df$Name <- aln_mat_names[aln_df$SeqNum]
  }
  aln_df$Name <- factor(aln_df$Name, levels = unique(aln_df$Name))
  aln_df[aln_df$SeqNum == ref, "MatchesReference"] <- FALSE
  aln_df$BaseMasked <- as.character(aln_df$Base)
  aln_df$BaseMasked[aln_df$MatchesReference] <- " "

  aln_df
}

# Take a general alignment data frame and style it up to be plotted.
prep_aln_df <- function(aln_df, style) {
  aln_df_for_plot <- aln_df

  fill_colors <- list(
    vsref = c(
      "A" = "#88ff88",
      "T" = "#ff8888",
      "C" = "#8888ff",
      "G" = "#ffff88",
      "-" = "#CCCCCC",
      " " = "#FFFFFF"
    ),
    dnaplotr = c(
      "A" = "#00ff00",
      "T" = "#ff0000",
      "C" = "#0000ff",
      "G" = "#ffff00",
      "-" = "#BEBEBE"
    )
  )

  if (style == "vsref") {
    aln_df_for_plot$plot_fill <-
      fill_colors[[style]][as.character(aln_df_for_plot$BaseMasked)]
    aln_df_for_plot$plot_text <- aln_df_for_plot[["Base"]]
    aln_df_for_plot$plot_text_alpha <- ifelse(
      aln_df_for_plot$MatchesReference, 0.5, 1)
  } else if (style == "dnaplotr") {
    aln_df_for_plot$plot_fill <- fill_colors[[style]][
      as.character(aln_df_for_plot$Base)]
    aln_df_for_plot$plot_text <- ""
    aln_df_for_plot$plot_text_alpha <- 1
  } else {
    stop("style \"", style, "\" not recognized")
  }

  # Set up a vertical axis for sequence labels
  if (! any(is.na(aln_df_for_plot$Name)) &&
      sum(diag(table(aln_df_for_plot$SeqNum, aln_df_for_plot$Name))) ==
      nrow(aln_df_for_plot)) {
    # every SeqNum has a unique name
    aln_df_for_plot$plot_row_label <- aln_df_for_plot$Name
  } else {
    # or, we'll make up a consistent label for each.  Important that it's a
    # factor for any potential use of discrete scale re-labelling to work as
    # expected.
    txt <- as.character(aln_df_for_plot$Name)
    txt <- ifelse(
      is.na(txt),
      aln_df_for_plot$SeqNum,
      paste(aln_df_for_plot$SeqNum, txt))
    aln_df_for_plot$plot_row_label <- factor(txt)
  }

  aln_df_for_plot
}

#' Plot Aligned Sequences
#'
#' @param aln_vec character vector of already-aligned sequences
#' @param ref index of reference within aln_vec
#' @param aln_df insted of aln_vec, use this already-prepared data frame
#' @param draw_text write text labels on plot?
#' @param style keyword for plotting style to use.  One of: dnaplotr, vsref
#'
#' @return ggplot2 plot object
plot_aln <- function(aln_vec, ref = 1, aln_df=NULL, draw_text=TRUE, style="vsref") {
  if (is.null(aln_df)) {
    aln_df <- make_aln_df(aln_vec, ref)
  }
  aln_df_for_plot <- prep_aln_df(aln_df, style)

  p <- ggplot2::ggplot(
    aln_df_for_plot,
    ggplot2::aes(x = Position, y = plot_row_label, fill = plot_fill)) +
    ggplot2::geom_tile()
  if (! any(is.na(aln_df_for_plot$Name))) {
    ylabels <- unique(subset(aln_df_for_plot, select = c(plot_row_label, Name)))
    p <- p + ggplot2::scale_y_discrete(
      labels = ylabels$Name,
      breaks = ylabels$plot_row_label)
  }
  if (draw_text){
    p <- p + ggplot2::geom_text(
      ggplot2::aes(
        x = Position,
        y = plot_row_label,
        label = plot_text,
        alpha = plot_text_alpha),
      size = 1.25,
      show.legend = FALSE)
  }

  fill_vec <- c(vsref="BaseMasked", dnaplotr="Base")[[style]]
  clrs <- unique(aln_df_for_plot[, c(fill_vec, "plot_fill")])
  clrs[, 1] <- factor(clrs[, 1])
  clrs <- clrs[order(clrs[, 1]), ]
  p <- p +
    ggplot2::scale_fill_identity(
      guide = "legend",
      labels = clrs[, 1],
      breaks = clrs[, 2],
      na.value = "white",
      name = "") +
    ggplot2::coord_cartesian(
      xlim = c(-1, max(aln_df_for_plot$Position) + 3), expand = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))

  if ("Facet" %in% colnames(aln_df_for_plot)) {
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(Facet),
      scales = "free",
      space = "free")
  }

  p
}
