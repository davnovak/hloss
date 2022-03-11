
.JaccardSimilarityHeatmap <- function(
  jt, counts, eval_bijective, eval_fixedc1, eval_fixedc2, unassigned, title, c1_name = 'c1', c2_name = 'c2'
) {

  ## Change Jaccard-table column order to get a nice diagonal from the 1:1 matching
  ## Move rows corresponding to unassigned's to the bottom

  matches_bijective   <- as.vector(stats::na.exclude(eval_bijective$Matches))
  unassigneds_present <- if (is.null(unassigned)) c() else intersect(unassigned, rownames(jt))

  order_c1 <- c(names(eval_bijective$Matches), unassigneds_present)
  order_c2 <- c(matches_bijective, names(eval_fixedc2$Matches)[!names(eval_fixedc2$Matches) %in% matches_bijective])

  row_ordering <- match(order_c1, rownames(jt))
  col_ordering <- match(order_c2, colnames(jt))

  if (length(col_ordering) > length(matches_bijective)) {
    ## Change Jaccard-table column order to get another potential diagonal right of the submatrix corresponding to the 1:1 matches
    ## This is achieved by ranking columns by the fixed-c2 matching

    matches_fixedc2       <- eval_fixedc2$Matches
    additional_c2_matches <- names(matches_fixedc2)[!names(matches_fixedc2) %in% matches_bijective]
    additional_c2_matches <- additional_c2_matches[!additional_c2_matches %in% names(matches_fixedc2)[is.na(matches_fixedc2)]]
    unmatched_c2          <- colnames(jt)[!colnames(jt) %in% matches_bijective & !colnames(jt) %in% additional_c2_matches]

    matched_c1 <- matches_fixedc2[match(additional_c2_matches, names(matches_fixedc2))]
    ranking <- match(matched_c1, rownames(jt)[row_ordering])
    col_ordering[(length(matches_bijective) + 1):length(col_ordering)] <- c(names(matched_c1)[order(ranking)], unmatched_c2)
  }

  jt     <- jt[row_ordering, col_ordering]
  counts <- counts[row_ordering, col_ordering]

  ## Retrieve Precision, Recall and F1 scores for each match for each matching approach
  ## Format them as sidebar annotations (data frames accepted by pheatmap)

  annot_row <- data.frame(
    Precision = eval_fixedc1$Precision.PerMatch,
    Recall = eval_fixedc1$Recall.PerMatch,
    F1 = eval_fixedc1$F1.PerMatch,
    'Precision_Bijective' = eval_bijective$Precision.PerMatch,
    'Recall_Bijective' = eval_bijective$Recall.PerMatch,
    'F1_Bijective' = eval_bijective$F1.PerMatch
  )
  rownames(annot_row) <- names(eval_fixedc1$Precision.PerMatch)

  annot_col <- data.frame(
    Precision = eval_fixedc2$Precision.PerMatch,
    Recall = eval_fixedc2$Recall.PerMatch,
    F1 = eval_fixedc2$F1.PerMatch
  )
  rownames(annot_col) <- names(eval_fixedc2$Precision.PerMatch)

  ## Produce the basic heatmap

  ph <- pheatmap::pheatmap(
    jt, main = title, legend = FALSE, display_numbers = counts, fontsize_number = 6, cluster_cols = FALSE, cluster_rows = FALSE,
    annotation_row = annot_row, annotation_col = annot_col, annotation_legend = FALSE, border_color = 'grey',
    gaps_row = min(nrow(jt) - length(unassigneds_present), sum(!is.na(matches_bijective))),
    gaps_col = length(matches_bijective), number_color = 'black', silent = TRUE,
    annotation_colors =
      list(
        Precision           = RColorBrewer::brewer.pal(6, 'Greens'),
        Recall              = RColorBrewer::brewer.pal(6, 'Blues'),
        F1                  = RColorBrewer::brewer.pal(6, 'Purples'),
        Precision_Bijective = RColorBrewer::brewer.pal(6, 'Greens'),
        Recall_Bijective    = RColorBrewer::brewer.pal(6, 'Blues'),
        F1_Bijective        = RColorBrewer::brewer.pal(6, 'Purples')
      )
  )

  ## Access and tweak graphical parameters of the grob corresponding to cells of the heatmap

  grob_classes <- lapply(ph$gtable$grobs, class)
  idx_maingrob <- which(sapply(grob_classes, function(cl) 'gTree' %in% cl))[1]
  grob_names   <- names(ph$gtable$grobs[[idx_maingrob]]$children)
  rects_name   <- grob_names[grep('rect', grob_names)]
  gp           <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp # graphical parameters of the rectangles

  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$col <- gp$fill
  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$lwd <-
    matrix(0, nrow = nrow(gp$fill), ncol = ncol(gp$fill), dimnames = list(rownames(gp$fill), colnames(gp$fill)))

  ## Highlight cells corresponding to columns chosen by fixed-c1 matching by drawing pink frames around them
  ## To that end, create a new 'overlay' grob with the pink frames to the grob tree

  overlay_grob         <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]
  overlay_grob$gp$fill <- NA
  m                    <- data.frame(Group.c1 = names(eval_fixedc1$Matches), Group.c2 = eval_fixedc1$Matches)
  m                    <- m[!is.na(m$Group.c2),]
  for (idx_match in seq_len(nrow(m))) {
    overlay_grob$gp$col[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- '#ff96e1'
    overlay_grob$gp$lwd[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- 2
  }
  ph$gtable$grobs[[idx_maingrob]]$children$overlay_fixedc1 <- overlay_grob

  ## Highlight cells corresponding to rows chosen by fixed-c2 matching by drawing light-green frames around them
  ## To that end, create a new 'overlay' grob with the light-green frames to the grob tree (smaller than the pink ones, to avoid overlaps)

  overlay_grob         <- ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]
  overlay_grob$gp$fill <- NA
  overlay_grob$height  <- overlay_grob$height - grid::unit(0.6, 'strheight', as.character(overlay_grob$height))
  overlay_grob$width   <- overlay_grob$width - grid::unit(0.05, 'strwidth', as.character(overlay_grob$width))
  m                    <- data.frame(Group.c2 = names(eval_fixedc2$Matches), Group.c1 = eval_fixedc2$Matches)
  m                    <- m[!is.na(m$Group.c1),]
  for (idx_match in seq_len(nrow(m))) {
    overlay_grob$gp$col[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- '#85ffa5'
    overlay_grob$gp$lwd[m$Group.c1[idx_match], m$Group.c2[idx_match]] <- 2
  }
  ph$gtable$grobs[[idx_maingrob]]$children$overlay_fixedc2 <- overlay_grob

  text_grob_name    <- ph$gtable$grobs[[idx_maingrob]]$childrenOrder['text']
  other_grobs_names <- ph$gtable$grobs[[idx_maingrob]]$childrenOrder[ph$gtable$grobs[[idx_maingrob]]$childrenOrder != text_grob_name]
  ph$gtable$grobs[[idx_maingrob]]$childrenOrder <- c(other_grobs_names, 'overlay_fixedc1', 'overlay_fixedc2', text_grob_name)
  ph$gtable$grobs[[idx_maingrob]]$children[[rects_name]]$gp$lwd <- 1

  ## Make the labels for row and column annotation (with stats for each match) prettier

  idx_grob <- which(sapply(ph$gtable$grobs, function(x) !is.null(x$label) && 'Precision_Bijective' %in% x$label))
  ph$gtable$grobs[[idx_grob]]$label   <- c(paste0('Precision (Fixed ', c1_name, ')'), paste0('Recall (Fixed ', c1_name, ')'), paste0('F1 (Fixed ', c1_name, ')'), 'Precision (bijective)', 'Recall (bijective)', 'F1 (bijective)')
  ph$gtable$grobs[[idx_grob]]$gp$font <- c(2L, 2L, 2L, 3L, 3L, 3L)

  idx_grob <- which(sapply(ph$gtable$grobs, function(x) !is.null(x$label) && 'Precision' %in% x$label))
  ph$gtable$grobs[[idx_grob]]$label   <- c(paste('Precision (Fixed ', c2_name, ')'), paste0('Recall (Fixed ', c2_name, ')'), paste0('F1 (Fixed ', c2_name, ')'))

  ## The labels of those groups in c1 which belong to unassigned will be in italics

  if (length(unassigneds_present) > 0) {
    idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(rownames(jt) %in% x$label)))
    font_vec <- rep(1L, length(ph$gtable$grobs[[idx_grob]]$label))
    font_vec[length(font_vec):(length(font_vec - length(unassigneds_present)) + 1)] <- 3L
    ph$gtable$grobs[[idx_grob]]$gp$font <- font_vec
  }

  ## Put c2 groups used in bijective matching in bold and rotate the labels of c2 groups

  idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(colnames(jt) %in% x$label)))
  if (length(idx_grob) > 1)
    idx_grob <- which(sapply(ph$gtable$grobs, function(x) all(rownames(jt) %in% x$label)))
  font_vec <- rep(1L, length(ph$gtable$grobs[[idx_grob]]$label))
  font_vec[1:length(eval_fixedc1$Matches)] <- 2L

  ph$gtable$grobs[[idx_grob]]$gp$font <- font_vec
  ph$gtable$grobs[[idx_grob]]$rot     <- 0
  ph$gtable$grobs[[idx_grob]]$hjust   <- 0.5
  ph$gtable$grobs[[idx_grob]]$vjust   <- 1.5

  cowplot::plot_grid(
    NULL,
    ph[[4]],
    grid::textGrob(
      paste0('One-to-one (bijective) mapping corresponds to matrix diagonal.\nMapping of groups from ', c2_name, ' to groups from ', c1_name, ' corresponds to green frames.\nMapping of groups from ', c1_name, ' to groups from ', c2_name, ' corresponds to pink frames.'),
      gp = grid::gpar(fontsize = 12)
    ),
    nrow = 3, rel_heights = c(0.2, 15, 2)
  )
}

