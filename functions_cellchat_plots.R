# CellChat plotting helpers


library(ggplot2)
library(dplyr)
library(circlize)
library(ComplexHeatmap)

theme_academic <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),

      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.3),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.x = element_blank(),

      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_blank(),
      legend.key = element_rect(color = "black", size = 0.3),
      legend.key.size = unit(4, "mm"),
      legend.key.width = unit(12, "mm"),
      legend.text = element_text(size = base_size + 1),
      legend.title = element_blank(),
      legend.margin = margin(5, 8, 5, 5),
      legend.box.spacing = unit(4, "mm"),

      plot.title = element_text(size = base_size + 1, face = "bold",
                                 hjust = 0.5, margin = margin(b = 8)),

      plot.margin = margin(10, 10, 5, 5)
    )
}

## per-cell-type outgoing/incoming counts and strengths from a CellChat object
extract_signaling_data <- function(cellchat_obj, group_name) {
  cell_types <- levels(cellchat_obj@idents)

  outgoing_counts <- rowSums(cellchat_obj@net$count)
  outgoing_strength <- rowSums(cellchat_obj@net$weight)

  incoming_counts <- colSums(cellchat_obj@net$count)
  incoming_strength <- colSums(cellchat_obj@net$weight)

  data.frame(
    cell_type = cell_types,
    outgoing_count = outgoing_counts[cell_types],
    incoming_count = incoming_counts[cell_types],
    outgoing_strength = outgoing_strength[cell_types],
    incoming_strength = incoming_strength[cell_types],
    group = group_name,
    stringsAsFactors = FALSE
  )
}

## grouped bar plot of outgoing/incoming count or strength, AD vs NC (Fig 5b/d)
create_cellchat_academic_plot <- function(combined_data, metric = "count", direction = "outgoing",
                                           group_colors = c("grey40", "#ac0000"),
                                           group_labels = c("NC", "AD"),
                                           plot_title = "",
                                           y_expand = 1.2) {

  if (metric == "count" && direction == "outgoing") {
    y_var <- "outgoing_count"
    y_label <- ""
    if (plot_title == "") plot_title <- "Number of Outgoing Signals"
  } else if (metric == "count" && direction == "incoming") {
    y_var <- "incoming_count"
    y_label <- ""
    if (plot_title == "") plot_title <- "Number of Incoming Signals"
  } else if (metric == "strength" && direction == "outgoing") {
    y_var <- "outgoing_strength"
    y_label <- ""
    if (plot_title == "") plot_title <- "Outgoing Signal Strength"
  } else if (metric == "strength" && direction == "incoming") {
    y_var <- "incoming_strength"
    y_label <- ""
    if (plot_title == "") plot_title <- "Incoming Signal Strength"
  }

  combined_data$group <- factor(combined_data$group, levels = c("NC", "AD"))

  max_value <- max(combined_data[[y_var]], na.rm = TRUE)
  y_limit <- c(0, max_value * y_expand)

  ggplot(combined_data, aes_string(x = "cell_type", y = y_var, fill = "group")) +
    geom_bar(stat = "identity", position = "dodge",
             width = 0.8, color = "black", size = 0.3) +
    scale_fill_manual(values = group_colors, labels = group_labels) +
    coord_cartesian(ylim = y_limit) +
    labs(y = y_label, title = plot_title) +
    theme_academic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1))
}

## classify cell types into Excitatory/Inhibitory/Non-neurons by E/I substring match
.classify_cell_types <- function(all_cell_types) {
  exc_neurons <- sort(all_cell_types[grepl("E", all_cell_types)])
  inh_neurons <- sort(all_cell_types[grepl("I", all_cell_types)])
  non_neurons <- sort(setdiff(all_cell_types, c(exc_neurons, inh_neurons)))

  cell_category <- c()
  for (cell in non_neurons) cell_category[cell] <- "Non-neurons"
  for (cell in exc_neurons) cell_category[cell] <- "Excitatory"
  for (cell in inh_neurons) cell_category[cell] <- "Inhibitory"
  cell_category
}

## horizontal diverging bar plot of log2FC(condition/ref)
plot_diff_barplot_2group <- function(df, all_cell_types, condition = "AD", ref = "Control",
                                      threshold = 0.5, x_limit = 1.5, title = "Differential output (log2FC)") {

  cell_category <- .classify_cell_types(all_cell_types)

  all_values <- c(df[[condition]], df[[ref]])
  pseudocount <- 0.1 * median(abs(all_values[all_values != 0])) # 0.1 * median non-zero value

  df_plot <- df %>%
    mutate(
      log_diff_raw = log2((get(condition) + pseudocount) / (get(ref) + pseudocount)),

      exceeds_positive = log_diff_raw > x_limit,
      exceeds_negative = log_diff_raw < -x_limit,

      log_diff_output = pmin(pmax(log_diff_raw, -x_limit), x_limit),
      CellType = sub("\\..*", "", CellType),
      Category = cell_category[CellType],

      bar_color = ifelse(abs(log_diff_raw) >= threshold, "#4472C4", "#D9D9D9"),
      label_text = case_when(
        exceeds_positive ~ paste0(sprintf("%.1f", x_limit), "+"),
        exceeds_negative ~ paste0(sprintf("%.1f", -x_limit), "-"),
        TRUE ~ ""
      )
    ) %>%
    filter(!is.na(log_diff_output) & !is.infinite(log_diff_output)) %>%
    arrange(factor(Category, levels = c("Excitatory", "Inhibitory", "Non-neurons")), CellType) %>%
    mutate(CellType = factor(CellType, levels = CellType))

  df_plot %>%
    ggplot(aes(x = log_diff_output, y = CellType)) +
    geom_col(aes(fill = bar_color), color = "black", linewidth = 0.3, width = 0.8) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed",
               color = "black", linewidth = 0.3, alpha = 0.7) +
    geom_text(data = filter(df_plot, exceeds_positive),
              aes(x = x_limit, y = CellType, label = label_text),
              hjust = -0.15, vjust = 0.5, size = 3.5, color = "black", fontface = "bold") +
    geom_text(data = filter(df_plot, exceeds_negative),
              aes(x = -x_limit, y = CellType, label = label_text),
              hjust = 1.15, vjust = 0.5, size = 3.5, color = "black", fontface = "bold") +
    labs(x = title, y = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 12, face = "italic", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_blank()
    ) +
    scale_x_continuous(limits = c(-x_limit * 1.15, x_limit * 1.15),
                        breaks = seq(-x_limit, x_limit, 0.5))
}

## horizontal grouped bar plot of raw signal strength, AD vs Control - paired with plot_diff_barplot_2group()
plot_signal_barplot_2group <- function(df, all_cell_types, condition = "AD", ref = "Control", xlab = "") {

  cell_category <- .classify_cell_types(all_cell_types)

  df_melt <- reshape2::melt(df, id.vars = "CellType", variable.name = "Group", value.name = "Weight")
  df_filtered <- df_melt[df_melt$Group %in% c(ref, condition), ]

  df_filtered <- df_filtered %>%
    mutate(Category = cell_category[CellType]) %>%
    filter(!is.na(Category)) %>%
    arrange(factor(Category, levels = c("Excitatory", "Inhibitory", "Non-neurons")), CellType) %>%
    mutate(CellType = factor(CellType, levels = unique(CellType)),
           Group = factor(Group, levels = c(condition, ref)))

  two_group_colors <- c("AD" = "#ac0000", "Control" = "grey40")

  ggplot(df_filtered, aes(x = Weight, y = CellType, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black", size = 0.3) +
    scale_fill_manual(values = two_group_colors, name = "Condition") +
    coord_cartesian(xlim = c(0, max(df_filtered$Weight) * 1.1)) +
    labs(x = xlab, y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 12, face = "italic", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title.x = element_text(size = 12, color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_blank(),

      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_blank(),
      legend.key = element_rect(color = "black", size = 0.3),
      legend.key.size = unit(4, "mm"),
      legend.key.width = unit(8, "mm"),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.margin = margin(5, 8, 5, 5),
      legend.box.spacing = unit(4, "mm")
    )
}

## differential chord diagram, adapted from https://github.com/mjgirgenti/PTSDsnDLPFC/blob/main/CCC/netVisual_chord.R
## red = higher in AD, blue = higher in Control
## range_value NULL auto-scales to max(abs(net$prob)) - that's what 11mo used; 8mo used a fixed 0.5, pass range_value = 0.5 to match
my_netVisual_chord <- function(net.diff, edge.transparency = FALSE, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                measure = c("weight", "count"),
                                layout = c("circle", "chord"),
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max = 8,
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8, title.cex = 1.1,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL, cell.order = NULL, small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20, title.name = NULL,
                                range_value = NULL,
                                ...) {
  my_netVisual_chord_cell_internal(net.diff, edge.transparency = edge.transparency, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                    group = group, cell.order = cell.order,
                                    lab.cex = vertex.label.cex, small.gap = small.gap, big.gap = big.gap,
                                    scale = scale, reduce = reduce, title.cex = title.cex,
                                    title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y,
                                    range_value = range_value)
}

my_netVisual_chord_cell_internal <- function(net, edge.transparency = FALSE, color.use = NULL, group = NULL, cell.order = NULL,
                                              sources.use = NULL, targets.use = NULL,
                                              lab.cex = 0.8, small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                              remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                              transparency = 0.4, link.border = NA, title.cex = 1,
                                              title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,
                                              range_value = NULL, ...) {
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source", "target")
  } else if (is.data.frame(net)) {
    if (all(c("source", "target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source, net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  if (!is.null(sources.use)) {
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  net <- subset(net, prob != 0)
  if (dim(net)[1] <= 0) { message("No interaction between those cells") }

  # keep cell types (sectors) with no interactions as isolated, near-zero nodes
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source, net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10 * sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if (nrow(net) > nrow(net.fake)) {
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      scale <- TRUE
    }
  }

  df <- net
  cells.use <- union(df$source, df$target)
  order.sector <- cell.levels[cell.levels %in% cells.use]

  if (is.null(color.use)) {
    color.use <- scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }

  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector

  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }

  # red/white/blue diverging color scale for the AD-minus-Control difference
  if (is.null(range_value)) {
    range_value <- max(abs(net$prob))
  }
  color_ramp <- colorRamp(c("#317EC2", "white", "#C03830"))
  color_ramp2 <- colorRamp2(c(-range_value, 0, range_value), c("#317EC2", "white", "#C03830"))

  edge.color <- rgb(color_ramp((net$prob + range_value) / (2 * range_value)), maxColorValue = 255)
  edge.color2 <- color_ramp2

  rgb_matrix <- col2rgb(edge.color) / 255
  if (edge.transparency) {
    max_prob <- max(abs(df$prob))
    for (i in 1:dim(df)[1]) {
      edge.color[i] <- rgb(rgb_matrix[1, i], rgb_matrix[2, i], rgb_matrix[3, i], abs(df$prob[i] / max_prob))
    }
  }

  link.arr.type <- if (directional == 0 | directional == 2) "triangle" else "big.arrow"

  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = link.arr.type,
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = TRUE,
               group = group,
               link.target.prop = link.target.prop,
               diffHeight = mm_h(4), target.prop.height = mm_h(2),
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = lab.cex)
  }, bg.border = NA)

  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(title = "prob", at = c(-round(range_value, 2), 0, round(range_value, 2)), col_fun = edge.color2)
    ComplexHeatmap::draw(lgd, x = unit(0.97, "npc"), y = unit(legend.pos.y, "mm"), just = c("right"))
  }

  if (!is.null(title.name)) {
    text(0, 0.9, title.name, cex = title.cex)
  }
  circos.clear()
  recordPlot()
}
