#' @title Unsupervised Clinical Text Mining (Latent Semantic Analysis)
#'
#' @description
#' Extracts hidden themes/topics from unstructured clinical free text
#' (e.g., EMR notes, CPPT, or patient satisfaction surveys) using
#' Latent Semantic Analysis (LSA) via Singular Value Decomposition (SVD).
#' Implemented entirely in Base R without external NLP dependencies.
#'
#' @param text_vector Character vector containing the clinical free text documents.
#' @param num_topics Integer. The number of latent topics to extract (default: 2).
#' @param top_terms Integer. The number of top words to display per topic.
#' @param stopwords Character vector of words to exclude (e.g., "dan", "yang", "di").
#' @param plot_result Logical. If TRUE, generates a 1x2 panel barplot of the top 2 topics.
#' @param colors Character vector of length 2 for the color gradient of barplots.
#' @param reverse_palette Logical. If TRUE, reverses the order of the color palette.
#' @param save_plot Logical. If TRUE, saves the generated plot as a high-resolution 4K PNG.
#' @param save_prefix Character. Optional prefix string for the saved plot filename.
#' @param footer Logical. If TRUE, adds a reproducible generation footer to the plot.
#' @param cex Numeric. Text and symbol scaling factor.
#' @param font Numeric. Font style (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' @param family Character. Font family to use for plotting (e.g., "sans", "serif").
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{topic_terms}: A list of data frames showing the top words for each extracted topic.
#'   \item \code{document_topics}: A matrix showing the topic distribution/weight for each document.
#'   \item \code{term_document_matrix}: The preprocessed term-document frequency matrix.
#' }
#'
#' @importFrom graphics par barplot mtext text
#' @importFrom grDevices colorRampPalette png dev.copy dev.off
#'
#' @export
#'
#' @examples
#' # Generate dummy clinical free text (Indonesian EMR & Surveys)
#' dummy_cppt <- c(
#'   "Pasien datang dengan keluhan demam tinggi sejak 3 hari lalu, nyeri sendi hebat, dan ruam merah di kulit. Suspek Dengue.",
#'   "Antrean di farmasi sangat lambat, saya menunggu obat sampai 2 jam. Pelayanan perawat kurang ramah.",
#'   "Demam hari ke-4, trombosit turun, ruam mulai muncul di lengan, pasien mengeluh pusing dan mual.",
#'   "Ruang tunggu poli anak panas sekali, AC mati, dan pelayanan obat lambat.",
#'   "Suhu tubuh 39 derajat, ada bintik merah, ptekie positif, curiga infeksi virus tropis.",
#'   "Sistem antrean BPJS membingungkan, loket pendaftaran antre panjang, obat lama keluar."
#' )
#'
#' # Define basic Indonesian clinical stopwords
#' id_stopwords <- c("dan", "di", "dengan", "ke", "dari", "ini", "itu", "yang",
#'                   "pada", "untuk", "saya", "ada", "sejak", "hari", "lalu")
#'
#' # Run Unsupervised Text Mining
#' res_text <- A5_unsupervised_textmining(
#'   text_vector = df_A5,
#'   num_topics = 2,
#'   top_terms = 5,
#'   stopwords = id_stopwords,
#'   plot_result = TRUE,
#'   save_plot = FALSE
#' )
#'
#' # View top terms for Topic 1
#' print(res_text$topic_terms[[1]])

A5_unsupervised_textmining <- function(text_vector,
                                       num_topics = 2,
                                       top_terms = 10,
                                       stopwords = c("the", "and", "in", "to", "of", "a", "is", "for"),
                                       plot_result = TRUE,
                                       colors = c("#ADD8E6", "#00008B"),
                                       reverse_palette = FALSE,
                                       save_plot = FALSE,
                                       save_prefix = "",
                                       footer = TRUE,
                                       cex = 1,
                                       font = 1,
                                       family = "sans") {

  # 1. Input Validation
  stopifnot("Error: 'text_vector' must be a character vector." = is.character(text_vector))
  stopifnot("Error: 'num_topics' must be at least 1." = num_topics >= 1)

  n_docs <- length(text_vector)
  if (n_docs < 2) {
    stop("Error: Need at least 2 documents to perform text mining.")
  }

  # 2. Native Base R Text Preprocessing
  message("Medical Insight: Preprocessing clinical text (lowercasing, punctuation removal)...")
  clean_text <- tolower(text_vector)
  clean_text <- gsub("[[:punct:]]", " ", clean_text)
  clean_text <- gsub("[[:digit:]]", " ", clean_text)

  # Split into words
  words_list <- strsplit(clean_text, "\\s+")

  # Flatten and filter stopwords/empty strings
  all_words <- unlist(words_list)
  all_words <- all_words[all_words != "" & !(all_words %in% stopwords)]
  unique_words <- unique(all_words)

  if (length(unique_words) < num_topics) {
    stop("Error: The number of unique words after cleaning is smaller than the requested num_topics.")
  }

  # 3. Build Term-Document Matrix (TDM)
  tdm <- matrix(0, nrow = length(unique_words), ncol = n_docs)
  rownames(tdm) <- unique_words

  for (i in seq_along(words_list)) {
    doc_words <- words_list[[i]]
    doc_words <- doc_words[doc_words != "" & !(doc_words %in% stopwords)]
    if (length(doc_words) > 0) {
      word_counts <- table(doc_words)
      tdm[names(word_counts), i] <- as.numeric(word_counts)
    }
  }

  # Remove extremely sparse terms (words that only appear once across all docs)
  # to improve Topic Modeling stability
  term_freq <- rowSums(tdm)
  tdm <- tdm[term_freq > 1, , drop = FALSE]

  if (nrow(tdm) < num_topics) {
    stop("Error: Not enough overlapping terms to build reliable topics. Try reducing num_topics or expanding your dataset.")
  }

  # 4. Core Algorithm: Singular Value Decomposition (SVD) for LSA
  # scale = FALSE to preserve sparse matrix characteristics
  s_decomp <- svd(tdm, nu = num_topics, nv = num_topics)

  # U matrix: Term-Topic associations
  term_topic_mat <- s_decomp$u
  rownames(term_topic_mat) <- rownames(tdm)

  # V matrix: Document-Topic associations
  doc_topic_mat <- s_decomp$v

  # Extract top terms for each topic
  topic_terms_list <- list()
  for (k in 1:num_topics) {
    # Get absolute weights to capture semantic magnitude
    term_weights <- abs(term_topic_mat[, k])
    sorted_idx <- order(term_weights, decreasing = TRUE)

    top_n_idx <- sorted_idx[1:min(top_terms, length(sorted_idx))]

    topic_terms_list[[k]] <- data.frame(
      Term = rownames(term_topic_mat)[top_n_idx],
      Weight = term_weights[top_n_idx],
      stringsAsFactors = FALSE
    )
  }
  names(topic_terms_list) <- paste0("Topic_", 1:num_topics)

  # 5. Strict Visualization (1x2 Panel Plot for Top 2 Topics)
  if (plot_result && num_topics >= 2) {

    if (reverse_palette) {
      colors <- rev(colors)
    }

    # Save old par
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))

    # Create global layout with oma for footer and prevent overlap
    graphics::par(mfrow = c(1, 2),
                  oma = c(4, 0, 2, 0), # Top margin for global title, bottom for footer
                  family = family, font = font, cex = cex)

    color_pal <- grDevices::colorRampPalette(c(colors[1], colors[2]))(top_terms)

    for (k in 1:2) { # We strictly plot the first 2 topics for visual balance
      # Dynamic margin for each subplot: Left margin (9) to fit long clinical words
      graphics::par(mar = c(5, 9, 2, 2) + 0.1)

      plot_df <- topic_terms_list[[k]]
      # Reverse order so highest is on top of barplot
      plot_df <- plot_df[order(plot_df$Weight, decreasing = FALSE), ]

      bp <- graphics::barplot(plot_df$Weight,
                              names.arg = plot_df$Term,
                              horiz = TRUE,
                              las = 1,
                              col = color_pal,
                              border = "white",
                              xlab = "Semantic Weight",
                              main = paste("Latent Topic", k))

      # Add text label for values
      graphics::text(x = 0, # align left
                     y = bp,
                     labels = paste0(" ", round(plot_df$Weight, 3)),
                     pos = 4,
                     col = "white", # White text inside the colored bar
                     cex = 0.8,
                     font = 2)
    }

    # Global Title
    graphics::mtext("Unsupervised Clinical Text Mining (LSA)", outer = TRUE, side = 3, line = 0, cex = 1.2, font = 2)

    # Global Footer
    if (footer) {
      footer_text <- paste0("Generated by R-Studio (", R.version.string, ") on ", format(Sys.time(), "%B %d, %Y at %H:%M:%S"))
      graphics::mtext(footer_text, side = 1, line = 2, outer = TRUE, adj = 0.5, cex = 0.8, col = "dimgray")
    }

    # High-Res 4K Export
    if (save_plot) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      clean_prefix <- ifelse(nchar(save_prefix) > 0, paste0("_", save_prefix), "")
      filename <- paste0("A5_unsupervised_textmining_", timestamp, clean_prefix, ".png")

      grDevices::dev.copy(grDevices::png, filename = filename, width = 3840, height = 2160, res = 300)
      grDevices::dev.off()

      message(paste("High-resolution 4K plot saved as:", filename))
    }
  } else if (plot_result && num_topics < 2) {
    warning("Medical Insight: At least 2 topics are required for the standard 1x2 panel plot visualization.")
  }

  # 6. Return Data
  result_list <- list(
    topic_terms = topic_terms_list,
    document_topics = doc_topic_mat,
    term_document_matrix = tdm
  )

  return(result_list)
}
