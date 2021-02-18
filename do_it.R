# Goal: find pattern in flanking regions presented in MHC-I elution study
#
# Unclear if this is for all epitopes or only TMH-derived epitopes
#


library(stringr, quietly = TRUE, warn.conflicts = FALSE)
library(english, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

#' Analyse the epitopes
#' @param epitope_sequences sequences of epitopes
#' @param positions string positions as used by \link[stringr]{str_sub}, e.g. '1' denotes the first character,
#'   '-1' means the last character, '-2' the one-but-last
analyse_epitopes <- function(
  epitope_sequences,
  positions = c(1, 2, 3, -3, -2, -1)
) {
  # Store sequences as tibble, as stringr needs it
  char_tibble <- tibble::tibble(
    epitope_sequence = epitope_sequences
  )

  # Collect all tibbles in this list
  tibbles <- list()

  # Count all chararacters at all positions
  for (i in seq_along(positions)) {
    pos <- positions[i]
    char_tibble$char <- stringr::str_sub(char_tibble$epitope_sequence, pos, pos)
    t <- dplyr::count(char_tibble, char)
    t$pos <- pos
    t <- dplyr::relocate(t, pos)
    tibbles[[i]] <- t
  }
  # Combine the reults per position
  t_all <- dplyr::bind_rows(tibbles)

  # Convert the data to wide format for Libreoffice Calc users
  t_wide <- tidyr::pivot_wider(t_all, names_from = pos, values_from = n) %>% dplyr::arrange(char)

  # Make the column names reader friendly
  names(t_wide) <- c(
    names(t_wide)[1],
    stringr::str_replace(english::as.english(as.numeric(names(t_wide)[2:7])), " ", "_")
  )

  return(t_wide)
}

t_all_unique_matches <- readr::read_csv("matches_1.csv") %>% filter(n_matches == 1)

png("all_presented_epitope_lengths.png")
hist(nchar(t_all_unique_matches$epitope_sequence), freq = TRUE, breaks = seq(6, 25), main = "Presented MHC-I epitope lengths")
dev.off()

t_wide <- analyse_epitopes(epitope_sequences = t_all_unique_matches$epitope_sequence)
t_wide
readr::write_csv(t_wide, "all.csv")

#
# For TMH epitopes
#

t_unique_tmhs_tmhmm_1 <- readr::read_csv("tmhs_tmhmm_1.csv") %>% filter(n_matches == 1)
png("tmhs_presented_epitope_lengths.png")
hist(nchar(t_unique_tmhs_tmhmm_1$epitope_sequence), freq = TRUE, breaks = seq(7, 24), main = "Presented MHC-I TMH epitopes' lengths")
dev.off()


t_wide <- analyse_epitopes(epitope_sequences = t_unique_tmhs_tmhmm_1$epitope_sequence)
t_wide
readr::write_csv(t_wide, "tmh.csv")


