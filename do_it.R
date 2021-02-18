# I love the descriptiveness of the name :-)


library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

t_all_unique_matches <- readr::read_csv("matches_1.csv") %>% filter(n_matches == 1)

png("all_presented_epitope_lengths.png")
hist(nchar(t_all_unique_matches$epitope_sequence), freq = TRUE, breaks = seq(6, 25), main = "Presented MHC-I epitope lengths")
dev.off()

positions <- c(1, 2, 3, -3, -2, -1)
tibbles <- list()
# I see the epitope length of 8 is most often present
# Flanking regions: pos 1, 2, 3, two but last, one but last, last
for (i in seq_along(positions)) {
  pos <- positions[i]
  t_all_unique_matches$char <- stringr::str_sub(t_all_unique_matches$epitope_sequence, pos, pos)
  t <- dplyr::count(t_all_unique_matches, char)
  t$pos <- pos
  t <- dplyr::relocate(t, pos)
  tibbles[[i]] <- t
}
t_all <- dplyr::bind_rows(tibbles)
t_wide <- tidyr::pivot_wider(t_all, names_from = pos, values_from = n) %>% dplyr::arrange(char)
names(t_wide) <- c(
  names(t_wide)[1],
  stringr::str_replace(english::as.english(as.numeric(names(t_wide)[2:7])), " ", "_")
)
t_wide
readr::write_csv(t_wide, "all.csv")

#
# For TMH epitopes
#

t_unique_tmhs_tmhmm_1 <- readr::read_csv("tmhs_tmhmm_1.csv") %>% filter(n_matches == 1)
png("tmhs_presented_epitope_lengths.png")
hist(nchar(t_unique_tmhs_tmhmm_1$epitope_sequence), freq = TRUE, breaks = seq(7, 24), main = "Presented MHC-I TMH epitopes' lengths")
dev.off()

# I see the epitope length of 8 is most often present
# Flanking regions: pos 1, 2, 3, two but last, one but last, last
positions <- c(1, 2, 3, -3, -2, -1)
tibbles <- list()
# I see the epitope length of 8 is most often present
# Flanking regions: pos 1, 2, 3, two but last, one but last, last
for (i in seq_along(positions)) {
  pos <- positions[i]
  t_unique_tmhs_tmhmm_1$char <- stringr::str_sub(t_unique_tmhs_tmhmm_1$epitope_sequence, pos, pos)
  t <- dplyr::count(t_unique_tmhs_tmhmm_1, char)
  t$pos <- pos
  t <- dplyr::relocate(t, pos)
  tibbles[[i]] <- t
}
t_all <- dplyr::bind_rows(tibbles)
t_wide <- tidyr::pivot_wider(t_all, names_from = pos, values_from = n) %>% dplyr::arrange(char)
names(t_wide) <- c(
  names(t_wide)[1],
  stringr::str_replace(english::as.english(as.numeric(names(t_wide)[2:7])), " ", "_")
)
t_wide
readr::write_csv(t_wide, "tmh.csv")


