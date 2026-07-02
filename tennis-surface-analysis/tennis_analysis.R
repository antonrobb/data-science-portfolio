# Tennis Majors 2013 (UCI) — Surface Effects on Playing Style

# Data source: UCI Machine Learning Repository —
# "Tennis Major Tournament Match Statistics" (2013 Grand Slams)

# I use the 2013 Grand Slam match statistics to explore how court
# surface, gender, and tournament round shape playing style.

# 0) Packages
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)

# 1) Data dictionary (abbreviations used throughout)
# I've only listed the variables I actually use below. If you're
# extending this analysis, run `names(tennis)` after loading to see
# the full column set and add to this table as needed.
data_dictionary <- tibble::tribble(
  ~Abbreviation, ~Meaning,
  "ACE",  "Aces served",
  "FSP",  "First serve percentage (%)",
  "FSW",  "First serve points won",
  "TPW",  "Total points won",
  "UFE",  "Unforced errors committed",
  "WNR",  "Winners hit"
)
print(data_dictionary)

# 2) Read + combine all CSV files
files <- c(
  "AusOpen-men-2013.csv",
  "AusOpen-women-2013.csv",
  "FrenchOpen-men-2013.csv",
  "FrenchOpen-women-2013.csv",
  "Wimbledon-men-2013.csv",
  "Wimbledon-women-2013.csv",
  "USOpen-men-2013.csv",
  "USOpen-women-2013.csv"
)

parse_file_info <- function(fname) {
  tour <- str_replace(fname, "-(men|women)-2013\\.csv$", "")
  gender <- str_match(fname, "-(men|women)-2013\\.csv$")[, 2]
  list(tournament = tour, gender = str_to_title(gender))
}

tennis <- lapply(files, function(f) {
  info <- parse_file_info(f)
  df <- read_csv(f, show_col_types = FALSE)
  df$tournament <- info$tournament
  df$gender <- info$gender
  df$file <- f
  df
}) %>% bind_rows()

# 3) Dataset size
cat(sprintf(
  "Dataset: %d matches, %d columns, across %d tournament-gender files.\n",
  nrow(tennis), ncol(tennis), length(files)
))
print(tennis %>% count(tournament, gender, name = "n_matches"))

# 4) Minimal cleaning + derived fields
# (a) Outcome label (if Result exists). Result == 1 means Player 1 won.
if ("Result" %in% names(tennis)) {
  tennis <- tennis %>%
    mutate(
      Result = as.integer(Result),
      Outcome = factor(Result, levels = c(0, 1), labels = c("P1 Loss", "P1 Win"))
    )
}

# (b) Factors
tennis <- tennis %>%
  mutate(
    tournament = factor(tournament, levels = c("AusOpen", "FrenchOpen", "USOpen", "Wimbledon")),
    gender = factor(gender, levels = c("Men", "Women"))
  )

# (c) Surface derived from tournament
tennis <- tennis %>%
  mutate(
    Surface = case_when(
      tournament == "FrenchOpen" ~ "Clay",
      tournament == "Wimbledon" ~ "Grass",
      tournament %in% c("AusOpen", "USOpen") ~ "Hard",
      TRUE ~ NA_character_
    ),
    Surface = factor(Surface, levels = c("Clay", "Grass", "Hard"))
  )

# (d) Round as ordered factor (if present)
if ("Round" %in% names(tennis)) {
  tennis <- tennis %>% mutate(Round = factor(Round, ordered = TRUE))
}

# 5) Row-level missingness
tennis <- tennis %>%
  mutate(missing_count = rowSums(is.na(.)))

# 6) Identify paired P1/P2 stat columns (e.g. ACE.1/ACE.2)
p1_cols <- names(tennis)[str_detect(names(tennis), "\\.1$|_1$")]
p2_cols <- names(tennis)[str_detect(names(tennis), "\\.2$|_2$")]
base_p1 <- str_replace(p1_cols, "(\\.1$|_1$)", "")
base_p2 <- str_replace(p2_cols, "(\\.2$|_2$)", "")
paired_bases <- intersect(base_p1, base_p2)

get_col <- function(df, base, suffix_num) {
  dot_name <- paste0(base, ".", suffix_num)
  us_name  <- paste0(base, "_", suffix_num)
  if (dot_name %in% names(df)) dot_name else us_name
}

# 7) Recode P1/P2 -> Winner/Loser
# "P1" and "P2" in the raw data are arbitrary labels that don't track
# who actually won. That means a naive column like WNR_diff (winners
# hit, P1-P2) mixes together "the winner's advantage" and "the loser's
# disadvantage" depending on who happened to be listed as P1 in that
# row. I recode so Player 1 = winner and Player 2 = loser, so every
# "_margin" column I build from here on consistently means "how much
# better the winner performed than the player they beat", which is a
# much more interpretable quantity than a raw P1-minus-P2 difference.
if ("Result" %in% names(tennis)) {
  for (b in paired_bases) {
    c1 <- get_col(tennis, b, 1)
    c2 <- get_col(tennis, b, 2)
    if (c1 %in% names(tennis) && c2 %in% names(tennis)) {
      v1 <- suppressWarnings(as.numeric(tennis[[c1]]))
      v2 <- suppressWarnings(as.numeric(tennis[[c2]]))
      tennis[[paste0(b, "_W")]] <- ifelse(tennis$Result == 1, v1, v2) # winner's stat
      tennis[[paste0(b, "_L")]] <- ifelse(tennis$Result == 1, v2, v1) # loser's stat
      tennis[[paste0(b, "_margin")]] <- tennis[[paste0(b, "_W")]] - tennis[[paste0(b, "_L")]]
    }
  }
} else {
  message("Result column not found — cannot recode to Winner/Loser; falling back to raw P1-P2 diffs.")
  for (b in paired_bases) {
    c1 <- get_col(tennis, b, 1)
    c2 <- get_col(tennis, b, 2)
    if (c1 %in% names(tennis) && c2 %in% names(tennis)) {
      tennis[[paste0(b, "_margin")]] <- suppressWarnings(as.numeric(tennis[[c1]]) - as.numeric(tennis[[c2]]))
    }
  }
}

# 8) Clean analysis frame (known Surface)
tennis_clean <- tennis %>% filter(!is.na(Surface))

# 9) Player-level long data for ACE/FSP/UFE/WNR (winner + loser rows)
player_stats <- c("ACE", "FSP", "UFE", "WNR")
w_cols <- paste0(player_stats, "_W")
l_cols <- paste0(player_stats, "_L")
player_cols <- c(rbind(w_cols, l_cols))
player_cols <- player_cols[player_cols %in% names(tennis_clean)]

tennis_long <- tennis_clean %>%
  pivot_longer(
    cols = all_of(player_cols),
    names_to = c("stat", "player"),
    names_pattern = "(.*)_(W|L)",
    values_to = "value"
  ) %>%
  mutate(
    player = factor(player, levels = c("W", "L"), labels = c("Winner", "Loser")),
    stat = factor(stat, levels = player_stats)
  )

# 10) Consistent colour palettes
surface_cols <- c("Clay" = "#F8766D", "Grass" = "#00BA38", "Hard" = "#619CFF")
tourn_cols <- c(
  "AusOpen" = "steelblue4",
  "FrenchOpen" = "firebrick3",
  "USOpen" = "goldenrod2",
  "Wimbledon" = "seagreen4"
)


# FIGURE 0a — How many matches per tournament/gender
fig0a_match_counts <- ggplot(
  tennis_clean %>% count(tournament, gender),
  aes(x = tournament, y = n, fill = gender)
) +
  geom_col(position = "dodge") +
  labs(
    title = "Number of Matches by Tournament and Gender",
    x = "Tournament", y = "Number of matches", fill = "Gender"
  ) +
  theme_minimal()

fig0a_match_counts

# FIGURE 0b — Overall distribution of a key stat (Aces)
fig0b_ace_hist <- ggplot(tennis_long %>% filter(stat == "ACE"), aes(x = value)) +
  geom_histogram(binwidth = 2, fill = "steelblue4", colour = "white") +
  labs(
    title = "Distribution of Aces per Player per Match (all surfaces)",
    x = "Aces", y = "Count"
  ) +
  theme_minimal()

fig0b_ace_hist


# FIGURE 1 — Missing data boxplot split by tournament
fig1_missing_by_tour <- ggplot(
  tennis_clean,
  aes(x = reorder(tournament, missing_count, FUN = median),
      y = missing_count,
      fill = tournament)
) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_manual(values = tourn_cols) +
  labs(
    title = "Distribution of Missing Fields per Match by Tournament",
    x = "Tournament",
    y = "Number of missing fields per match"
  ) +
  theme_minimal() +
  guides(fill = "none")

fig1_missing_by_tour


# FIGURE 2 — Playing style by surface & gender (violin + box)
# Winner vs Loser instead of arbitrary P1/P2.
# I trim the y-axis to the 1st-99th percentile per stat for
# readability, there were a few extreme outliers were stretching the
# violins out so far that the bulk of the distribution was hard to see.
# This only affects the plotted range, not the data used elsewhere.
trim_limits <- tennis_long %>%
  group_by(stat) %>%
  summarise(
    lo = quantile(value, 0.01, na.rm = TRUE),
    hi = quantile(value, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

fig2_violin_style_surface_gender <- ggplot(tennis_long, aes(x = Surface, y = value, fill = Surface)) +
  geom_violin(trim = TRUE, alpha = 0.55, colour = NA) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.2, fill = "white") +
  facet_grid(stat ~ gender, scales = "free_y") +
  scale_fill_manual(values = surface_cols, guide = "none") +
  labs(
    title = "Playing Style Statistics by Surface and Gender (Winner + Loser)",
    subtitle = "Y-axis trimmed to 1st-99th percentile per statistic for readability",
    x = "Court Surface",
    y = "Count (per player per match)"
  ) +
  theme_minimal()

fig2_violin_style_surface_gender <- fig2_violin_style_surface_gender +
  geom_blank(data = trim_limits %>% tidyr::pivot_longer(c(lo, hi), values_to = "value") %>%
               mutate(Surface = tennis_long$Surface[1]),
             aes(x = Surface, y = value))

fig2_violin_style_surface_gender


# FIGURE 3a — Winning margin: winners vs errors
# One overall panel, no facets, no annotation. I like to look at the
# headline relationship before splitting it apart by tournament.
if (all(c("WNR_margin", "UFE_margin") %in% names(tennis_clean))) {
  
  fig3a_winners_errors_simple <- ggplot(tennis_clean, aes(x = WNR_margin, y = UFE_margin)) +
    geom_point(alpha = 0.4, colour = "steelblue4") +
    geom_smooth(method = "lm", se = FALSE, colour = "black") +
    labs(
      title = "Winner's Margin in Winners vs Unforced Errors (all matches)",
      x = "Winners margin (winner - loser)",
      y = "Unforced errors margin (winner - loser)"
    ) +
    theme_minimal()
  
  fig3a_winners_errors_simple
  
  # FIGURE 3b — Same relationship, now faceted by tournament with
  # slope / R^2 / n annotated per panel, so I can compare how tightly
  # winners and errors track each other on each surface.
  slopes_df <- tennis_clean %>%
    filter(!is.na(WNR_margin), !is.na(UFE_margin), !is.na(tournament)) %>%
    group_by(tournament) %>%
    do({
      m <- lm(UFE_margin ~ WNR_margin, data = .)
      tibble(
        slope = coef(m)[["WNR_margin"]],
        r2 = summary(m)$r.squared,
        n = nrow(.)
      )
    }) %>%
    ungroup()
  
  label_pos <- tennis_clean %>%
    filter(!is.na(WNR_margin), !is.na(UFE_margin), !is.na(tournament)) %>%
    group_by(tournament) %>%
    summarise(
      x = min(WNR_margin, na.rm = TRUE),
      y = max(UFE_margin, na.rm = TRUE),
      .groups = "drop"
    )
  
  slopes_annot <- slopes_df %>%
    left_join(label_pos, by = "tournament") %>%
    mutate(label = paste0("slope = ", round(slope, 3), "\nR\u00b2 = ", round(r2, 3), "\nn = ", n))
  
  fig3b_winners_errors_by_tour <- ggplot(tennis_clean, aes(x = WNR_margin, y = UFE_margin)) +
    geom_point(alpha = 0.5, colour = "steelblue4") +
    geom_smooth(method = "lm", se = FALSE, colour = "black") +
    facet_wrap(~ tournament) +
    geom_text(
      data = slopes_annot,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3.5,
      colour = "black"
    ) +
    labs(
      title = "Winner's Margin in Winners vs Unforced Errors, by Tournament",
      subtitle = "Facet labels show slope (gradient), R\u00b2, and n. Positive slope = winners hit more\nwinners AND made fewer errors than the player they beat.",
      x = "Winners margin (winner - loser)",
      y = "Unforced errors margin (winner - loser)"
    ) +
    theme_minimal()
  
  fig3b_winners_errors_by_tour
  
} else {
  message("Fig 3 skipped: WNR_margin and/or UFE_margin not found. Check column names.")
}


# FIGURE 4 — Match advantages by round and tournament
# I drop the US Open here: its round information is much sparser than
# the other three majors (see the missingness section above), and
# including it produces thin, unstable-looking boxplots that don't
# add much.
trend_vars <- c("ACE_margin", "FSP_margin", "FSW_margin", "TPW_margin", "UFE_margin", "WNR_margin")
trend_vars <- trend_vars[trend_vars %in% names(tennis_clean)]

tennis_round <- tennis_clean %>% filter(tournament != "USOpen")

if ("Round" %in% names(tennis_round) && length(trend_vars) > 0) {
  trend_long <- tennis_round %>%
    select(tournament, Round, all_of(trend_vars)) %>%
    pivot_longer(cols = all_of(trend_vars), names_to = "metric", values_to = "margin_value") %>%
    mutate(metric = factor(metric, levels = trend_vars))
  
  round_labels <- function(x) {
    dplyr::recode(x, "5" = "Quarter", "6" = "Semi", "7" = "Final", .default = x)
  }
  
  fig4_advantages_by_round <- ggplot(trend_long, aes(x = Round, y = margin_value)) +
    geom_boxplot(aes(colour = tournament), fill = NA, outlier.alpha = 0.3) +
    facet_grid(metric ~ tournament, scales = "free_y") +
    scale_colour_manual(values = tourn_cols, guide = "none") +
    scale_x_discrete(labels = round_labels) +
    labs(
      title = "How the Winner's Margin Varies by Round (US Open excluded)",
      x = "Round",
      y = "Winner's margin (winner - loser)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  fig4_advantages_by_round
  
  # FIGURE 4b — Standard deviation of winner's margin by round
  # The boxplots above hint that matches get closer in later rounds,
  # this makes that trend explicit as a single summary line per metric.
  sd_by_round <- trend_long %>%
    group_by(tournament, Round, metric) %>%
    summarise(sd_value = sd(margin_value, na.rm = TRUE), .groups = "drop")
  
  fig4b_sd_by_round <- ggplot(sd_by_round, aes(x = Round, y = sd_value, colour = tournament, group = tournament)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_colour_manual(values = tourn_cols) +
    scale_x_discrete(labels = round_labels) +
    labs(
      title = "Standard Deviation of Winner's Margin by Round (US Open excluded)",
      subtitle = "A falling line means match outcomes become more consistent (less variable) in later rounds",
      x = "Round",
      y = "SD of winner's margin",
      colour = "Tournament"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  fig4b_sd_by_round
  
} else {
  message("Fig 4 skipped: Round not found or no *_margin columns available.")
}


# FIGURE 5 — Mean player stats by surface and gender (Winner + Loser)
mean_df <- tennis_long %>%
  group_by(gender, Surface, stat, player) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

fig5_mean_player_stats <- ggplot(mean_df, aes(x = Surface, y = mean_value, group = interaction(gender, player), colour = gender, linetype = player)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ stat, scales = "free_y") +
  labs(
    title = "Mean Player Statistics by Surface, Gender, and Winner/Loser",
    x = "Court Surface",
    y = "Mean Count",
    colour = "Gender",
    linetype = "Player"
  ) +
  theme_minimal()

fig5_mean_player_stats


# Save
# ggsave("Fig0a_match_counts.png", fig0a_match_counts, width = 8, height = 4, dpi = 300)
# ggsave("Fig0b_ace_histogram.png", fig0b_ace_hist, width = 8, height = 4, dpi = 300)
# ggsave("Fig1_missing_by_tournament.png", fig1_missing_by_tour, width = 12, height = 4, dpi = 300)
# ggsave("Fig2_violin_style_surface_gender.png", fig2_violin_style_surface_gender, width = 10, height = 8, dpi = 300)
# ggsave("Fig3a_winners_errors_simple.png", fig3a_winners_errors_simple, width = 7, height = 5, dpi = 300)
# ggsave("Fig3b_winners_errors_by_tournament.png", fig3b_winners_errors_by_tour, width = 10, height = 5, dpi = 300)
# ggsave("Fig4_advantages_by_round.png", fig4_advantages_by_round, width = 12, height = 8, dpi = 300)
# ggsave("Fig4b_sd_by_round.png", fig4b_sd_by_round, width = 10, height = 6, dpi = 300)
# ggsave("Fig5_mean_player_stats.png", fig5_mean_player_stats, width = 10, height = 5, dpi = 300)