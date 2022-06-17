library(dplyr)
library(ggplot2)
library(stringr)
df <- read.csv("./data/dbaasp_amp_acp.csv")

df %>% filter(ACP = TRUE) -> df_acp


# first approach - AA frequences in a concatenated sequence
df_acp[["SEQUENCE"]]
table(strsplit(paste0(df[["SEQUENCE"]], collapse=""), split=""))
freqs <- table(strsplit(toupper(paste0(df[["SEQUENCE"]], collapse="")), split=""))
freqs <- data.frame(freqs)
ggplot(data=freqs, aes(x=Var1, y=Freq)) + geom_bar(stat='identity')

# X, przecinki, spacje
spaces_in_seq <- unlist(lapply(df[["SEQUENCE"]], function(x) grepl(" ", x)))
commas_in_seq <- unlist(lapply(df[["SEQUENCE"]], function(x) grepl(",", x)))
head(df[spaces_in_seq,])
head(df[commas_in_seq,])

num_unique_chars <- unlist(lapply(df[["SEQUENCE"]], 
                                  function(x) length(table(strsplit(x, "")))))
head(df[num_unique_chars == 1, ])
