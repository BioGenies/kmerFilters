library(dplyr)

df <- read.csv("./data/dbaasp_amp_acp.csv")

df %>% filter(ACP = TRUE) -> df_acp


# first approach - AA frequences in a concatenated sequence
df_acp[["SEQUENCE"]]
st(paste0(df_acp[["SEQUENCE"]], collapse=""))

