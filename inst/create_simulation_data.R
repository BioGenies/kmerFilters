#' function generates set of RDS files containing n-gram matrices
#' @param replications number of matrices' generation for each set of parameters
#' @param seq_nums vector of sequences' numbers
#' @param seq_lengths vector of sequences' lengths
#' @param num_motifs vector of number of motifs to be injected
#' @param motif_set_size number of motifs to be created
#' @param alphabet elements used to build both sequence and motif
#' @param path file path of RDS files
#' @param title simulation name
#' @param save_files triggers if RDS will be saved
#' @param motifProbs alphabet probabilites for motifs
#' @param seqProbs alphabet probabilites for sequences
#' @param fraction fraction of positive sequences
#' @param n maximum number of alphabet elements in motif
#' @param d maximum number of gaps in motif
#' @return data frame with RDS files' details
#' @importFrom utils write.csv
#' @importFrom progress progress_bar
#' @importFrom tools md5sum
#' @examples
#' alph <- letters[1:4]
#' reps <- 10
#' n_seq <- c(10)
#' l_seq <- c(10)
#' n_motifs <- c(1, 2)
#' path = "./"
#' results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alph, path, "SEQ",FALSE)
#'
#' alph <- letters[1:4]
#' reps <- 3
#' n_seq <- c(10)
#' l_seq <- c(10)
#' n_motifs <- 1
#' results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alph,
#'                                  path, "SEQ", FALSE,
#'                                  motifProbs = c(0.7, 0.1, 0.1, 0.1),
#'                                  seqProbs = c(0.7, 0.1, 0.1, 0.1),
#'                                  n = 4, d = 4)
#' @export
create_simulation_data <- function(replications,
                                   seq_nums,
                                   seq_lengths,
                                   motif_set_size,
                                   num_motifs,
                                   alphabet,
                                   path,
                                   title,
                                   save_files = FALSE,
                                   motifProbs = NULL,
                                   seqProbs = NULL,
                                   fraction = 0.5,
                                   n = 4,
                                   d = 6) {
  
  if (!file.exists(path))
    stop("Output directory does not exist.")
  
  df <- data.frame(replication = integer(),
                   n_seq = integer(),
                   l_seq = integer(),
                   n_motifs = integer(),
                   n = integer(),
                   d = integer(),
                   seqProbs = character(),
                   motifProbs = character(),
                   path = character(),
                   md5sum = character())
  
  totalIterations = replications * length(seq_nums) * length(seq_lengths) *  length(num_motifs)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = totalIterations)
  pb$tick(0)
  
  for (replication in 1:replications) {
    for (n_motifs in num_motifs) {
      for (n_seq in seq_nums) {
        for (sequence_length in seq_lengths) {
          
          pb$tick(1)
          set.seed(replication)
          
          # 90% of all combinations of 
          #possible_90perc <- FALSE
          
          #while (!possible_90perc) {
            
            motifs <- lapply(1L:motif_set_size, 
                             function(dummy) generate_motif(alphabet, n, d, motifProbs))     
            
            motifs_grid <- do.call(expand.grid, rep(list(1:motif_set_size), n_motifs))
            
            possible_motifs_grid <- motifs_grid[apply(motifs_grid, 1, 
                                                      function(x) length(unique(unlist(x))) == n_motifs), , drop=FALSE]
            
            #check <- table(apply(possible_motifs_grid, 1, function(x) validate_motifs(motifs[unlist(x)], 10)))
            
            #if (length(check) == 1) {
            #  possible_90perc <- TRUE
            #} else if (check["FALSE"] * 9 <= check["TRUE"]) {
            #  possible_90perc <- TRUE              
            #}
          #}
          
          dat <- generate_kmer_data(n_seq,
                                    sequence_length,
                                    alphabet,
                                    motifs,
                                    n_motifs,
                                    fraction = fraction,
                                    seqProbs = seqProbs,
                                    n = n,
                                    d = d)
          
          filePath = paste0(path, title,
                            "_rep_", replication,
                            "_seqNum_", n_seq,
                            "_seqLen_", l_seq,
                            "_nMotifs_", n_motifs,
                            ".Rds")
          
          if (save_files)
            saveRDS(dat, file = filePath)
          
          md5 <- ifelse(save_files, md5sum(filePath), "none")
          
          newdf <- data.frame(replication = replication,
                              n_seq = n_seq,
                              l_seq = l_seq,
                              n_motifs = n_motifs,
                              n = n,
                              d = d,
                              seqProbs = paste(motifProbs, collapse = ";"),
                              motifProbs = paste(motifProbs, collapse = ";"),
                              path = filePath,
                              md5sum = md5)
          
          df <- rbind(df, newdf)
          
        }
      }
    }
  }
  
  if (save_files)
    write.csv(df, file = paste0(path, title,".csv"))
  df
}