count_parsimony <- function(tree, asr, occurrence = TRUE, transition_n = FALSE, type = c("all", "tips", "inner_nodes")){
  if(any(apply(asr, 1, sum) != 1)){
    stop("One or more states are ambigous!")
  }
  if(is.matrix(asr)){
    asr <- as.data.frame(asr)
  }else{
    if(!is.data.frame(asr)){
      stop("asr must be a matrix or a data.frame")
    }
  }
  asr_nstates <- ncol(asr)
  states <- 1:asr_nstates
  asr_nrow <- nrow(asr)
  if(asr_nrow != max(tree$edge)){
    stop("asr must be with a row number = # of tree nodes + # of tree tips")
  }
  asr$states <- rep(NA, asr_nrow)
  for(i in 1:asr_nstates){
    asr$states[asr[,i] == 1] <- i
  }
  asr$nodes <- 1:asr_nrow
  tip_num <- 1:length(tree$tip.label)
  switch(type,
         all = {node_df <- as.data.frame(tree$edge)},
         tips = {node_df <- as.data.frame(tree$edge[tree$edge[,2] %in% tip_num,])},
         inner_nodes = {node_df <- as.data.frame(tree$edge[!tree$edge[,2] %in% tip_num,])})
  node_df$left <- rep(NA, nrow(node_df))
  node_df$right <- rep(NA, nrow(node_df))
  node_df$left <- asr$states[match(node_df$V1, asr$nodes)]
  node_df$right <- asr$states[match(node_df$V2, asr$nodes)]
  transitions <- expand.grid(states, states, stringsAsFactors = FALSE)
  transitions <- transitions[transitions[,1] != transitions[,2],]
  for(i in 1:nrow(transitions)){
    trans <- node_df$left == transitions[i, 1] & node_df$right == transitions[i, 2]
    node_df <- cbind(node_df, trans)
  }
  colnames(node_df)[5:(4+nrow(transitions))] <- str_c(transitions[,1], "_", transitions[,2])
  transition_count <- apply(node_df[,5:(4+nrow(transitions))], 2, sum)
  names(transition_count) <- str_c(transitions[,1], "_", transitions[,2])
  occurrence_count <- vector()
  for(i in 1:asr_nstates){
    occurrence_count[i] <- sum(transition_count[which(transitions[,2] == i)])
  }
  names(occurrence_count) <- str_c("state_", 1:asr_nstates)
  if(occurrence && transition_n == FALSE){
    return(occurrence_count)
  }
  if(occurrence && transition_n){
    asr_result <- list(occurrence_count, transition_count)
    return(asr_result)
  }
  if(occurrence == FALSE && transition_n){
    return(transition_count)
  }
}
