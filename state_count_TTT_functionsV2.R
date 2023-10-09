TTT_plot <- function(corHMM = NULL, window_width = 12, step = 6, trait = NULL) {
  
  ace_ML <- corhmm_to_ace(corHMM)
  asr <- as.data.frame(ace_ML)
  asr_nstates <- ncol(asr)
  states <- 1:asr_nstates
  asr_nrow <- nrow(asr)
  asr$states <- rep(NA, asr_nrow)
  for(i in 1:asr_nstates){
    asr$states[asr[,i] == 1] <- i
  }
  asr$nodes <- 1:asr_nrow
  
  tree <- corHMM$phy
  node_df <- as.data.frame(tree$edge)
  node_df$left <- rep(NA, nrow(node_df))
  node_df$right <- rep(NA, nrow(node_df))
  node_df$left <- asr$states[match(node_df$V1, asr$node)]
  node_df$right <- asr$states[match(node_df$V2, asr$node)]
  transitions <- expand.grid(states, states, stringsAsFactors = FALSE)
  transitions <- transitions[transitions[,1] != transitions[,2],]
  for(i in 1:nrow(transitions)){
    trans <- node_df$left == transitions[i, 1] & node_df$right == transitions[i, 2]
    node_df <- cbind(node_df, trans)
  }
  colnames(node_df)[5:(4+nrow(transitions))] <- str_c(transitions[,1], "_", transitions[,2])
  node_df$transition <- NA
  
  for(i in 1:nrow(node_df)){
    trans_v <- unlist(node_df[i,5:(nrow(transitions)+4)])
    if(any(trans_v)){
      if(sum(trans_v) > 1){
        cat("more than 1 transition!")
      }else{
        node_df$transition[i] <- names(trans_v)[trans_v]
      }
    }else{
      node_df$transition[i] <- "No_change"
    }
  }
  
  # assign ages to nodes
  tree2 <- tree[!names(tree) %in% "node.label"]  # delete node labels to get node indices for the branching times
  class(tree2) <- "phylo"
  
  ages <- branching.times(tree2)
  node_df$time_left <- ages[match(node_df$V1, names(ages))]
  node_df$time_right <- ages[match(node_df$V2, names(ages))]
  
  node_df$time_left[is.na(node_df$time_left)] <- 0
  node_df$time_right[is.na(node_df$time_right)] <- 0
  
  
  # create time intervals with sliding window
  max_time <- max(c(node_df$time_left, node_df$time_right))
  
  i <- 0
  intervals <- list()
  count <- 1
  window <- step
  while(i <= max_time){
    interval_temp <- seq(from = i, to = max_time, by = window_width)
    intervals[[count]] <- interval_temp
    i <- i + window
    count <- count + 1
  }
  
  idx <- sapply(intervals, function(x) length(x) < 2)
  intervals <- intervals[!idx]
  
  # calculate mean number of transitions in each of the time interval
  # I consider that a transition can happen any time along a branch
  trans_df_all <- as.data.frame(matrix(ncol = 4))
  for(i in 1:length(intervals)){
    interval_v <- intervals[[i]]
    trans_df <- as.data.frame(matrix(nrow = length(interval_v) - 1, ncol = 4))
    for (j in 2:length(interval_v)) {
      x <- j - 1
      start <- interval_v[x]
      end <- interval_v[j]
      temp_df <- node_df[start >= node_df$time_right & start <= node_df$time_left,] 
      temp_df <- rbind(temp_df, node_df[end >= node_df$time_right & end <= node_df$time_left,] )
      temp_df <- rbind(temp_df, node_df[node_df$time_right >= start & node_df$time_right <= end,] )
      n <- nrow(temp_df)
      transition <- sum(!temp_df$transition %in% "No_change")
      mean_trans <- transition / n
      trans_df[x, 1] <- end
      trans_df[x, 2] <- n
      trans_df[x, 3] <- transition
      trans_df[x, 4] <- mean_trans
    }
    trans_df_all <- rbind(trans_df_all, trans_df)
  }
  trans_df_all <- na.omit(trans_df_all)
  colnames(trans_df_all) <- c("time", "n", "transition", "mean_trans")
  if(!is.null(trait)){
    trans_df_all$trait <- trait
  }
  return(as.data.frame(trans_df_all))
}

corhmm_to_ace <- function(corhmm){
  tree <- corhmm$phy
  nodes_df <- tree$edge
  states_data <- data.frame(nodes = min(nodes_df) : max(nodes_df),
                            sates = c(corhmm$tip.states, corhmm$states),
                            stringsAsFactors = FALSE)
  states_data_n <- nrow(states_data)
  states <- unique(states_data[,2])
  states <- states[order(states)]
  for(i in 1:length(states)){
    states_data <- cbind(states_data, rep(NA, states_data_n))
  }
  colnames(states_data) <- c(colnames(states_data)[1:2], states)
  for(i in 1:length(states)){
    state <- states[i]
    states_data[,(i+2)] <- unlist(lapply(states_data[,2], function(x) if(x %in% state){return(1)}else{0}))
  }
  states_data <- states_data[,-c(1,2)]
  return(as.matrix(states_data))
}


Calculate_cladogenic <- function(tree, asr, occurrence = TRUE, transition_n = FALSE, type = c("all", "tips", "inner_nodes")){
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
  switch(type,
         all = {node_df <- as.data.frame(tree$edge)},
         tips = {node_df <- as.data.frame(tree$edge[1:length(tree$tip.label),])},
         inner_nodes = {node_df <- as.data.frame(tree$edge[(length(tree$tip.label) + 1):nrow(tree$edge),])})
  node_df$left <- rep(NA, nrow(node_df))
  node_df$right <- rep(NA, nrow(node_df))
  node_df$left <- asr$states[match(node_df$V1, asr$node)]
  node_df$right <- asr$states[match(node_df$V2, asr$node)]
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

Calculate_section <- function(corHMM = NULL, min_time = 145, max_time = 155, trait = NULL) {
  
  ace_ML <- corhmm_to_ace(corHMM)
  asr <- as.data.frame(ace_ML)
  asr_nstates <- ncol(asr)
  states <- 1:asr_nstates
  asr_nrow <- nrow(asr)
  asr$states <- rep(NA, asr_nrow)
  for(i in 1:asr_nstates){
    asr$states[asr[,i] == 1] <- i
  }
  asr$nodes <- 1:asr_nrow
  
  tree <- corHMM$phy
  node_df <- as.data.frame(tree$edge)
  node_df$left <- rep(NA, nrow(node_df))
  node_df$right <- rep(NA, nrow(node_df))
  node_df$left <- asr$states[match(node_df$V1, asr$node)]
  node_df$right <- asr$states[match(node_df$V2, asr$node)]
  transitions <- expand.grid(states, states, stringsAsFactors = FALSE)
  transitions <- transitions[transitions[,1] != transitions[,2],]
  for(i in 1:nrow(transitions)){
    trans <- node_df$left == transitions[i, 1] & node_df$right == transitions[i, 2]
    node_df <- cbind(node_df, trans)
  }
  colnames(node_df)[5:(4+nrow(transitions))] <- str_c(transitions[,1], "_", transitions[,2])
  node_df$transition <- NA
  
  for(i in 1:nrow(node_df)){
    trans_v <- unlist(node_df[i,5:(nrow(transitions)+4)])
    if(any(trans_v)){
      if(sum(trans_v) > 1){
        cat("more than 1 transition!")
      }else{
        node_df$transition[i] <- names(trans_v)[trans_v]
      }
    }else{
      node_df$transition[i] <- "No_change"
    }
  }
  
  # assign ages to nodes
  tree2 <- tree[!names(tree) %in% "node.label"]  # delete node labels to get node indices for the branching times
  class(tree2) <- "phylo"
  
  ages <- branching.times(tree2)
  node_df$time_left <- ages[match(node_df$V1, names(ages))]
  node_df$time_right <- ages[match(node_df$V2, names(ages))]
  
  node_df$time_left[is.na(node_df$time_left)] <- 0
  node_df$time_right[is.na(node_df$time_right)] <- 0
  
  
  interval <- c(min_time, max_time)
  
  # I consider that a transition can happen any time along a branch
  
  trans_df <- as.data.frame(matrix(nrow = length(interval) - 1, ncol = 4))
  start <- interval[1]
  end <- interval[2]
  temp_df <- node_df[start >= node_df$time_right & start <= node_df$time_left,] 
  temp_df <- rbind(temp_df, node_df[end >= node_df$time_right & end <= node_df$time_left,] )
  temp_df <- rbind(temp_df, node_df[node_df$time_right >= start & node_df$time_right <= end,] )
  n <- nrow(temp_df)
  transition <- sum(!temp_df$transition %in% "No_change")
  mean_trans <- transition / n
  trans_df[1, 1] <- end
  trans_df[1, 2] <- n
  trans_df[1, 3] <- transition
  trans_df[1, 4] <- mean_trans
  
  trans_df_all <- na.omit(trans_df_all)
  colnames(trans_df) <- c("Time", "n_branches", "n_transitions", "mean_transitions")
  if(!is.null(trait)){
    trans_df$trait <- trait
  }
  return(as.data.frame(trans_df))
}
