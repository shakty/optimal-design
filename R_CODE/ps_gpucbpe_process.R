##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###       Maximizing Information Gain      ###
###        code by: Brennan Klein          ###
###                                        ###
##############################################

## Initializing the landscape and search process.
randomRows <- function(datagrid,n_init=25,sobol="Sobol"){
    ## During initialization, this function is called to generate
    ## the initial samples of points on the information surface.
    if (sobol=="Sobol") {
        return(getSobol(datagrid, n_init))
    }
    samp <- datagrid[sample(nrow(datagrid), n_init),]
    if (any(duplicated(samp)==T)){
        samp <- datagrid[sample(nrow(datagrid), n_init),]
    }
    return(samp)
}

getSobol <- function(datagrid, n_init) {
    ## Generates a random Sobol Sequence...
    ## this is a quasi-random sequence, with a random seed.

    df <- rbind(c(0,0), sobol(n_init, 2, scrambling=3, seed=sample(c(1:1000), 1)))
    colnames(df) <-c('PI', 'A')
    df <- data.frame(df)
    sampy <- min(datagrid$A)  + df$A * (max(datagrid$A) - min(datagrid$A))
    sampx <- min(datagrid$PI) + df$PI* (max(datagrid$PI)- min(datagrid$PI))
    out <- data.frame("PI"=sampx,"A"=sampy)
    out <- out[-1,]
    rownames(out) <- 1:nrow(out)
    return(out)
}

initialize_process <- function(posteriors, n_init, expParam1_range,
                               expParam2_range, n_dim, sampling_method,
                               n_players, n_rounds, n_samples, sobol, sampling,
                               asymmetric, L, models_used, current_top_model=1){

    ## First: Call this function to get the initial dataframe and the
    ## corresponding information numbers for each point in the surface.

    datagrid <- expand.grid(A=seq(expParam1_range[1], expParam1_range[2],
                                  length.out=n_dim),
                            PI=seq(expParam2_range[1],
                                   expParam2_range[2], length.out=n_dim))
    history <- randomRows(datagrid, n_init, sobol)
    history$I <- 0
    history$Stop <- 2
    history$Round <- c(1:dim(history)[1])


    str <- paste0("Initializing %02i %s seeds with %05s sampled histories per ",
                  "search ... Current Round is: %02i ... Current Time is: %s")

    for (row in 1:(nrow(history)-L+1)) {
        a <- history[row,]$A
        p <- history[row,]$PI

        history[row,]$I <- getInfoNumber(a, p, sampling_method, posteriors,
                                         n_samples, asymmetric, sampling,
                                         models_used, n_players, n_rounds,
                                         current_top_model)
        print(sprintf(str, n_init, sobol, n_samples, row, Sys.time()))
    }

    GPredict_prev <- run_gp(history[1:(nrow(history)-L+1),], datagrid)
    for (row in (nrow(history)-L+2):nrow(history)) {
        a <- history[row,]$A
        p <- history[row,]$PI
        history[row,]$I <- getInfoNumber(a, p, sampling_method, posteriors,
                                         n_samples, asymmetric, sampling,
                                         models_used, n_players, n_rounds,
                                         current_top_model)
        GPredict_curr <- run_gp(history[1:row,], datagrid)
        history[row,]$Stop <- stopping_criteria(GPredict_prev, GPredict_curr)
        GPredict_prev <- GPredict_curr

        print(sprintf(str, n_init, sobol, n_samples, row, Sys.time()))
    }
    return(history)
}


###################################
### Gaussian Process functions. ###
###################################
                                        # Gaussian process

run_algorithm_once <- function(history, datagrid, posteriors, rho0, failsafe, L, n_samples, k, asymmetric,
                               sampling, search, n_init, sampling_method, sobol, seed, models_used, current_top_model=1){
                                        # Iterates through the Gaussian Process functions, searching more and more spaces.
    current_round <- dim(history)[1]
    if(search=="GPUCBPE"){
        history <- add_gpucb_pe_point(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, sampling_method, models_used, current_top_model)
        current_round <- dim(history)[1]
        while(current_round < failsafe){
            history <- add_gpucb_pe_point(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, sampling_method, models_used, current_top_model)
            current_round <- dim(history)[1]
            stop <- all(tail(history,L)$Stop > (1-rho0))
            if(stop){
                print(sprintf("Current Round is: %03i (GPUCBPE) ... Stop is: %04f ... Currently stopping at: %s", current_round, round(tail(history,1)$Stop, 5), Sys.time()))
                history$sampling_method <- sampling_method; history$Search <- search; history$n_init <- n_init
                history$n_samples <- n_samples; history$Sobol <- sobol; history$Sampling <- sampling
                history$csvname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="CSV")
                history$pngname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="FIG")
                gp_seq <- c("UCB",rep("PE", k))
                times <- trunc(dim(history)[1]/length(gp_seq))
                plus <- dim(history)[1] %% length(gp_seq)
                history$SearchLabel <- c(rep(gp_seq, times),head(gp_seq,plus))
                return(history)
            }
        }
        print(sprintf("Current Round is: %03i (GPUCBPE) ... Stop is: %04f ... Currently stopping at: %s", current_round, round(tail(history,1)$Stop, 5), Sys.time()))
        history$sampling_method <- sampling_method; history$Search <- search; history$n_init <- n_init
        history$n_samples <- n_samples; history$Sobol <- sobol; history$Sampling <- sampling
        history$csvname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="CSV")
        history$pngname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="FIG")
        gp_seq <- c("UCB",rep("PE", k))
        times <- trunc(dim(history)[1]/length(gp_seq))
        plus <- dim(history)[1] %% length(gp_seq)
        history$SearchLabel <- c(rep(gp_seq, times),head(gp_seq,plus))
        return(history)
    }
    if(search=="Random"){
        current_round <- dim(history)[1]
        history <- add_random_point(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, search, n_init, sampling_method, models_used, current_top_model)
        current_round <- dim(history)[1]
        while(current_round < failsafe){
            history <- add_random_point(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, search, n_init, sampling_method, models_used, current_top_model)
            current_round <- dim(history)[1]
            stop <- all(tail(history,L)$Stop > (1-rho0))
            if(stop){
                print(sprintf("Current Round is: %03i (Random) ... Stop is: %04f ... Currently stopping at: %s", current_round, round(tail(history,1)$Stop, 5), Sys.time()))
                history$sampling_method <- sampling_method; history$Search <- search; history$n_init <- n_init
                history$n_samples <- n_samples; history$Sobol <- sobol; history$Sampling <- sampling
                history$csvname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="CSV")
                history$pngname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="FIG")
                history$SearchLabel <- "Rand"
                return(history)
            }
        }
        print(sprintf("Current Round is: %03i (Random) ... Stop is: %04f ... Currently stopping at: %s", current_round, round(tail(history,1)$Stop, 5), Sys.time()))
        history$sampling_method <- sampling_method; history$Search <- search; history$n_init <- n_init
        history$n_samples <- n_samples; history$Sobol <- sobol; history$Sampling <- sampling
        history$csvname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="CSV")
        history$pngname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="FIG")
        history$SearchLabel <- "Rand"
        return(history)
    }
    return(history)
}

run_algorithm_grid <- function(datagrid, posteriors, n_samples, asymmetric, sampling, sampling_method, seed, models_used=c(1,2,3,4), current_top_model=1){
    k <- 0; n_init <- 0; search <- "Grid"; sobol <- "Grid"
    history <- datagrid
    history$I <- 0; history$Round <- c(1:dim(datagrid)[1]); history$Stop <- 2
    for(row in history$Round){
        A <- history[row,]$A
        PI <- history[row,]$PI
        history[row,]$I <- getInfoNumber(A, PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model)
        print(sprintf("Current Round is: %03i (Grid) ... Just dug at: (%04f, %04f) ... Current Time is: %s", row, round(history[row,]$PI, 5), round(history[row,]$A, 5), Sys.time()))
    }
    history$sampling_method <- sampling_method; history$Search <- search; history$n_init <- n_init
    history$n_samples <- n_samples; history$Sobol <- sobol; history$Sampling <- sampling
    history$csvname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="CSV")
    history$pngname <- make_filename(asymmetric, sampling, n_init, sampling_method, n_samples, search, seed, sobol, models_used, current_top_model, dtype="FIG")
    history$SearchLabel <- "Grid"
    return(history)
}

add_random_point <- function(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, search, n_init, sampling_method, current_top_model=1){
                                        # Adds one new random point to history.
    GPredict_prev <- run_gp(history, datagrid)
    new_point <- data.frame(PI=runif(1,.2,.8), A=runif(1,2,6))
    new_point$I <- getInfoNumber(new_point$A, new_point$PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model)
    new_point$Round <- current_round+1; new_point$Stop <- 2;
    history_curr <- rbind(history, new_point)
    GPredict_curr <- run_gp(history_curr, datagrid)
    current_round <- dim(history_curr)[1]
    history_curr$Stop[current_round] <- stopping_criteria(GPredict_prev, GPredict_curr)
    history 	  <- history_curr
    GPredict_prev <- GPredict_curr
    print(sprintf("Current Round is: %03i (Random) ... Just dug at: (%04f, %04f) ... Stop is: %04f ... Current Time is: %s",
                  current_round, round(history[current_round-1,]$PI, 5), round(history[current_round-1,]$A, 5), round(tail(history,1)$Stop, 5), Sys.time()))
    return(history)
}

add_gpucb_pe_point <- function(history, datagrid, current_round, n_samples, posteriors, k, L, asymmetric, sampling, rho0, sampling_method, models_used, current_top_model=1){
                                        # Does two things: First, a point is added that is the max(y_hat)+UCB. Then it adds k PE points.
### UCB step
    GPredict_prev <- run_gp(history, datagrid)
    ucb <- add_ucb_point(GPredict_prev, history, datagrid, current_round, n_samples, posteriors, asymmetric, sampling, sampling_method, models_used, current_top_model)
    history 	  <- ucb[[1]]
    GPredict_prev <- ucb[[2]]
    current_round <- ucb[[3]]
    print(sprintf("Current Round is: %03i (UCB) ... Just dug at: (%04f, %04f) ... Stop is: %04f ... Current Time is: %s",
                  current_round, round(history[current_round,]$PI, 5), round(history[current_round,]$A, 5), round(tail(history,1)$Stop, 5), Sys.time()))
    stop <- all(tail(history,L)$Stop > (1-rho0))
    if(stop){
        return(history)
    }

### PE step
    for(i in 1:k){
        pe <- add_pe_point(GPredict_prev, history, datagrid, current_round, n_samples, posteriors, asymmetric, sampling, sampling_method, models_used, current_top_model)
        history 	  <- pe[[1]]
        GPredict_prev <- pe[[2]]
        current_round <- pe[[3]]
        kind <- pe[[4]]
        print(sprintf("Current Round is: %03i %s... Just dug at: (%04f, %04f) ... Stop is: %04f ... Current Time is: %s",
                      current_round, kind, round(history[current_round,]$PI, 5), round(history[current_round,]$A, 5), round(tail(history,1)$Stop, 5), Sys.time()))
        stop <- all(tail(history,L)$Stop > (1-rho0))
        if(stop){
            return(history)
        }
    }
    return(history)
}

jitter_points <- function(coords, datagrid, history){
    jitterA  <- (unique(datagrid$A)[2]  - unique(datagrid$A )[1])/2
    jitterPI <- (unique(datagrid$PI)[2] - unique(datagrid$PI)[1])/2
    if(typeof(coords)!="list"){
        newA  <- coords[1] + runif(1, -jitterA,  jitterA)
        newPI <- coords[2] + runif(1, -jitterPI, jitterPI)
        coords <- c(newA, newPI)
    } else {
        newA  <- coords$A + runif(1, -jitterA,  jitterA)
        newPI <- coords$PI + runif(1, -jitterPI, jitterPI)
        for(round in c(1:dim(history)[1])){
            if(newPI[round] > max(datagrid$PI)){
                newPI[round] <- max(datagrid$PI)
            } else if(newPI[round] < min(datagrid$PI)){
                newPI[round] <- min(datagrid$PI)
            }
            if(newA[round] > max(datagrid$A)){
                newA[round] <- max(datagrid$A)
            } else if(newA[round] < min(datagrid$A)){
                newA[round] <- min(datagrid$A)
            }
        }
        coords <- data.frame("A"=newA, "PI"=newPI)
    }
    return(coords)
}

add_ucb_point <- function(GPredict_prev, history, datagrid, current_round, n_samples, posteriors, asymmetric, sampling, sampling_method, models_used, current_top_model=1){
                                        # Adds one new UCB point to history.
    y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
    UCB_point <- GPredict_prev$complete_data[which.max(y_plus_sigma),]
    coords <- as.numeric(UCB_point[1:2])
    new_point <- data.frame(PI=coords[2], A=coords[1])
    new_point$I <- getInfoNumber(new_point$A, new_point$PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model)
    new_point$Round <- current_round+1; new_point$Stop <- 2
    history_curr <- rbind(history, new_point)
    GPredict_curr <- run_gp(history_curr, datagrid)
    current_round <- dim(history_curr)[1]
    history_curr$Stop[current_round] <- stopping_criteria(GPredict_prev, GPredict_curr)
    return(list(history_curr, GPredict_curr, current_round))
}

add_pe_point <- function(GPredict_prev, history, datagrid, current_round, n_samples, posteriors, asymmetric, sampling, sampling_method, models_used, current_top_model=1){
                                        # Adds one new PE point to history.
    y_minus_sigma <- GPredict_prev$Y_hat - 2*(GPredict_prev$MSE)**.5
    y_plus_sigma  <- GPredict_prev$Y_hat + 2*(GPredict_prev$MSE)**.5
    region <- which(y_plus_sigma > max(y_minus_sigma))
    if(length(region) > 1){
        space_to_sample <- GPredict_prev$complete_data[region,]
        PE_point <- space_to_sample[which.max(space_to_sample[,4]),]
        coords <- as.numeric(PE_point[1:2])
        new_point <- data.frame(PI=coords[2], A=coords[1])
        new_point$I <- getInfoNumber(new_point$A, new_point$PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model)
        new_point$Round <- current_round+1; new_point$Stop <- 2
        history_curr <- rbind(history, new_point)
        GPredict_curr <- run_gp(history_curr, datagrid)
        current_round <- dim(history_curr)[1]
        history_curr$Stop[current_round] <- stopping_criteria(GPredict_prev, GPredict_curr)
        kind <- "(PE)  "

    } else {
        PE_point <- GPredict_prev$complete_data[which.max(GPredict_prev$MSE),]
        coords <- as.numeric(PE_point[1:2])
        new_point <- data.frame(PI=coords[2], A=coords[1])
        new_point$I <- getInfoNumber(new_point$A, new_point$PI, sampling_method, posteriors, n_samples, asymmetric, sampling, models_used, n_players, n_rounds, current_top_model)
        new_point$Round <- current_round+1; new_point$Stop <- 2
        history_curr <- rbind(history, new_point)
        GPredict_curr <- run_gp(history_curr, datagrid)
        current_round <- dim(history_curr)[1]
        history_curr$Stop[current_round] <- stopping_criteria(GPredict_prev, GPredict_curr)
        kind <- "(PER) "
    }
    return(list(history_curr, GPredict_curr, current_round, kind))
}

run_gp <- function(history, datagrid){
                                        # Takes the current history and creates a landscape using predict and GP_fit.
    temp <- history
                                        # temp$A <- round(temp$A, 3); temp$PI <- round(temp$PI, 3)
    temp$A <- round(temp$A, 6); temp$PI <- round(temp$PI, 6)
    history_curr <- suppressWarnings(aggregate(temp, by=list(temp$A, temp$PI), FUN=mean))
    anorm <- (history_curr$A  - min(datagrid$A))  / (max(datagrid$A)  - min(datagrid$A))
    pnorm <- (history_curr$PI - min(datagrid$PI)) / (max(datagrid$PI) - min(datagrid$PI))
    evidence <- data.frame("A"= anorm, "PI"=pnorm, "I"= history_curr$I)
    gp_model <- GP_fit(evidence[,c("A", "PI")], evidence$I)
    gp_model$X <- history_curr[,c("A","PI")]
    GPredict <- predict(gp_model, datagrid)
    return(GPredict)
}

stopping_criteria <- function(GPredict_prev, GPredict_curr){
                                        # Stop when the procedure ceases to learn about the landscape, when comparing
                                        # the global changes in mu between two successive iterations.
    prev_df <- data.frame(GPredict_prev$complete_data)
    curr_df <- data.frame(GPredict_curr$complete_data)
    prev_df$num <- c(1:dim(prev_df)[1])
    curr_df$num <- c(1:dim(curr_df)[1])
    pi_prev_sorted <- prev_df[order(prev_df$Y_hat),]
    pi_curr_sorted <- curr_df[order(curr_df$Y_hat),]
    pi_prev_rank <- pi_prev_sorted$num
    pi_curr_rank <- pi_curr_sorted$num
    numerator <- discounted_rank_dissimilarity(pi_curr_rank, pi_prev_rank)
    denominator <- get_max_distance(length(pi_curr_rank))
    rhoXv <- 1 - numerator/denominator
    return(rhoXv)
}

get_rank_distances <- function(pi_t1, pi_t0){
                                        # Returns the distance between two ordered lists.
    return(c(pi_t1 - pi_t0)**2)
}

discounted_rank_dissimilarity <- function(pi_t1, pi_t0){
                                        # Calculates the rank dissimilarity.
    numerator <- get_rank_distances(pi_t1, pi_t0)
    denominator <- pi_t1**2
    dists <- numerator/denominator
    d <- sum(dists)
    return(d)
}

get_max_distance <- function(nv){
                                        # Find the maximum distance.
    curr_max <- nv
    for(i in 1:nv){
        normal <- c(1:i); revers <- c(i:1)
        d <- discounted_rank_dissimilarity(normal,revers)
        if(d > curr_max){ curr_max <- d }
    }
    return(curr_max)
}
