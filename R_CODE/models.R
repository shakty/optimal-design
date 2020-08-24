##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###   Maximizing Information Gain          ###
###   code by: B. Klein and S.Balietti     ###
###                                        ###
##############################################
## Four models that generate datasets based on different principles
## For use in the Parameter-Sampled GPUCB-PE

model_1n <- function(pi_per, A, epsilon, alpha, action, playerdata) {
    ## See El-Gamal & Palfrey (1995) for more details
    pi_per_hat <- A/(2+A)
    epsil1_Case_A <- (2*pi_per_hat) * (1-pi_per)/(pi_per_hat + pi_per - 2*pi_per_hat*pi_per)
    epsil2_Case_A <- 2/A

    epsilon_conditionA <- epsilon <= epsil1_Case_A & epsilon <= epsil2_Case_A

    et <- epsilon
    etm <- 1 - epsilon
    et2 <- epsilon/2
    if (pi_per > pi_per_hat & epsilon_conditionA) {
        if (action=='a') {
            p_a <- A*(1-pi_per)/(2*pi_per) -
                ((A*pi_per+2*pi_per - A)*(epsilon/2))/((2*pi_per-2*pi_per*epsilon))
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action=='q') {
            q <- 1/(1-epsilon)*((A-1)/A - epsilon/2)
            q <- etm * q + et2
            return(q)
        }
    }
    epsilon_conditionB <- epsilon > epsil1_Case_A & epsilon <= epsil2_Case_A
    if (pi_per > pi_per_hat & epsilon_conditionB) {
        if (action=='a') {
            p_a <- 0
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action=='q') {
            q <- 1
            q <- etm * q + et2
            return(q)
        }
    }
    epsil1_Case_C <- ((2*pi_per) * (1-pi_per_hat))/(pi_per_hat + pi_per - 2*pi_per_hat*pi_per)
    if (pi_per <= pi_per_hat & epsilon <= epsil1_Case_C) {
        if (action=='a') {
            p_a <- 1
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b') {
            p_b <- (2*pi_per) / (A*(1-pi_per)) + (((A+2)*pi_per-A)*epsilon/2) / ((1-epsilon)*A*(1-pi_per))
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action=='q') {
            q <- 1/2
            q <- etm * q + et2
            return(q)
        }
    }
    epsil1_Case_D <- 2*pi_per * (1-pi_per_hat) / (pi_per_hat + pi_per - 2*pi_per_hat * pi_per)
    if (pi_per <= pi_per_hat & epsilon > epsil1_Case_D) {
        if (action=='a') {
            p_a <- 1
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b') {
            p_b <- 0
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action=='q') {
            q <- 0
            q <- etm * q + et2
            return(q)
        }
    }
    if (pi_per > pi_per_hat) {
        if (action=='a') {
            p_a <- 1
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action=='q') {
            q <- 1
            q <- etm * q + et2
            return(q)
        }
    }
}

model_2n <- function(pi_per, A, epsilon, alpha, action, playerdata) {
    L <- 2*pi_per
    R <- A*(1-pi_per)

    et <- epsilon
    etm <- 1 - epsilon
    et2 <- epsilon/2

    if (L > R) {
        if (action == 'a') {
            x <- epsilon*A/2
            if (x > 1) {
                p_a <- 1
            } else if (x == 1) {
                p_a <- .5
            } else if (x < 1) {
                p_a <- 0
            }
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action == 'b') {
            x <- 2*(1-(epsilon/2))
            if (x > 1) {
                p_b <- 1
            } else if (x == 1) {
                p_b <- .5
            } else if (x < 1) {
                p_b <- 0
            }
            p_b <- etm * p_b + et2
            return(p_b)
        }
        if (action == 'q') {
            q <- 1
            q <- etm * q + et2
            return(q)
        }
    }
    if (L < R) {
        if (action=='a') {
            x <- (1-(epsilon/2))*A
            if (x > 1) {
                p_a <- 1
            } else if (x == 1) {
                p_a <- .5
            } else if (x < 1) {
                p_a <- 0
            }
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action=='b' || action == 'q') {
            q <- 0
            q <- etm * q + et2
            return(q)
        }
    }
    if (L == R) {
        if (action=="a") {
            p_a <- 1
            p_a <- etm * p_a + et2
            return(p_a)
        }
        if (action == "b" || action == "q") {
            q <- 1/2
            q <- etm * q + et2
            return(q)
        }
    }
}

model_3n <- function(pi_per, A, epsilon, alpha, action, playerdata, emp = 0.5) {

    et <- epsilon
    etm <- 1 - epsilon
    et2 <- epsilon/2

    empa <- empb <- emp <- 0.5
    n_left_go  <- str_count(playerdata, "L")
    n_go       <- str_count(playerdata, "G")
    n_gamea    <- str_count(playerdata, "a")
    n_gameb    <- str_count(playerdata, "b")
    n_go_gamea <- str_count(playerdata, "aG")
    n_go_gameb <- str_count(playerdata, "bG")
    emp_t  <- (emp+n_left_go) / (n_go+1)
    emp_at <- (empa+n_go_gamea) / (n_gamea+1)
    emp_bt <- (empb+n_go_gameb) / (n_gameb+1)
    if (action == "a") {
        l = (1-emp_t) * A
        if      (l > 1) { p_a <- 1   }
        else if (l ==1) { p_a <- 0.5 }
        else if (l < 1) { p_a <- 0   }
        p_a <- etm * p_a + et2
        return(p_a)
    }
    if (action == 'b') {
        l = emp_t * 2
        if      (l > 1) { p_b <- 1 }
        else if (l ==1) { p_b <- 0.5 }
        else if (l < 1) { p_b <- 0 }
        p_b <- etm * p_b + et2
       	return(p_b)
    }
    if (action == 'q') {
        if (pi_per == 0) {
            emp_pi_t = 0
        } else { emp_pi_t <- (emp_at*pi_per) / (emp_at*pi_per + emp_bt*pi_per) }
        l = 2*emp_pi_t
        r = (1-emp_pi_t)*A
        if      ( l > r ) { q <- 1 }
        else if ( l ==r ) { q <- 0.5 }
        else if ( l < r ) { q <- 0 }
        q <- etm * q + et2
      	return(q)
    }
}

model_4n <- function(pi_per, A, epsilon, alpha, action, playerdata) {

    et <- epsilon
    etm <- 1 - epsilon
    et2 <- epsilon/2
    payoff <- data.frame("aSL"=c(1,1), "aSR"=c(1,1),
                         "bSL"=c(1,1), "bSR"=c(1,1),
                         "aGL"=c(0,2), "aGR"=as.numeric(c(A,0)),
                         "bGL"=c(2,0), "bGR"=as.numeric(c(0,A)))

    s <- pi_per * 10
    if (s == 0) s = 0.0000000000000001

    propA0 <- c(1,1); propB0 <- c(1,1); propQ0 <- c(1,1)
    refPayoff <- 0
    minPropensity = 0.00001

    propA <- propA0[1] * s; propB <- propB0[1] * s; propQ <- propQ0[1] * s
    propNotA <- propA0[2] * s; propNotB <- propB0[2] * s; propNotQ <- propQ0[2] * s
    pa <- propA / (propA + propNotA); pb <- propB / (propB + propNotB); pq <- propQ / (propQ + propNotQ)
    pNotA <- 1 - pa; pNotB <- 1 - pb; pNotQ <- 1 - pq

    ## Round 1
    if( nchar(playerdata) == 0 ) {
        if (action == "a") {
            pa <- etm * pa + et2
            return(pa)
        } else if (action == "b") {
            pb <- etm * pb + et2
            return(pb)
        } else {
            pq <- etm * pq + et2
            return(pq)
        }
    }

    playerNum <- 1
    if( action=='q' ) playerNum <- 2
    history <- str_split(playerdata, "-")[[1]]
    history <- history[history!=""]
    if(nchar(history[1])==4){
        history <- substr(history, 2, 4)
    }
    tmp <- c()
    for(hist in history){
        tmp <- c(tmp, payoff[playerNum, hist])
    }

    ## Round 2
    if( length(history) == 1 ) {
        Rx1 <- tmp[1] - refPayoff
        update <- Rx1 * alpha

        if (action == "a") {
            if (substr(history[1], 2, 2) == "S") {
                propNotA <- propNotA + update
                if (propNotA <= 0) propNotA <- minPropensity
            } else {
                propA <- propA + update
                if (propA <= 0) propA <- minPropensity
            }
            pa <- propA / (propA + propNotA)
            pa <- etm * pa + et2
            return(pa)
        } else if (action == "b") {
            if (substr(history[1], 2, 2) == "S") {
                propNotB <- propNotB + update
                if (propNotB <= 0) propNotB <- minPropensity
            } else {
                propB <- propB + update
                if (propB <= 0) propB <- minPropensity
            }
            pb <- propB / (propB + propNotB)
            pb <- etm * pb + et2
            return(pb)
        } else {
            if (substr(history[1], 3, 3) == "R") {
                propNotQ <- propNotQ + update
                if (propNotQ <= 0) propNotQ <- minPropensity
            } else {
                propQ <- propQ + update
                if (propQ <= 0) propQ <- minPropensity
            }
            pq <- propQ / (propQ + propNotQ)
            pq <- etm * pq + et2
            return(pq)
        }
    }

    ## Round 3
    if( length(history) == 2 ) {
        Rx2 <- sum(tmp) - refPayoff
        update <- Rx2 * alpha

        if (action == "a") {
            if (substr(history[2], 2, 2) == "S") {
                propNotA <- propNotA + update
                if (propNotA <= 0) propNotA <- minPropensity
            } else {
                propA <- propA + update
                if (propA <= 0) propA <- minPropensity
            }
            pa <- propA / (propA + propNotA)
            pa <- etm * pa + et2
            return(pa)
        } else if (action == "b") {
            if (substr(history[2], 2, 2) == "S") {
                propNotB <- propNotB + update
                if (propNotB <= 0) propNotB <- minPropensity
            } else {
                propB <- propB + update
                if (propB <= 0) propB <- minPropensity
            }
            pb <- propB / (propB + propNotB)
            pb <- etm * pb + et2
            return(pb)
        } else {
            if (substr(history[2], 3, 3) == "R") {
                propNotQ <- propNotQ + update
                if (propNotQ <= 0) propNotQ <- minPropensity
            } else {
                propQ <- propQ + update
                if (propQ <= 0) propQ <- minPropensity
            }
            pq <- propQ / (propQ + propNotQ)
            pq <- etm * pq + et2
            return(pq)
        }
    }
}

