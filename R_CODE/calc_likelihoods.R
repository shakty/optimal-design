##############################################
###                                        ###
###   Adaptive Design of Experiments for   ###
###       Maximizing Information Gain      ###
###        code by: Brennan Klein          ###
###                                        ###
##############################################
# Four models that generate likelihoods of each dataset
# For use in calculating likelihoods when spanning all datasets
#######################
### Likelihoods

## a is the prob that player1 chooses GO in a world that is A.
## b is the prob that player1 chooses GO in a world that is B.
## q is the prob that player2 chooses LEFT.
roles <- c('r','l')
n_players <- 2

getLikParamsForModel <- function(A, PI, thresh, model, batch, epsilons=NA, alphas=NA, deltas=NA,
                                 priorsEpsilons=NA, priorsAlphas=NA, priorsDeltas=NA) {
  	if (any(is.na(epsilons))) epsilons <- seq(0, 1, thresh)
  	lenEpsilon <- length(epsilons)
  	if (any(is.na(priorsEpsilons))) priorsEpsilons <- rep(1/lenEpsilon, lenEpsilon)
  	## Alphas.
  	if (any(is.na(alphas))) alphas <- seq(0, 1, thresh)
  	lenAlpha <- length(alphas)
  	if (any(is.na(priorsAlphas))) priorsAlphas <- rep(1/lenAlpha, lenAlpha)
  	## Deltas.
  	if (any(is.na(deltas))) deltas <- seq(0, .2, thresh)
  	lenDelta <- length(deltas)
  	if (any(is.na(priorsDeltas))) priorsDeltas <- rep(1/lenDelta, lenDelta)

	gameStats <- NA
	if (model == 1) {
    	lik <- innerFunMod1
  	} else if (model == 2) {
    	lik <- innerFunMod2
  	} else if (model == 3) {
    	gameStats <- list()
    for (p in 1:n_players) {            
      gameStats[[p]] <- getGameStats(batch[p,]) 
    }
    lik <- innerFunMod3
  	} else if (model == 4) {
    	lik <- innerFunModReinforcement
  	} else {
    	lik <- innerFunModRandom
  	}
  	out <- 0
	for (idxEpsilon in 1:lenEpsilon) {
		epsilon <- epsilons[idxEpsilon]
	    priorEpsilon <- priorsEpsilons[idxEpsilon]
	    for (idxAlpha in 1:lenAlpha) {
	    	alpha <- alphas[idxAlpha]
	    	priorAlpha <- priorsAlphas[idxAlpha]
	    	for (idxDelta in 1:lenDelta) {
	       		delta <- deltas[idxDelta]
	        	priorDelta <- priorsDeltas[idxDelta]          
		        pipers <- seq(max(0, PI - delta), min(1, PI + delta), thresh)
		    	for (piper in pipers) {
	          		res <- lik(piper, A, epsilon, alpha, batch, gameStats)
 		        	posterior <- res * priorEpsilon * priorAlpha * priorDelta
		        	out <- out + posterior
	            }
	        }     
	    }
	}
	return(out)
}

innerFunMod1 <- function(pi_per, A, e, alpha, batch, gameStats = NA) {
   
  ## Total likelihood.
  res <- 1
    
  for (round in 1:3) {
    # Tremble rate.
    et <- e * exp(- alpha * round)
    etm <- 1 - et 
    et2 <- et/2
    for (p in 1:n_players) {
      g <- batch[p,round]
      ## Red.
      if (roles[p] == 'r') {
        world <- substr(g, 1, 1)
        likAction <- model_1(pi_per, A, et, world)
        if (substr(g, 2, 2) == "S") {
            likAction <- 1 - likAction
        }
      } else {
        ## Blue.
        likAction <- model_1(pi_per, A, et, "q")
        if (substr(g, 3, 3) == "R") {
            likAction <- 1 - likAction
        }
      }
      likActionTremble <- etm * likAction + et2
      res <- res * likActionTremble
    }
  }
  return(res)
}

innerFunMod2 <- function(pi_per, A, e, alpha, batch, gameStats = NA) {
   
  ## Total likelihood.
  res <- 1
    
  for (round in 1:3) {
    # Tremble rate.
    et <- e * exp(- alpha * round)
    etm <- 1 - et 
    et2 <- et/2
    for (p in 1:n_players) {
      g <- batch[p,round]
      ## Red.
      if (roles[p] == 'r') {
        world <- substr(g, 1, 1)
        likAction <- model_2(pi_per, A, et, world)
        if (substr(g, 2, 2) == "S") {
            likAction <- 1 - likAction
        }
      } else {
        ## Blue.
        likAction <- model_2(pi_per, A, et, "q")
        if (substr(g, 3, 3) == "R") {
            likAction <- 1 - likAction
        }
      }
      likActionTremble <- etm * likAction + et2
      res <- res * likActionTremble
    }
  }
  return(res)
}

innerFunMod3 <- function(pi_per, A, e, alpha, batch, gameStats) {
 
  ## Total likelihood.
  res <- 1

  for (round in 1:3) {
    ## Tremble rate.
    et <- e * exp(- alpha * round)
    etm <- 1 - et 
    et2 <- et/2
    for (p in 1:n_players) {
      g <- batch[p, round]
      ## Cumulative history.
      gameStat <- gameStats[[p]][round,]
      ## Red.
      if (roles[p] == 'r') {
        world <- substr(g, 1, 1)
        likAction <- model_3(pi_per, A, et, world, gameStat, round)
        if (substr(g, 2, 2) == "S") {
          likAction <- 1 - likAction
        }
      } else {
        ## Blue.
        likAction <- model_3(pi_per, A, et, "q", gameStat, round)
        if (substr(g, 3, 3) == "R") {
          likAction <- 1 - likAction
        }
      }
      likActionTremble <- etm * likAction + et2
      res <- res * likActionTremble
    }
  }
  return(res)
}

innerFunModReinforcement <- function(pi_per, A, e, alpha, batch, gameStats = NA) {    
  s <- pi_per * 10
  if (s == 0) s = 0.0000000000000001
	               	   
  payoff <- list(aSL=c(1,1), aSR=c(1,1), aGL=c(0,2), aGR=as.numeric(c(A,0)), 
                 bSL=c(1,1), bSR=c(1,1), bGL=c(2,0), bGR=as.numeric(c(0,A)))
  res <- 1
  propA0 <- c(1,1); propB0 <- c(1,1); propQ0 <- c(1,1)
  refPayoff <- 0
  minPropensity = 0.00001
	
  propA <- propA0[1] * s; propB <- propB0[1] * s; propQ <- propQ0[1] * s
  propNotA <- propA0[2] * s; propNotB <- propB0[2] * s; propNotQ <- propQ0[2] * s
  pa <- propA / (propA + propNotA); pb <- propB / (propB + propNotB); pq <- propQ / (propQ + propNotQ)
  pNotA <- 1 - pa; pNotB <- 1 - pb; pNotQ <- 1 - pq
    
  for (round in 1:3) {
    ## Tremble rate.
    et <- e * exp(- alpha * round)
    etm <- 1 - et 
    et2 <- et/2
    for (p in 1:n_players) {
      ## Player history.    
      g <- batch[p, 1:round]
      ## Red.
      role <- roles[p]
      likAction <- model_4(pi_per, payoff, role, g, round, alpha,
                           propA0=propA0, propB0=propB0, propQ0=propQ0,
                           refPayoff=refPayoff)
      if (likAction <= 0) {
#            print(payoff); print(g); print(likAction)
            stop(paste('AAAAAAAh ', likAction, payoff, role, g, round, alpha, "\n"))
      }
      if (roles[p] == 'r') {
	    if (substr(g, 2, 2)[round] == "S") {
	        likAction <- 1 - likAction
	    }    
	  }                           
      if (roles[p] == 'r') {
	    if (substr(g, 2, 2)[round] == "S") {
	        likAction <- 1 - likAction
	    }    
	  }                           
      if (roles[p] == 'l') {
        if (substr(g, 3, 3)[round] == "R") {
          likAction <- 1 - likAction
        }
	  }                           
      
      likActionTremble <- etm * likAction + et2
      res <- res * likActionTremble
    }
  }
  return(res)
}

innerFunModRandom <- function(pi_per, A, epsilon, alpha, batch, gameStats = NA) {
    ## Probability of any action.
    likAction = 0.5
    return((likAction*likAction)^3)
}

nA_idx= 1
nB_idx = 2
nGoLeft_idx = 3
nGoRight_idx = 4
nAGo_idx = 5
nBGo_idx = 6
nAStop_idx = 7
nBStop_idx = 8
nGo_idx = 9
nStop_idx = 10
nStopRight_idx = 11
nStopLeft_idx = 12

# Returns a matrix cumulative stats for each round.
getGameStats <- function(games, toRound = NA) {
  ## Init.
  nA = 0
  nB = 0
  nGoLeft = 0
  nGoRight = 0
  nAGo = 0
  nBGo = 0
  nAStop = 0
  nBStop = 0
  nGo = 0
  nStop = 0
  nStopRight = 0
  nStopLeft = 0


  ## Plus 1 because round 1 has all zeros.
  nGames <- length(games) + 1

  ## If not specified we compute stats for all rounds.
  if (is.na(toRound)) toRound <- nGames

  # Init matrix.
  res <- matrix(nrow=toRound, ncol=12)
    
  for (r in 1:toRound) {
    if (r > 1) {
      ## Round 1 has no history.
      g <- games[r-1]
      world <- substr(g, 1, 1)
      redMove <- substr(g, 2, 2)
      ## World.
      if (world == 'a') nA = nA + 1
      else nB = nB + 1
      ## Red.
      if (redMove == 'G') {
        ## Go.
        nGo = nGo + 1
        if (world == 'a') nAGo = nAGo + 1
        else nBGo = nBGo + 1
        ## Blue.
        blueMove <- substr(g, 3, 3)
        if (blueMove == 'L') nGoLeft = nGoLeft + 1
        else  nGoRight = nGoRight + 1
      } else {
        ## Stop
        nStop = nStop + 1
        if (world == 'a') nAStop = nAStop + 1
        else nBStop = nBStop + 1
        ## Blue (supports decision after Stop or not).
        if (nchar(g) == 3) {
          blueMove <- substr(g, 3, 3)
          if (blueMove == 'L') nStopLeft = nStopLeft + 1
          else  nStopRight = nStopRight + 1
        }
      }
    }
    res[r, nA_idx] <- nA   
    res[r, nB_idx] <- nB
    res[r, nGoLeft_idx] <- nGoLeft
    res[r, nGoRight_idx] <- nGoRight
    res[r, nAGo_idx] <- nAGo
    res[r, nBGo_idx] <- nBGo
    res[r, nAStop_idx] <- nAStop
    res[r, nBStop_idx] <- nBStop
    res[r, nGo_idx] <- nGo
    res[r, nStop_idx] <- nStop
    res[r, nGo_idx] <- nGo
    res[r, nStopRight_idx] <- nStopRight
    res[r, nStopLeft_idx] <- nStopLeft
  }
  return(res)
}

model_1 <- function(pi_per, A, epsilon, action) {
    if (length(pi_per) > 1) {
        stop(paste("something wrong: pi_per more than one element",
                   paste(pi_per, collapse=", ")))
    }
	
    ## Player 1 and Player 2 both play the Bayes Nash equilibrium of the
    ## game, defined by (pi_per, A, epsilon)

    pi_per_hat <- A/(2+A)

    epsil1_Case_A <- (2*pi_per_hat) *
        (1-pi_per)/(pi_per_hat + pi_per - 2*pi_per_hat*pi_per)

    epsil2_Case_A <- 2/A
    epsilon_conditionA <- epsilon <= epsil1_Case_A & epsilon <= epsil2_Case_A
    
    if (pi_per > pi_per_hat & epsilon_conditionA) {
        
      if (action=='a') {
          p_a <- A*(1-pi_per)/(2*pi_per) -
            ((A*pi_per+2*pi_per - A)*(epsilon/2))/((2*pi_per-2*pi_per*epsilon))
            
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            return(p_b)
        }
        if (action=='q') {
            q <- 1/(1-epsilon)*((A-1)/A - epsilon/2)            
            return(q)
        }
    }
    
    epsilon_conditionB <- epsilon > epsil1_Case_A & epsilon <= epsil2_Case_A
    
    if (pi_per > pi_per_hat & epsilon_conditionB) {
        if (action=='a') {
            p_a <- 0
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            return(p_b)
        }
        if (action=='q') {
            q <- 1
            return(q)
        }
    }
    
    epsil1_Case_C <- ((2*pi_per) *
                      (1-pi_per_hat))/(pi_per_hat + pi_per - 2*pi_per_hat*pi_per)

    if (pi_per <= pi_per_hat & epsilon <= epsil1_Case_C) {
        if (action=='a') {
            p_a <- 1
            return(p_a)
        }
        if (action=='b') {
          p_b <- (2*pi_per) / (A*(1-pi_per)) + 
                (((A+2)*pi_per-A)*epsilon/2) / 
                ((1-epsilon)*A*(1-pi_per))                        
            return(p_b)
        }
        if (action=='q') {
            q <- 1/2
            return(q)
        }
    }
    
    epsil1_Case_D <- 2*pi_per * (1-pi_per_hat) / (pi_per_hat + pi_per - 2*pi_per_hat * pi_per)
    
    if (pi_per <= pi_per_hat & epsilon > epsil1_Case_D) {
        if (action=='a') {
            p_a <- 1
            return(p_a)
        }
        if (action=='b') {
            p_b <- 0
            return(p_b)
        }
        if (action=='q') {
            q <- 0
            return(q)
        }
    }
    
    if (pi_per > pi_per_hat) {
        if (action=='a') {
            p_a <- 1
            return(p_a)
        }
        if (action=='b') {
            p_b <- 1
            return(p_b)
        }
        if (action=='q') {
            q <- 1
            return(q)
        }
    }
}

model_2 <- function(pi_per, A, epsilon, action) {
    ## Player 2, after Player 1 plays 'go', does not update their pi_per,
    ## and this is assumed to be known by both players

    L <- 2*pi_per
    R <- A*(1-pi_per)
    
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
            return(p_b)
        }

        if (action == 'q') {
            q <- 1
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
            return(p_a)
        }
        if (action=='b' || action == 'q') {
            ## p_b <- q <- 0
            return(0)
        }
    }

    if (L == R) {
        if (action=="a") {
            p_a <- 1
            return(p_a)
        }
        if (action == "b" || action == "q") {
            return(1/2)
        }
    }

    stop("What??")
}


model_3 <- function(pi_per, A, epsilon, action, roundStat, round, emp = 0.5) {
    ## Individuals use fictucious play to construct beliefs about opponents play
    ## emp is the proportion of times that player1 thinks player2 is gonna go left
    ## empa and empb are the opposite, but wheter player1 chooses go under
    ## A or B state of the world.
    
    empa <- emp
    empb <- emp
    n_left_go <- roundStat[nGoLeft_idx]
    n_go <- roundStat[nGo_idx]

    n_gamea <- roundStat[nA_idx]
    n_gameb <- roundStat[nB_idx]
    n_go_gamea <- roundStat[nAGo_idx]
    n_go_gameb <- roundStat[nBGo_idx]

    emp_t <- (emp+n_left_go) / (n_go+1)
    emp_at <- (empa+n_go_gamea) / (n_gamea+1)
    emp_bt <- (empb+n_go_gameb) / (n_gameb+1)

    if (action == "a") {
        l = (1-emp_t)*A
        if (l > 1) {
            p_a <- 1
        }
        else if (l ==1) {
            p_a <- 0.5
        }
        else if (l < 1) {
            p_a <- 0
        }
        return(p_a)
    }

    if (action == 'b') {
        l = emp_t*2
        if (l > 1) {
            p_b <- 1
        }
        else if (l == 1) {
            p_b <- 0.5
        }
        else if (l < 1) {
            p_b <- 0
        }
        return(p_b)
    }

    if (action == 'q') {
        if (pi_per == 0) {
            emp_pi_t = 0
        } else {
            emp_pi_t <- (emp_at*pi_per) / (emp_at*pi_per + emp_bt*pi_per)
        }
        ## This is the updated belief that the state of the world is a,
        ## given that player1 chooses 'go'.
        l = 2*emp_pi_t
        r = (1-emp_pi_t)*A
        if ( l > r ) {
            q <- 1
        }
        else if ( l == r ) {
            q <- 0.5
        }
        else if ( l < r ) {
            q <- 0
        }
        return(q)
    }
}

model_4 <- function(pi_per, payoff, role, history=c("aGL","bGR"),
                    round, alpha, propA0=c(1,1), propB0=c(1,1),
                    propQ0=c(1,1), refPayoff=0) {
        
    ## INITIALIZE: equal probability of each action.

    s <- pi_per * 10

    ## Fix?
    if (s == 0) s = 0.0000000000000001
    
    propA <- propA0[1] * s
    propB <- propB0[1] * s
    propQ <- propQ0[1] * s
    propNotA <- propA0[2] * s
    propNotB <- propB0[2] * s
    propNotQ <- propQ0[2] * s

    pa <- propA / (propA + propNotA)
    pb <- propB / (propB + propNotB)
    pq <- propQ / (propQ + propNotQ)

    pNotA <- 1 - pa
    pNotB <- 1 - pb
    pNotQ <- 1 - pq

    if (role == "r") {
        finalAction <- substr(history[round], 1, 1)
    } else {
        finalAction <- "q"
    }
    
    ###########
    ## ROUND 1: We are done
    ###########
        
    if (round == 1) {
        if (finalAction == "a") {
			if(substr(history[round], 2, 2) == "G"){
	            return(pa)			
			} else{
				return(1-pa)
			}
        } else if (finalAction == "b") {
			if(substr(history[round], 2, 2) == "G"){
	            return(pb)			
			} else{
				return(1-pb)
			}
        } else {
			if(substr(history[round], 3, 3) == "L"){
	            return(pq)			
			} else{
				return(1-pq)
			}
        }
    }  
    
    ###########
    ## ROUND 2: 
    ###########

    if (round == 2 || role != "r") {
        action2 <- finalAction
    } else {
        action2 <- substr(history[2], 1, 1)
    }
    
    ## History of round 1.
    h1 <- history[1]
    ## Payoff of whatever was played.
    tmp <- payoff[[h1]]

    if (action2 == "a") {
        ## Player 1 payoff.
        Rx1 <- tmp[1] - refPayoff
        update <- Rx1 * alpha
        
        if (substr(h1, 2, 2) == "S") {
            propNotA <- propNotA + update
            if (propNotA <= 0) propNotA <- minPropensity
        } else {
            propA <- propA + update
            if (propA <= 0) propA <- minPropensity
        }
        pa <- propA / (propA + propNotA)
        currProb <- pa
        
    } else if (action2 == "b") {
        ## Player 1 payoff.
        Rx1 <- tmp[1] - refPayoff
        update <- Rx1 * alpha
        
        if (substr(h1, 2, 2) == "S") {
            propNotB <- propNotB + update
            if (propNotB <= 0) propNotB <- minPropensity
        } else {
            propB <- propB + update
            if (propB <= 0) propB <- minPropensity
        }
        pb <- propB / (propB + propNotB)
        currProb <- pb
        
    } else {
        ## Player 2 payoff.
        Rx2 <- tmp[2] - refPayoff
        update <- Rx2 * alpha
        
        if (substr(h1, 3, 3) == "R") {
            propNotQ <- propNotQ + update
            if (propNotQ <= 0) propNotQ <- minPropensity
        } else {
            propQ <- propQ + update
            if (propQ <= 0) propQ <- minPropensity
        }
        pq <- propQ / (propQ + propNotQ)
        currProb <- pq

    }

    ## Return results if round == 2.

    if (round == 2) {
        return(currProb)
    }

    ###########
    ## ROUND 3: 
    ###########

    ## History of round 2.
    h2 <- history[2]

    ## Payoff of whatever was played.
    tmp <- payoff[[h2]]

    if (finalAction == "a") {
        ## Player 1 payoff.
        Rx1 <- tmp[1] - refPayoff
        update <- Rx1 * alpha

        if (substr(h1, 2, 2) == "S") {
            propNotA <- propNotA + update
            if (propNotA <= 0) propNotA <- minPropensity
        } else {
            propA <- propA + update
            if (propA <= 0) propA <- minPropensity
        }
        pa <- propA / (propA + propNotA)
        return(pa)

    } else if (finalAction == "b") {
        ## Player 1 payoff.
        Rx1 <- tmp[1] - refPayoff
        update <- Rx1 * alpha

        if (substr(h1, 2, 2) == "S") {
            propNotB <- propNotB + update
            if (propNotB <= 0) propNotB <- minPropensity
        } else {
            propB <- propB + update
            if (propB <= 0) propB <- minPropensity
        }
        pb <- propB / (propB + propNotB)
        return(pb)

    } else {
        ## Player 2 payoff.
        Rx2 <- tmp[2] - refPayoff

        if (substr(h1, 3, 3) == "R") {
            propNotQ <- propNotQ + update
            if (propNotQ <= 0) propNotQ <- minPropensity
        } else {
            propQ <- propQ + update
            if (propQ <= 0) propQ <- minPropensity
        }   
        pq <- propQ / (propQ + propNotQ)
        return(pq)
    }
}
