## History.R

# Set to TRUE to compile all program functions with cmpfun.
COMPILE <- TRUE
# Check that cmpfun exists.
if (!exists('cmpfun')) {
  COMPILE <- FALSE
}

## library('combinat')
library('gtools')

source('R_CODE/matches.R')
#source('matches.R')

# All combinations in one round.
combinations <- expand.grid(world=c("a","b"),
                            p1=c("G", "S"),
                            p2=c("L", "R"))
combinationsV <- apply(combinations,1,paste,collapse="")

## All possible outcomes for the first round.

## You need to adjust this value to make functions
## makeDatasets2Players vs makeDatasets work.

nPlayers <- 2
## How many matches?
if (nPlayers == 10) {
  playersPerRole <- 5
} else if (nPlayers == 2) {
  playersPerRole <- 1
} else {
  stop(paste('invalid nPlayers:', nPlayers))
}

outcomesR1 <- combinations(8, playersPerRole, combinationsV)
N_HISTORIES_ROUND_1 <- nrow(outcomesR1)
allRows <- 1:N_HISTORIES_ROUND_1

## Returns a matrix where each row contains a possible dataset.
##
## A dataset consists in 10 entries of the the type aSR-aSR-bGL,
## representing the history of a player. The first 5 entries
## are for the Red players, and second 5 for the Blue players.
## If collapse is TRUE the entries are collapsed in one string
## and are separated by @.
makeDatasets <- function(collapse=TRUE) {
    counter <- 1
    nCols <- ifelse(collapse, 1, playersPerRole*2)
    datasets <- matrix(nrow=N_HISTORIES_ROUND_1^3, ncol=nCols)
    # Round 1.
    for (row in allRows) {
        ## Round 1 Blue players and Red players have same history.
        redGame <- outcomesR1[row, ]
        blueGame <- redGame
        ## Make a vector of lenght 10, we will add to each element later.
        game1 <- c(redGame, blueGame)
        # Round 2.
        for (row2 in allRows) {
            ## Round 2 we need to check how Red and Blue players are matched.
            ## Red players have "default" history.
            redGame2 <- outcomesR1[row2, ]
            blueGame2 <- matrix(nrow=playersPerRole)
            for (bluePlayer in 1:playersPerRole) {
                ## Blue players have index 6-10.
                matchedRed <- getMatchFor((playersPerRole+bluePlayer), 2)
                redAction <- redGame2[matchedRed]
                blueGame2[bluePlayer] = redAction
            }
            game2 <- c(redGame2, blueGame2)
            dataset2 <- paste(game1, game2, sep="-")
            # Round 3.
            for (row3 in allRows) {
                redGame3 <- outcomesR1[row3, ]
                blueGame3 <- matrix(nrow=playersPerRole)
                for (bluePlayer in 1:playersPerRole) {
                    ## Blue players have index 6-10.
                    matchedRed <- getMatchFor((playersPerRole+bluePlayer), 3)
                    redAction <- redGame3[matchedRed]
                    blueGame3[bluePlayer] = redAction
                }
                game3 <- c(redGame3, blueGame3)
                dataset3 <- paste(dataset2, game3, sep="-")
                if (collapse) dataset3 <- paste(dataset3, collapse="@")
                datasets[counter,] <- dataset3
                counter <- counter + 1
            }
        }
    }
    return(datasets)
}

# Returns a list of matrices.
makeDatasets2Players <- function() {
    counter <- 1
    nCols <- playersPerRole*2
    datasets <- list()
    dataset <- matrix(nrow=N_HISTORIES_ROUND_1^3, ncol=nCols)
    # Round 1.
    for (row in allRows) {
        ## Round 1 Blue players and Red players have same history.
        redGame <- outcomesR1[row, ]
        blueGame <- redGame
        ## Make a vector of lenght 10, we will add each element later.
        game1 <- c(redGame, blueGame)
        # Round 2.
        for (row2 in allRows) {
            ## Round 2 we need to check how Red and Blue players are matched.
            ## Red players have "default" history.
            redGame2 <- outcomesR1[row2, ]
            blueGame2 <- matrix(nrow=playersPerRole)
            for (bluePlayer in 1:playersPerRole) {
                matchedRed <- 1
                redAction <- redGame2[matchedRed]
                blueGame2[bluePlayer] = redAction
            }
            game2 <- c(redGame2, blueGame2)
            dataset2 <- c(game1, game2)
            # Round 3.
            for (row3 in allRows) {
                redGame3 <- outcomesR1[row3, ]
                blueGame3 <- matrix(nrow=playersPerRole)
                for (bluePlayer in 1:playersPerRole) {
                    matchedRed <- 1
                    redAction <- redGame3[matchedRed]
                    blueGame3[bluePlayer] = redAction
                }
                game3 <- c(redGame3, blueGame3)
                dataset3 <- c(dataset2, game3)
                # if (collapse) dataset3 <- paste(dataset3, collapse="@")
                datasets[[counter]] <- matrix(ncol=3, nrow=2, dataset3,
                                              byrow=TRUE)
                counter <- counter + 1
            }
        }
    }
    return(datasets)
}


## Returns 8^3*8^2 different games of where 1 Red player meets 3 Blue
## players, and at round 3 the game for Blue and Red is the same.
makeDatasets1FullGame <- function(collapse=TRUE) {
  pp <- permutations(8, 3, repeats.allowed=TRUE,
                     v=c("aGL", "bGL", "aGR", "bGR",
                       "aSL", "bSL", "aSR", "bSR"))
  len <- dim(pp)[1]
  allTheseRows <- 1:len
  if (collapse) {
    datasets <- matrix(nrow=N_HISTORIES_ROUND_1^5, ncol=1)
  } else {
    datasets <- list()
  }
  counter <- 1
  for (row in allTheseRows) {
    ## print(paste0('row1 ', row))
    redGame <- pp[row,]
    for (row2 in allTheseRows) {
      ## print(paste0('row2 ', row2))
      blueGame <- pp[row2,]
      ## The third round must be in common.
      if (blueGame[3] == redGame[3]) {
        if (collapse) {
          datasets[counter] <- paste(paste(redGame , collapse='-'),
                                     paste(blueGame, collapse='-'),
                                     sep='@')
        } else {
          datasets[[counter]] <- matrix(ncol=3, nrow=2,
                                        c(redGame, blueGame),
                                        byrow=TRUE)
        }
        counter <- counter + 1        
        ## print(paste0('counter ', counter))
      }
    }
  }
  return(datasets)
}




# Test
# datasets <- makeDatasets()
# cc <- makeDatasets2Players()
# for (row in 1:length(cc)) {
#   m1 <- cc[[row]]
#   for (row2 in 1:length(cc)) {
#     if (row != row2) {
#       m2 <- cc[[row2]]
#       if (identical(m1, m2)) {
#         stop('AH')
#       }
#     }
#   }
# }

## Test.
# write.csv(datasets, file="/tmp/aa.csv")


# Compile, if requested.
if (exists('COMPILE') && COMPILE) {
  makeDatasets <- cmpfun(makeDatasets)
  makeDatasets2Players <- cmpfun(makeDatasets2Players)
  makeDatasets1FullGame <- cmpfun(makeDatasets1FullGame)
}
