#matches <- read.csv('../DATA_INPUT/matches.csv', header=T, strip.white=TRUE)
matches <- read.csv('../DATA_INPUT/matches.csv', header=T, strip.white=TRUE)

N_MATCHES <- length(unique(matches$match))

matches.matrix <- as.matrix(matches)

# Precompute matches rounds.

# Empty matches at round 1.
r1Matches <- matches[matches$round == 1,]
r2Matches <- matches[matches$round == 2,]
r3Matches <- matches[matches$round == 3,]

# Make only a vector for each round.
r1MM <- as.matrix(r1Matches[order(r1Matches$p1),4])
r2MM <- as.matrix(r2Matches[order(r2Matches$p1),4])
r3MM <- as.matrix(r3Matches[order(r3Matches$p1),4])

## Returns the matched player for a given round
getMatchFor <- function(player, round) {
    if (round == 1) return(r1MM[player])
    if (round == 2) return(r2MM[player])
    if (round == 3) return(r3MM[player])
    stop('wtf')
}

# Compile, if requested.
if (exists('COMPILE') & COMPILE) {
    getMatchFor <- cmpfun(getMatchFor)
}
