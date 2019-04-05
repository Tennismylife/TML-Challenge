SimulateChallenge <- function(n, weights1, weights2, matches_per_round, prob_bet_success=0.5, prob_walkover=0.025)
{

## prepare the containers for the results
res1 <- res2 <- matrix(NA, ncol=n, nrow=length(matches_per_round))
colnames(res1) <- paste0("p", seq(1,n))
colnames(res2) <- paste0("p", seq(1,n))
## container for the walkovers
walkovers <- vector(mode="list", length=length(matches_per_round))

## start the simulation round by round
for (round in seq_along(matches_per_round)) {
if (round == 1) {
## the first run can't have WO
walkovers[[round]] <- rep(1, matches_per_round[round])
} else {
walkovers[[round]] <- 1 - rbinom(matches_per_round[round], 1, prob_walkover)
}
## we simulate the outcomes as a binomial distribution for the matches actually played, and
## multiply the result for the weight of each round
bets <- rbinom(n, sum(walkovers[[round]]), prob_bet_success)
res1[round,] <- bets*weights1[round]
res2[round,] <- bets*weights2[round]
}
points_1 <- colSums(res1)
points_2 <- colSums(res2)
count_1 <- diff(rev(sort(points_1, decreasing=TRUE)[1:2])) > weights1[length(weights1)]
count_2 <- diff(rev(sort(points_2, decreasing=TRUE)[1:2])) > weights2[length(weights2)]
different_lead <- names(count_1)!=names(count_2)
## return a list
ret <- list(points_1=points_1, points_2=points_2, wo=sapply(walkovers, function(x) sum(1-x)),
count_1=count_1, count_2=count_2, different_lead=different_lead)

return(ret)
}

## simulations for Miami and IW
weights_top <- c(1, 1, 2, 4, 8, 16, 32)
weights_nd <- c(1, 1, 2, 3, 6, 11, 20)
mpr_miami <- c(32, 32, 16, 8, 4, 2)

miami <- replicate(10000, SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_miami,
prob_bet_success=0.55, prob_walkover=0.025), simplify=FALSE)

colSums(t(sapply(miami, function(x) return(c(x$count_1, x$count_2)), USE.NAMES = FALSE)))
sum(sapply(miami, function(x) return(x$different_lead)))

### simulations for the slams

weights_top <- c(1, 2, 4, 8, 16, 32, 64)
weights_nd <- c(1, 2, 3, 6, 11, 20, 37)
mpr_slam <- c(64, 32, 16, 8, 4, 2)

slams <- replicate(10000, SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_slam,
prob_bet_success=0.55, prob_walkover=0.025), simplify=FALSE)

colSums(t(sapply(slams, function(x) return(c(x$count_1, x$count_2)), USE.NAMES = FALSE)))
sum(sapply(slams, function(x) return(x$different_lead)))

