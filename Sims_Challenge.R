## Code to emulate and extend Top's java simulation
 
SimulateChallenge <- function(n, weights1, weights2, matches_per_round, prob_bet_success=0.5, prob_walkover=0.025) 
{
    ## prepare the containers for the results
    res1 <- res2 <- matrix(NA, ncol=n, nrow=length(matches_per_round))
    colnames(res1)   <- paste0("p", seq(1,n))
    colnames(res2)   <- paste0("p", seq(1,n))
    ## container for the walkovers
    walkovers <- vector(mode="list", length=length(matches_per_round))

    ## start the simulation round by round
    for (round in seq_along(matches_per_round)) {
        if (round == 1) {
            ## the first run can't have WO
            walkovers[[round]] <- rep(1, matches_per_round[round])
        } else {
            walkovers[[round]] <- 1 -  rbinom(matches_per_round[round], 1,  prob_walkover)
        }
        ## we simulate the outcomes as a binomial distribution for the matches actually played 
        ## imposing a constant probability for all of them, 
        ## and multiply the result for the weight of each round
        bets <-  rbinom(n, sum(walkovers[[round]]), prob_bet_success)
        res1[round,] <- bets*weights1[round]
        res2[round,] <- bets*weights2[round] 
    }
    ## the row of the matrices are the total points per round of each player/column
    ## so the colSums are the total points after the last simulated run
    points_1 <- colSums(res1)
    points_2 <- colSums(res2)
    sp_1 <- sort(points_1, decreasing=TRUE)
    sp_2 <- sort(points_2, decreasing=TRUE)

    ## counting the bad cases
    bc_1 <- diff(rev(sp_1[1:2])) >  weights1[length(weights1)]
    bc_2 <- diff(rev(sp_2[1:2])) >  weights2[length(weights2)]
    ## counting the differences in lead with the two weight systems (as TRUE/FALSE)
    different_lead <- names(bc_1)!=names(bc_2)
    ## total walkovers in each round this simulation, as vector
    total_walkovers <- sapply(walkovers, function(x) sum(1-x))
    ## position of first player excluded from win
    pos_win_1 <-  min(which(abs(sp_1 -sp_1[1]) > weights1[length(weights1)], useNames=FALSE))
    pos_win_2 <-  min(which(abs(sp_2 -sp_2[1]) > weights2[length(weights2)], useNames=FALSE))
    ## position of first player excluded from top ten
    pos_ten_1 <-  min(which(abs(sp_1[10:length(sp_1)] -sp_1[10]) > weights1[length(weights1)], useNames=FALSE))+9
    pos_ten_2 <-  min(which(abs(sp_2[10:length(sp_2)] -sp_2[10]) > weights2[length(weights2)], useNames=FALSE))+9
    ## return a list 
    ret <- list(points_1=points_1, points_2=points_2, wo=total_walkovers,
                count_1=bc_1, count_2=bc_2, 
                different_lead=different_lead,
                pos_win_1=pos_win_1, pos_win_2=pos_win_2,
                pos_ten_1=pos_ten_1, pos_ten_2=pos_ten_2)

     return(ret)
}

## simulations for Miami and IW
weights_top <- c(1, 1, 2, 4, 8, 16, 32)
weights_nd <- c(1, 1, 2, 3, 6, 11,  20)
mpr_miami <- c(32, 32, 16, 8, 4, 2)


miami <- replicate(10000, SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_miami, 
                                            prob_bet_success=0.55, prob_walkover=0.025), simplify=FALSE)
 
colSums(t(sapply(miami, function(x) return(c(x$count_1, x$count_2)), USE.NAMES = FALSE))) 
sum(sapply(miami, function(x) return(x$different_lead)))
a <- t(sapply(miami, function(x) return(c(x$pos_win_1, x$pos_win_2, x$pos_ten_1, x$pos_ten_2 ))))


### simulations for the slams

weights_top <- c(1, 2, 4, 8, 16, 32, 64)
weights_nd <- c(1,  2, 3, 6, 11, 20, 37)
mpr_slam <- c(64, 32, 16, 8, 4, 2)

slams <- replicate(10000, SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_slam, 
                                                prob_bet_success=0.55, prob_walkover=0.025), simplify=FALSE)

## parallel version
ParSim <- function(x) SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_slam, prob_bet_success=0.55, prob_walkover=0.025)
library(parallel)
slams <- mclapply(as.list(1:10000), ParSim, mc.cores=7)

colSums(t(sapply(slams, function(x) return(c(x$count_1, x$count_2)), USE.NAMES = FALSE))) 
sum(sapply(slams, function(x) return(x$different_lead)))



### Grid of probabilities
probs_success <- seq(0.4,0.8, by=0.05)
probs_wo <- seq(0.005, 0.05, by=0.005)
exp <- expand.grid(pb=probs_success, pw=probs_wo)

library(foreach)
library(doMC)
## I have 8 cores, I'll split 2 "outer" processes each of which will use 4 cores
registerDoMC(2)

## wrap the function in a parallel-friendly fashion
ParSimMulti <- function(x, y) SimulateChallenge(30, weights1=weights_top, weights2=weights_nd, matches_per_round=mpr_slam, prob_bet_success=exp[y,1], prob_walkover=exp[y,2])

## simulations with 10k realizations for each probability combinations
prob_res <- foreach (i=seq(1, nrow(exp))) %dopar% {
       s <- mclapply(as.list(1:10000), ParSimMulti, y=i, mc.cores=4)
       cat(paste("Simulation",i,"done\n"))
       return(s)
}

ExtractRes <- function(ls) {
    a <- colSums(t(sapply(ls, function(x) return(c(x$count_1, x$count_2)), USE.NAMES = FALSE))) 
    b <- sum(sapply(ls, function(x) return(x$different_lead)))/length(ls)
    pm <- t(sapply(ls, function(x) return(c(x$pos_win_1, x$pos_win_2, x$pos_ten_1, x$pos_ten_2 ))))
    pm[pm>30] <- 30
    p <- colMeans(pm)
    ret <- c(a,b, p)
    names(ret) <- c("bad1", "bad2","dl","pw1", "pw2", "pt1","pt2")
    return(ret)
}
        
ss <- cbind(exp, t(sapply(prob_res, ExtractRes)))

summary(ss$bad1/10000*100)
summary(ss$bad2/10000*100)
summary(ss$dl)

ts <- ss
ts$bad1 <- ts$bad1/10000*100
ts$bad2 <- ts$bad2/10000*100

p1 <- contourplot(bad1 ~ pb * pw, data = ts,
            at=seq(0,0.16,by=0.02), 
            colorkey=list(at=seq(0,0.16,by=0.02)), 
            col.regions=terrain.colors(9),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Bad cases with Top's system, percent")
p2 <- contourplot(bad2 ~ pb * pw, data = ts,
            at=seq(0,0.7,by=0.1), 
            colorkey=list(at=seq(0,0.7,by=0.1)), 
            col.regions=heat.colors(8),outer=TRUE,
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Bad cases with ND's system, percent")

p3 <- contourplot(I(dl*100) ~ pb * pw, data = ts,
            at=seq(10,15,by=1), 
            colorkey=list(at=seq(10,15,by=1)), 
            col.regions=heat.colors(5),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Different lead between the two systems, percent")

p4 <- contourplot( pw1 ~ pb * pw, data = ts,
            at=seq(10,30,by=2), 
            colorkey=list(at=seq(10,30,by=2)), 
            col.regions=heat.colors(10),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Rank of the first player who can't win, Top")
p5 <- contourplot( pw2 ~ pb * pw, data = ts,
            at=seq(10, 30, by=2), 
            colorkey=list(at=seq(10,30,by=2)), 
            col.regions=heat.colors(10),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Rank of the first player who can't win, ND")

p6 <- contourplot( pt1 ~ pb * pw, data = ts,
            at=seq(24,30,by=0.5), 
            colorkey=list(at=seq(24,30,by=0.5)), 
            col.regions=heat.colors(14),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Rank of the first player who can't get top ten, Top")

p7 <- contourplot( pt2 ~ pb * pw, data = ts,
            at=seq(24,30,by=0.5), 
            colorkey=list(at=seq(24,30,by=0.5)), 
            col.regions=heat.colors(13),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Rank of the first player who can't get top ten, Top")
          
            
       
p8 <- contourplot( I(pw1-pw2) ~ pb * pw, data = ts,
            col.regions=terrain.colors(15),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Mean rank difference between excluded from victory, Top - ND")

p9 <- contourplot( I(pt1-pt2) ~ pb * pw, data = ts,
            col.regions=terrain.colors(15),
            panel = latticeExtra::panel.2dsmoother,
            region = TRUE, ##cuts = 9, 
            xlab = "P of betting success",
            ylab = "P of walkovers",
            main="Mean rank difference between excluded from top ten, Top - ND")
           
require(gridExtra)
cairo_pdf("SimulationResults.pdf", width=15, height=24, family="serif")
grid.arrange(p1, p2, p4, p5, p6, p7, p8, p9, ncol=2)
dev.off()      

cairo_pdf("SimulationResults_different_lead.pdf", width=7.5, height=6, family="serif")
p3
dev.off()      
         
