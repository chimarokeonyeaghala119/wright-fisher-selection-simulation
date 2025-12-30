install.packages("R6")
library(R6)


Population <- R6Class("Population",
                      public = list(
                        nInd = NULL,
                        mut.rate = NULL,
                        initial_population = NULL,
                        # selection.coeff = NULL,
                        history = NULL,
                        traj = NULL,
                        alleles =NULL,
                        selection.coeff.for.new.mutations= NULL,
                        
                        
                        
                        initialize = function(nInd, mut.rate, selection.coeff) {
                          self$nInd <- nInd
                          self$mut.rate <- mut.rate
                          # self$selection.coeff = selection.coeff
                          
                          # Create initial population (0s)
                          ancestor <- 0.5
                          self$alleles<- c(1:10000)
                          self$selection.coeff.for.new.mutations<- selection.coeff
                          pop <- rep(0, nInd)
                          self$initial_population <- pop
                          
                          # Store history as a list
                          self$history <- list(self$initial_population)
                        },
                        
                        generation = function() {
                          current_gen <- self$history[[length(self$history)]]
                          fitness <- rep(1, self$nInd)
                          for (i in 1:length(fitness)) {
                            if (current_gen[i] == 0) {
                              fitness[i] <- 1
                            } else {
                              fitness[i] <- 1 + self$selection.coeff.for.new.mutations[current_gen[i]]
                            }
                          }
                          
                          # ancestor<- 0.5
                          # fitness<- fitness- ancestor
                          # fitness <- rep(1, self$nInd)
                          # fitness[current_gen == 1] <- 1 + self$selection.coeff[1]
                          # fitness[current_gen == 2] <- 1 + self$selection.coeff[2]
                          # fitness[current_gen == 3] <- 1 + self$selection.coeff[3]
                          # fitness[current_gen == 4] <- 1 + self$selection.coeff[4]
                          # 
                          
                          probs <- fitness / sum(fitness)
                          nextgen <- sample(current_gen, size = self$nInd, replace = TRUE, prob = probs)
                          for (i in 1: length(self$initial_population)) {
                            x <-runif(1, min = 0, max = 1)
                            if(x < mut.rate){
                              nextgen[i] <- sample(self$alleles, size = 1)
                            }
                          }
                          
                          return(nextgen)
                        },
                        
                        evolve = function(nGen) {
                          for (i in 1:nGen) {
                            new_gen <- self$generation()
                            self$history[[length(self$history) + 1]] <- new_gen
                          }
                          self$getTraj()
                         },
                        
                        getTraj = function() {
                          hist_matrix <- do.call(rbind, self$history)
                          alleles <- sort(unique(as.vector(hist_matrix)))
                          freq_matrix <- matrix(0, nrow = nrow(hist_matrix), ncol = length(alleles))
                          colnames(freq_matrix) <- alleles
                          for (gen in 1:nrow(hist_matrix)) {
                            counts <- table(factor(hist_matrix[gen, ], levels = alleles))
                            freqs <- counts / self$nInd
                            freq_matrix[gen, ] <- freqs
                          }
                          
                          self$traj <- freq_matrix
                          return(freq_matrix)
                        },
                        
                        plotTraj = function() {
                          if (is.null(self$traj)) {
                            self$getTraj()
                          }
                          
                          matplot(self$traj, type = "l", lty = 1, lwd = 2,
                                  col = rainbow(ncol(self$traj)),
                                  xlab = "Generation", ylab = "Allele Frequency",
                                  main = paste("Allele Freq.Traj for", nInd , "individuals", " mut.rate of", mut.rate, "and", nGen, 'gens'))
                          
                          legend("topright", legend = colnames(self$traj),
                                 col = rainbow(ncol(self$traj)), lty = 1, lwd = 2,  cex = 0.6)
                        }
                      )                        
)
###runing the sims
nInd = 1000
nGen = 100
nRuns = 100
mut.rate = 0.1
ancestor = 0.5
set.seed(123)
selection_coefficients<- rgamma(10000, 7, 7.5) - ancestor

pops <- lapply(1:nRuns, function(i) {
  Population$new(nInd = nInd,  mut.rate = mut.rate, selection.coeff = selection_coefficients)
})
for(pop in pops){
  pop$evolve(nGen)
}
traj<-pop$getTraj()
pop$plotTraj()


top.alleles.every.gen<- data.frame(generation = integer(),
                                   top.allele = numeric(),
                                   frequency = numeric())

gen.check <- seq(1, nrow(traj))
for(g in gen.check){
  top.alleles <- colnames(traj)[which.max(traj[g,])]
  freq<- max(traj[g,])
  top.alleles.every.gen<- rbind(top.alleles.every.gen, 
                       data.frame(generation = g,
                                  top.allele = top.alleles,
                                  frequency = freq))
  
}
print(top.alleles.every.gen)
