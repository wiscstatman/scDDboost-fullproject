
# Jan 4, 2019
# returning to Madison from 12 days in Big Sky

# Xiuyu and I having some good ideas on how to make a version of EBarrays or EBseq
# for a large number of groups, where, ordinarily, the exponential number of partitions
# would stymie our efforts

# idea

# given some distance matrix between groups (e.g. dist on average profiles)

# use Xiuyu's distance randomization method to get randomized distances

# From each randomized distance matrix, do a hierarchical clustering, and
# record the hierarchical set of partitions from each one.   
# Next, build up a probability distribution over partitions from the marginal
# record of all these hierarchical sets.   We'll always have the null partition
# (one block) and always have the complete partition (every profile its own block),
# and we'll have some distribution of relevant partitions from the full set
# of possible partitions.   For the EB* computations, the idea is to use the top
# say 1-epsilon fraction of this distribution, ranked by frequency, as input into
# the set of plausible patterns


# To test the idea, I first want to see if it's easy to extract partitions from
# a hierarchical set.

set.seed(75751)

x <- rnorm(20)

dd <- dist(x)

h <- hclust(dd, method="average" )

parts <- vector("list",length=(length(x)))

for( k in 2:(length(x)-1) )
 {
  r <- rect.hclust( h, k=k )
  parts[[k]] <- r 
 }

## it is...so let's try the randomization


set.seed(75751)

x <- rnorm(20)

dd <- as.matrix( dist(x) )

B <- 1000
parts <- vector("list",length=(B*(length(x)-2)) )

counter <- 1

for( b in 1:B )
 {
   e <- rgamma( length(x), shape=1/2, scale=1 )
   w <- 1/( e %o% e )

   ddstar <-  as.dist(dd*w)

   hstar <- hclust(ddstar, method="average" )

   for( k in 2:(length(x)-1) )
    {
     r <- rect.hclust( hstar, k=k )
     parts[[counter]] <- r
     counter <- counter + 1
    }

  }

## slow owing to plotting of `rect.hclust`...surely a faster way!

