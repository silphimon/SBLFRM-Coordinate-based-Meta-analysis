############################################
##  Code to produce plots as in the paper ##
############################################
rm(list = ls())
# Set seed 
set.seed(50)
# Set working directory 
setwd("...")       
require(oro.nifti)
norm_vec <- function(x) sqrt(sum(x^2))  

#########################
##  Read in the data   ##
#########################
data <- read.table('foci.txt', header = F, sep = "\t")  
data <- data[, c(2:4, 1)]
names(data) <- c('X','Y', 'Z', 'study')          # Coordinates and study number
attach(data)
data[1:10, ]

N <- length(unique(data[, 4]))
ni <- tabulate(data$study)    		  # Number of foci per study
exclude <- which(ni == 0)           # Contrasts to exclude because no foci are left
ni <- ni[which((tabulate(data$study) != 0) == TRUE)]
ID <- rep(1:N, ni)
d <- dim(data)[2] - 1  				      # Number of dimensions
data <- cbind(data[, c(1:d)], ID)

##################
##  Create axes ##
##################
mask_full = readNIfTI('brainmask.nii') # Finer grid used for plotting: 2 x 2 x 2 mm brain mask (voxels 91 x 109 x 91)
# attributes(mask_full)
# mask_centre <- c(90, -126, -72)
Ix <- seq(1, 182, by = 2)        
Iy <- seq(1, 218, by = 2)      
Iz <- seq(1, 182, by = 2)       
A <- (Iy[2] - Iy[1])*(Ix[2] - Ix[1])*(Iz[2] - Iz[1])       # Volume of each voxel

# Axes to reproduce plots similar to those in the paper (make it a squared figure)
Ixplotting <- seq(-109, 107, by = 2)
Iyplotting <- seq(-127, 89, by = 2)
Izplotting <- seq(-91, 125, by = 2)

Ixresc <- seq(-91, 89, by = 2)
Iyresc <- seq(-127, 89, by = 2)
Izresc <- seq(-73, 107, by = 2)




########################################################
##  Define radial cubic B-splines or Gaussian kernels ##
########################################################
# ND: mask is 91 x 109 x 91, but mask is null above z = 82 
# Discard slices above z = 79 provided there are no foci falling outside the mask				
Grid <- c()  	
grid_slice <- as.matrix(expand.grid(Ix,Iy))
for (n in 17:48){	# 1:91
  map = mask_full[,,n]
  map = map > 0.5
  msk = as.vector(map)
  Grid <- rbind(Grid, cbind(grid_slice[msk, ], Iz[n]))
}
dim(Grid)
rm(grid_slice)

##################################################
##   Restrict computation to points in the mask ##
##    to avoid bias from off-mask estimation    ##
##################################################
xx <- seq(40, 145, length.out = 8)
yy <- seq(40, 180, length.out = 8)
zz <- seq(38, 90, length.out = 7)
kernels <- expand.grid(xx, yy, zz)
knots <- kernels

# At every odd numbered z slice, shift the x coordinate to re-create a chess-like kernel grid
for (n in 1:floor(length(zz)/2)){
  kernels[which(kernels[, 3] == zz[n*2]), 1] <- kernels[which(kernels[, 3] == zz[2*n]), 1] + 10
}

knots <- kernels[- which(kernels[, 1] > 145), ] 

nu <- 1/(2*256)
d <- 3
B.pred <- matrix(0, dim(Grid)[1], dim(knots)[1]) # Observed matrix of basis functions
for (i in 1:dim(B.pred)[1]){
  obs.knot <- matrix(as.numeric(Grid[i, 1:d]), dim(knots)[1], d, byrow = TRUE) - knots
  B.pred[i, ] <- as.matrix(exp(-nu * apply(obs.knot, 1, norm_vec)^2))
}
B.pred <- cbind(1, B.pred)    				 # Dense basis for integral evaluation

for (i in 1:dim(B.pred)[2]){
  ind = which(B.pred[, i] < 1e-35)
  if (length(ind) > 0) B.pred[ind, i] = 0
}
p = dim(B.pred)[2]

#######################################
##      Define global constants      ##
#######################################
nrun = 50000;
burn = 25000;
thin = 50;
every = 100;
start = 250;    						    # Starting iteration to update the stepsize
sp = (nrun - burn)/thin;				# Number of posterior samples
maxk = 50                       # Must correspond to the same value Matlab code

######################################
##      Reading covariates in       ##
######################################
Z = read.table('covariates.txt');
r = dim(Z)[2];          							    # Number of covariates per study

######################################
##      Reading study-type in       ##
######################################
outcome <- read.table('studytype.txt'); Y = outcome
anger <- which(outcome[, 1] == 1)
disgust <- which(outcome[, 1] == 2)
fear <- which(outcome[, 1] == 3)
happy <- which(outcome[, 1] == 4)
sad <- which(outcome[, 1] == 5)

# binary <- matrix(NA, length(ID), 1)
# for (i in 1:N) binary[which(ID == i),1] <- outcome[i, 1]
# anger_studies <- data[binary == 1, 1:d]
# disgust_studies <- data[binary == 2, 1:d] 
# fear_studies <- data[binary == 3, 1:d]
# happy_studies <- data[binary == 4, 1:d] 
# sad_studies <- data[binary == 5, 1:d] 

train <- as.matrix(as.numeric(read.table(file = "train.txt", header = FALSE, sep = ",")))
test <- setdiff(1:N, train)


##################################################
##    Code for 3D plotting of the intensities   ##
##  Need to extract the basis evaluated at the  ##
##               matching z-slice               ##
##################################################
library(fields)
theta <- as.matrix(read.table(file = "thetacbma_post.txt", header = FALSE, sep = " "))





##############################################
## Plotting posterior means log intensities ##
##############################################
type1 <- which(Y == 1); 
type1_cbma <- matrix(0, sp, p)
for (i in 1:length(type1)){
  j = type1[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type1_cbma <- type1_cbma + theta[, ind] 
}
type1_cbma <- type1_cbma/length(type1)
intensity1_cbma <- apply((B.pred %*% t(type1_cbma)), 1, mean)

type2 <- which(Y == 2); 
type2_cbma <- matrix(0, sp, p)
for (i in 1:length(type2)){
  j = type2[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type2_cbma <- type2_cbma + theta[, ind] 
}
type2_cbma <- type2_cbma/length(type2)

intensity2_cbma <- apply((B.pred %*% t(type2_cbma)), 1, mean)

type3 <- which(Y == 3); 
type3_cbma <- matrix(0, sp, p)
for (i in 1:length(type3)){
  j = type3[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type3_cbma <- type3_cbma + theta[, ind] 
}
type3_cbma <- type3_cbma/length(type3)

intensity3_cbma <- apply((B.pred %*% t(type3_cbma)), 1, mean)

type4 <- which(Y == 4); 
type4_cbma <- matrix(0, sp, p)
for (i in 1:length(type4)){
  j = type4[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type4_cbma <- type4_cbma + theta[, ind] 
}
type4_cbma <- type4_cbma/length(type4)

intensity4_cbma <- apply((B.pred %*% t(type4_cbma)), 1, mean)

type5 <- which(Y == 5); 
type5_cbma <- matrix(0, sp, p)
for (i in 1:length(type5)){
  j = type5[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type5_cbma <- type5_cbma + theta[, ind] 
}
type5_cbma <- type5_cbma/length(type5)
intensity5_cbma <- apply((B.pred %*% t(type5_cbma)), 1, mean)


## Choose some axial slices ##
set = c( 20, 24, 28, 32, 38, 45)

for (i in 1:length(set)){
  n = set[i]
  # This is the z coordinate in mms of what I'm plotting (the center of the voxel)
  map = mask_full[,,n]
  map = map > 0.5
  msk = as.vector(map)
  mask.pos <- which(map == TRUE, arr.ind = TRUE) 
  int1 <- intensity1_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int1[f]) 
  }
  image_matrix1 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale1 <- range((image_matrix1), na.rm = TRUE)
  
  int2 <- intensity2_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int2[f]) 
  }
  image_matrix2 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale2 <- range((image_matrix2), na.rm = TRUE)
  
  int3 <- intensity3_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int3[f]) 
  }
  image_matrix3 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale3 <- range((image_matrix3), na.rm = TRUE)
  
  int4 <- intensity4_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int4[f]) 
  }
  image_matrix4 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale4 <- range((image_matrix4), na.rm = TRUE)
  
  int5 <- intensity5_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int5[f]) 
  }
  image_matrix5 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale5 <- range((image_matrix5), na.rm = TRUE)
  
  out_cbma <- c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5)), 
                max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5)))
  
  pdf(file = paste("Z = ", Izresc[n], ".pdf"), width = 14, height = 4.5)
  m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.20, 0.90), c(0, 1, 0, 0.20))
  split.screen(m)
  screen(1)
  par(mar = c(0,1,1,1), oma = c(0, 0, 1, 0))
  title(main = paste("Post mean log intensities Z = ", Izresc[n]))
  
  split.screen(c(1,5), screen=2)-> ind2
  j = 1
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  # par(mfrow = c(2,3), pty = "s", mar = c(2,4,2,3))
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix1), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "Anger", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
  j = 2
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix2), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "Disgust", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
  
  j = 3
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix3), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "Fear", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
  
  j = 4
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix4), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "Happy", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
  
  j = 5
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix5), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "Sad", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
  
  screen(3)
  par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
  image.plot(zlim = out_cbma, horizontal = T, legend.only = TRUE, 
             legend.width = .6, col=gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2))
  close.screen( all=TRUE)
  
  dev.off()
}






# # # # # # # # # # # # # # # # # 
# # -- Examine learnt bases -- # 
# # # # # # # # # # # # # # # #
factor <- as.matrix(read.table(file = "Factor.txt", header = FALSE, sep = " ", skip = burn/every))
Lambdaout <- as.matrix(read.table(file = "Lambda.txt", header = FALSE, sep = " "))

j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out = out + B.pred[ind, ] %*% lambda_i
}
out <- out / sp

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out2 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out2 = out2 + B.pred[ind, ] %*% lambda_i
}
out2 <- out2 / sp

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out3 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out3 = out3 + B.pred[ind, ] %*% lambda_i
}
out3 <- out3 / sp

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out4 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out4 = out4 + B.pred[ind, ] %*% lambda_i
}
out4 <- out4 / sp

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out5 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out5 = out5 + B.pred[ind, ] %*% lambda_i
}
out5 <- out5 / sp

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out6 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out6 = out6 + B.pred[ind, ] %*% lambda_i
}
out6 <- out6 / sp

pdf(file = paste("learnt1.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[1])))
# close.screen(all = TRUE)

split.screen(c(1,6), screen = 2)-> ind2

l = 1
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out1range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

g = 1
screen(ind2[g])
j = 1
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 2
screen(ind2[g])
j = 2
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 3
screen(ind2[g])
j = 3
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 4
screen(ind2[g])
j = 4
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 5
screen(ind2[g])
j = 5
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 6
screen(ind2[g])
j = 6
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

screen(3) 
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out1range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))
close.screen(all = TRUE)
dev.off()


pdf(file = paste("learnt2.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[2])))
# close.screen(all = TRUE)

split.screen(c(1,6), screen = 2)-> ind2
l = 2
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3) 
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))
close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt3.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[3])))

split.screen(c(1,6), screen = 2)-> ind2
l = 3
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt4.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.80, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[4])))

split.screen(c(1,6), screen = 2)-> ind2
l = 4
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt5.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.80, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[5])))

split.screen(c(1,6), screen = 2)-> ind2
l = 5
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE)  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()









#######################################################################
##       Create axial-slice-specific matrices of covariates          ##
## of dimension (# of voxels at slice z) x (# of saved factors (50)) ##
#######################################################################
## Load covariates ##
iota <- read.table('iota_PBIN.txt')     # iota is size r x k
Lambdaout <- as.matrix(read.table(file = "Lambda.txt", header = FALSE, sep = " "))

set = c( 20, 24, 28, 32, 38, 45)  # Choose some axial slices 
counter = 1                       # Change counter to select a different covariate (emotion meta-analysis: counter = 1:5)
pdf(file = paste("fMRIpet.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = "Modality: fMRI (0) vs PET (1)")

split.screen(c(1,6), screen = 2)-> ind2

j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  # REMOVE INTERCEPT FROM PLOTS!!!
  out = out + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out <- out / sp

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out2 <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out2 = out2 + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out2 <- out2 / sp

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out3 <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out3 = out3 + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out3 <- out3 / sp

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out4 <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out4 = out4 + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out4 <- out4 / sp

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out5 <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out5 = out5 + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out5 <- out5 / sp

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out6 <- matrix(0, length(ind), 1)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out6 = out6 + (B.pred[ind, 2:p] %*% lambda_i[2:p, ]) %*% as.matrix(iota_i[counter, ])
}
out6 <- out6 / sp

image_matrix1 = list()
scale1 = matrix(NA, nrow = 6, ncol = 2)
usethis = list(out, out2, out3, out4, out5, out6)

for (l in 1:6){
  ind = which(Grid[, 3] == Iz[set[l]])
  map = mask_full[,,set[l]]
  map = map > 0.5
  msk = as.vector(map)
  mask.pos <- which(map == TRUE, arr.ind = TRUE) 
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (j in 1:dim(mask.pos)[1]){
    pos <- mask.pos[j, ]
    spec = usethis[[l]]
    image_matrix[pos[1], pos[2]] <- spec[j, ]
  }
  image_matrix1[[l]] <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
  scale1[l, ] = range(image_matrix1[[l]], na.rm = T)
}
out1range = c(min(scale1[, 1]), max(scale1[, 2]))

for (g in 1:6) {
  screen(ind2[g])
  par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
  image(c(Ixplotting), c(Iyplotting), image_matrix1[[g]], xlab = ' ', ylab = ' ', axes = F,
        col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[g]]))
}
screen(3) 
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out1range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))
close.screen(all = TRUE)
dev.off()






########################
## Boxplot covariates ##
########################
out <- matrix(0, sp, r)
for(i in 1:sp){
  iota_i = matrix(as.numeric(iota[i, ]), nrow = r, ncol = maxk, byrow = F)
  out[i, ] = apply(iota_i, 1, norm_vec)
}

pdf(file = "Covariates.pdf", width = 10, height = 10)
boxplot(out, outline = F, names = c("Modality", "# of Subjects", "Inference method", "Stimulus", "p-correction"), main = " ")
dev.off()






###################################
# # Examine traceplots of theta # #
###################################
pdf(file = "coeff_cbma.pdf", width = 12, height = 8)
par(mfrow = c(2,3), pty = "s")
for (l in 1:6){ 
  i <- sample(1:N, 1)
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + i
  }
  th <- theta[, ind]
  j <- sample(1:p, 1)
  plot(1:sp, th[, j], type = "l", xlab = "Posterior samples", ylab = expression(theta[j]))
  title(main = paste("Study ", i))
}
dev.off()

V = dim(B.pred)[1]
pdf(file = "int_paper_cbma.pdf", width = 12, height = 16)
par(mfrow = c(4,3), pty = "s")
for (l in 1:12){ 
  i <- sample(1:N, 1)
  ind <- c()
  for (g in 1:p){
    k = g - 1
    ind[g] <- k*N + i
  }
  th <- theta[, ind]
  # est_int_cbma <- exp(B.pred %*% t(th))
  est_int_cbma <- B.pred %*% t(th)
  j <- sample(1:V, 1)
  plot(1:sp, est_int_cbma[j, ], type = "l", xlab = "Posterior samples", ylab = " ")
  title(main = paste("Log int at voxel ", j, ", study ", i))
}
dev.off()








########################################
##    Check Gelman-Rubin statistics   ## 
########################################
library("coda")
# -- # Requirement: several runs initialised at over-dispersed starting values # -- #
factor <- as.matrix(read.table(file = "Factor.txt", header = FALSE, sep = " ", skip = burn/every))
factor2 <- as.matrix(read.table(file = "~/.../Factor.txt", header = FALSE, sep = " ", skip = burn/every))
factor3 <- as.matrix(read.table(file = "~/.../Factor.txt", header = FALSE, sep = " ", skip = burn/every))
factor4 <- as.matrix(read.table(file = "~/.../Factor.txt", header = FALSE, sep = " ", skip = burn/every))
out <- as.mcmc((factor))
out2 <- as.mcmc(factor2); out3 <- as.mcmc(factor3); out4 <- as.mcmc(factor4)
outlist <- mcmc.list(out, out2, out3, out4)
gelman.diag(outlist, confidence = 0.975)
gelman.plot(outlist)




#####################################
##    Classification performance   ## 
#####################################
pred_cat <- read.table("pred_cat.txt")    # Category predicted for 
pred_prob <- read.table("pred_prob.txt")

library(xtable)
ytest <- Y[test, ]
ind1 <- which(ytest == 1)                 # Change this from 1:5 to assess prediction of the types
out <- 0; out2 <- 0; out3 <- 0; out4 <- 0; out5 <- 0
for (i in 1:dim(pred_cat)[1]){
  out <- out + length(which(pred_cat[i, ind1] == 1)) / length(ind1)   # Prop of "1" classified as "1"
  out2 <- out2 + length(which(pred_cat[i, ind1] == 2)) /length(ind1)  # Prop of "1" classified as "2"
  out3 <- out3 + length(which(pred_cat[i, ind1] == 3)) /length(ind1)  # Prop of "1" classified as "3"
  out4 <- out4 + length(which(pred_cat[i, ind1] == 4)) /length(ind1)  # Prop of "1" classified as "4"
  out5 <- out5 + length(which(pred_cat[i, ind1] == 5)) /length(ind1)  # Prop of "1" classified as "5"
}
# Averages over iterations
out <- out/dim(pred_cat)[1]
out2 <- out2/dim(pred_cat)[1]
out3 <- out3/dim(pred_cat)[1]
out4 <- out4/dim(pred_cat)[1]
out5 <- out5/dim(pred_cat)[1]
xtable(data.frame(c(out, out2, out3, out4, out5)))

#            Correct Classifications Rates
# Truth     Anger   Disgust   Fear    Happy    Sad  
# Anger     .....    .....    ....    ....     ...
# Disgust   .....    .....    ....    ....     ...
# Fear      .....    .....    ....    ....     ... 
# Happy     .....    .....    ....    ....     ...  
# Sad       .....    .....    ....    ....     ... 




#####################################
##    Study-specific intensities   ## The following lines convert data to/from different spaces 
#####################################
# -- # Convert coordinates to voxel space # -- #
origin = c(90,-126,-72)
data$X = round((origin[1]-data$X)/2)
data$Y = round((data$Y-origin[2])/2)
data$Z = round((data$Z-origin[3])/2)

# foci that fall outside [1,91]x[1,109]x[1,91] are removed
keep = (data$X>=1)&(data$X<=91)&(data$Y>=1)&(data$Y<=109)&(data$Z>=1)&(data$Z<=91)
# Restrict data to the amygdalae
# keep = (data$X>=1)&(data$X<=91)&(data$Y>=1)&(data$Y<=109)&(data$Z>=20)&(data$Z<=45)
data = data[keep,]

# Eliminate foci that fall outside of the mask
mask = readNIfTI('brainmask.nii')
mask = (mask>0)*1
for (i in dim(data)[1]:1) {
  if (mask[ data[i,1] , data[i,2] , data[i,3] ]==0){
    data = data[-i,]
  }
}

# Convert data to mm space 
xaxis <- seq(1, 182, by = 2)
yaxis <- seq(1, 218, by = 2)
zaxis <- seq(1, 182, by = 2)
data$X = xaxis[data$X]
data$Y = yaxis[data$Y]
data$Z = zaxis[data$Z]

# Rescale data: NOTE: might have to change z-coordinate if plotting something other than axial slices
resdata <- matrix(NA, sum(ni), 3)
for (i in 1:sum(ni)){
  ind <- c(which(Ix == data[i, 1]), which(Iy == data[i, 2]), which(Iz == data[i, 3]))
  resdata[i, ] <- c(Ixplotting[ind[1]+9], Iyresc[ind[2]], Izresc[ind[3]])
}
data <- cbind(resdata, ID)
data <- as.data.frame(data)
names(data) <- c('x','y', 'z', 'study')          # Coordinates of the foci and study number
attach(data)

# -- # Choose a study # -- #
i <- sample(1:N, 1)
th <- matrix(0, sp, p)
ind <- c()
for (l in 1:p){
  k = l - 1
  ind[l] <- k*N + i
}
th <- theta[, ind] 

# -- # Choose a slice: this is the z coordinate in mms of what I'm plotting (the center of the voxel) # -- #
n <- 40
grid_slice <- as.matrix(expand.grid(Ix,Iy))
map = mask_full[,,n]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE)  
mask_plot <- grid_slice[msk, ]
coords <- cbind(mask_plot, Iz[n])

B.plot <- matrix(0, dim(mask_plot)[1], dim(knots)[1]) # Observed matrix of basis functions
for (j in 1:dim(B.plot)[1]){
  obs.knot <- matrix(as.numeric(coords[j, 1:d]), dim(knots)[1], d, byrow = TRUE) - knots
  B.plot[j, ] <- as.matrix(exp(-nu * apply(obs.knot, 1, norm_vec)^2)) 
}
B.plot <- cbind(1, B.plot)
intensity <- apply(exp(B.plot %*% t(th)), 1, mean)
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))

for (f in 1:dim(mask.pos)[1]){
  pos <- mask.pos[f, ]
  image_matrix[pos[1], pos[2]] <- intensity[f] 
}
image_matrix2 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
out <- range(image_matrix2, na.rm = TRUE)

pp <- data[ID == i, 1:d]; pp <- pp[order(pp$z), ]
image.plot(c(Ixplotting), c(Iyplotting), image_matrix2, nlevel = 64, xlab = 'x', ylab = 'y', zlim = out,
           legend.shrink = 1, legend.width = 1, col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)))
if (sum(which(pp[, 3] == Izresc[n])) > 0)  points(pp[which(pp[, 3] == Izresc[n]), 1:2], pch = 16)
title(main = paste("Est int study", i, ", slice z = ", Izresc[n]+1, "mm"))


