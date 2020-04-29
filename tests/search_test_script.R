#### 2D TEST  ####
rm(list=ls())
library(fdaPDE)

GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic'

data(peak2Ddata)
mesh <- create.mesh.2D(nodes, order=2)


# plot(mesh)

# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh)

set.seed(5847947)

# Exact data - Nodes locations
fs.test <- function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = TRUE) # Ramsay test function
{

  q <- pi * r/2
  a <- d <- x * 0

  ind <- x >= 0 & y > 0
  a[ind] <- q + x[ind]
  d[ind] <- y[ind] - r

  ind <- x >= 0 & y <= 0
  a[ind] <- (-q - x[ind])
  d[ind] <- -r - y[ind]

  ind <- x < 0
  a[ind] <- -atan(y[ind]/x[ind]) * r
  d[ind] <-( sqrt(x[ind]^2 + y[ind]^2) - r )* (y[ind]/r0*(as.numeric(abs(y[ind])<=r0 & x[ind]>-0.5))+(as.numeric(abs(y[ind])>r0 || x[ind]<(-0.5))))

  ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - r0)^2)

  f <- a * b + d^2

  if (exclude)
    f[ind] <- NA

  attr(f, "exclude") <- ind
  f
}
# Exact data - Nodes locations
data_exact = fs.test(mesh$nodes[,1], mesh$nodes[,2], exclude = FALSE)

# Perturbed data - nodes locations
sd = 5*abs(max(data_exact) - min(data_exact))/100
data = data_exact + rnorm(length(data_exact), sd = sd)

# Generate locations different from nodes
loc1=cbind(runif(50,min=0.5,max=2.5),runif(50,min=0.1,max=0.9)) #braccio superiore
loc2=cbind(runif(50,min=0.5,max=2.5),runif(50,min=-0.9,max=-0.1)) # braccio inferiore
loc3=cbind(runif(50,min=-0.7,max=-0.1),runif(50,min=-0.5,max=0.5)) #circonferenza grande
noditest=mesh$nodes[1:50,] # alcune loc coincidenti con nodi
loc=rbind(loc1,loc2,loc3,noditest)

# Exact data - Locations different from nodes
data_exact2 = fs.test(loc[,1], loc[,2], exclude = FALSE)

# Perturbed data - Locations different from nodes
sd2 = 5*abs(max(data_exact2) - min(data_exact2))/100
data2 = data_exact2 + rnorm(length(data_exact2), sd = sd2)

# Set smoothing coefficient lambda
lambda = 10^-2

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 = smooth.FEM(observations = data,
                         FEMbasis = FEMbasis,
                         lambda = lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)
output_CPP1$fit.FEM$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points1)))
length(which(is.na(points2)))
diff
# above results should be all 0

### 2) Tree mesh: through smooth.FEM ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM$FEMbasis$mesh
output_CPP_tree<-smooth.FEM(observations=data,
                            FEMbasis=FEMbasis,
                            lambda=lambda,
                            GCV=GCVFLAG,
                            GCVmethod = GCVMETHODFLAG)

#reuse the tree information
points2=eval.FEM(output_CPP_tree$fit.FEM,
                 locations = mesh$nodes)

len = length(points1)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0

### 3) Tree mesh: through create.FEM.basis ####

#create tree mesh information in advance (default is FALSE):
#check tree mesh information under FEMbasis$mesh
#check tree information between "output_CPP_tree" and "FEMbasis2" are the same
FEMbasis2 = create.FEM.basis(mesh, saveTree=TRUE)
dim(FEMbasis2$mesh$node_box)
dim(output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box
compare2=FEMbasis2$mesh$node_box


diff = 0
for (i in 1:row) {
  for (j in 1:col) {
    diff = diff + abs(compare1[i,j]-compare2[i,j])
  }

}

diff
# above results should be all 0

#####4) Tree vs naive: locations diff from mesh nodes  ######
output_CPP1 = smooth.FEM(locations=loc,
                         observations = data2,
                         FEMbasis = FEMbasis,
                         lambda = lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0


### 5) Barycenter : through smooth.FEM ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM(locations=loc,
                         observations = data2,
                         FEMbasis=FEMbasis,
                         lambda = lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM(output_CPP_bary$fit.FEM,
                 bary.locations = output_CPP_bary$bary.locations)

#C++ requres to calcuate the barycenter
points3=eval.FEM(output_CPP_bary$fit.FEM,
                 locations = loc)

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points3[i]-points2[i])
}

length(which(is.na(points2)))
length(which(is.na(points3)))
diff
# above results should be all 0


### 6) Location comparison in locations vs bary.locations ####
# If locations in 'bary.locations' and 'locations' are the same, stops
loc2=loc
loc2[1,]=c(0,0)
output_CPP2_bary = smooth.FEM(locations=loc2,
                              bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM(output_CPP_bary$fit.FEM,
                locations=loc2,
                 bary.locations = output_CPP_bary$bary.locations)
# give error of 'Locations are not same as the one in barycenter information.'

### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM(bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)

#####.##########################################################
#### 2.5D TEST ####
library(fdaPDE)
rm(list=ls())

GCVFLAG=TRUE
GCVMETHODFLAG='Stochastic' # Stochastic GCV (default option)

data(hub25Ddata)
mesh <- create.mesh.2.5D(nodes = nodes,triangles = triangles)

# Creat FEMbasis object
FEMbasis <- create.FEM.basis(mesh)

# Set seed
set.seed(5847947)

# Exact test function - locations at nodes
nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
}
fem.test = FEM(func_evaluation, FEMbasis)


# Create observations for locations
x <- seq(-3,3,by=0.6)
y <- seq(-3,3,by=0.6)
z <- seq(-3,3,by=0.6)
grid = expand.grid(x=x, y=y, z=z)
loc = points.projection.2.5D(mesh, grid)

npoints = dim(loc)[1]
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
data2 = numeric(npoints)

for (i in 0:(npoints-1)){
  data2[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) +1
}


# Set smoothing parameter lambda
lambda=c(0.00375)

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 =smooth.FEM(observations = func_evaluation,
                        FEMbasis = FEMbasis,
                        lambda = lambda,
                        GCV = GCVFLAG,
                        GCVmethod = GCVMETHODFLAG)
output_CPP1$fit.FEM$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points1)))
length(which(is.na(points2)))
diff
# above results should be all 0

### 2) Tree mesh: through smooth.FEM ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM$FEMbasis$mesh
data = func_evaluation
output_CPP_tree<-smooth.FEM(observations=data,
                            FEMbasis=FEMbasis,
                            lambda=lambda,
                            GCV=GCVFLAG,
                            GCVmethod = GCVMETHODFLAG)

#reuse the tree information
points2=eval.FEM(output_CPP_tree$fit.FEM,
                 locations = mesh$nodes)

len = length(points1)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0

### 3) Tree mesh: through create.FEM.basis ####

#create tree mesh information in advance (default is FALSE):
#check tree mesh information under FEMbasis$mesh
#check tree information between "output_CPP_tree" and "FEMbasis2" are the same
FEMbasis2 = create.FEM.basis(mesh, saveTree=TRUE)
dim(output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box)
dim(FEMbasis2$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box
compare2=FEMbasis2$mesh$node_box


diff = 0
for (i in 1:row) {
  for (j in 1:col) {
    diff = diff + abs(compare1[i,j]-compare2[i,j])
  }

}

diff
# above results should be all 0


#####4) Tree vs naive: locations diff from mesh nodes   ######
output_CPP1 = smooth.FEM(locations=loc,
                         observations = data2,
                         FEMbasis = FEMbasis,
                         lambda = lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0

### 5) Barycenter : through smooth.FEM ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM(locations=loc,
                             observations = data2,
                             FEMbasis=FEMbasis,
                             lambda = lambda,
                             GCV=GCVFLAG,
                             GCVmethod = GCVMETHODFLAG)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM(output_CPP_bary$fit.FEM,
                 bary.locations = output_CPP_bary$bary.locations)

#C++ requres to calcuate the barycenter
points3=eval.FEM(output_CPP_bary$fit.FEM,
                 locations = loc)

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points3[i]-points2[i])
}

length(which(is.na(points2)))
length(which(is.na(points3)))
diff
# above results should be all 0


### 6) Location comparison in locations vs bary.locations ####
# If locations in 'bary.locations' and 'locations' are the same, stops
loc2=loc
loc2[1,]=c(0,0,0)
output_CPP2_bary = smooth.FEM(locations=loc2,
                              bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM(output_CPP_bary$fit.FEM,
                locations=loc2,
                bary.locations = output_CPP_bary$bary.locations)
# give error of 'Locations are not same as the one in barycenter information.'


### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM(bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)

#####.##########################################################
#### 3D TEST ####
library(fdaPDE)
rm(list=ls())

GCVFLAG=TRUE

GCVMETHODFLAG='Stochastic'

# Select mesh: sphere
selected_mesh="sphere"
data("sphere3Ddata")
mesh<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)

# Creat FEMbasis object
FEMbasis <- create.FEM.basis(mesh)

# Locations at nodes
nodesLocations=mesh$nodes

# Set seed
set.seed(5847947)

# Exact test function - locations at nodes
nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)

func_evaluation = numeric(nnodes)

for (i in 0:(nnodes-1)){
  func_evaluation[i+1] = a1* sin(2*pi*mesh$nodes[i+1,1]) +  a2* sin(2*pi*mesh$nodes[i+1,2]) +  a3*sin(2*pi*mesh$nodes[i+1,3]) +1
}
ran = range(func_evaluation)

# Perturbed data - Locations at nodes (sd is set according to range(func_evaluation))
data=func_evaluation+rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))

# Set smoothing parameter
lambda=c(10^-2)

set.seed(5847947) # necessario per la rimozione dei punti fuori dala mesh C

# Build locations different from nodes
if(selected_mesh=="sphere"){
  nloc = 1000
  loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T) # randomly generated points

  ind=NULL
  for(row in 1:nloc){
    normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
    if(normvec>0.975)   # remove points that fall outside the sphere
      ind = c(ind,row)
  }

  loc=loc[-ind,]
  nloc=dim(loc)[1]

}else if(selected_mesh=="C_3D"){
  loc1=cbind(runif(750,min=0.95,max=6), runif(750,min=-0.95,max=0.95),runif(750,min=0.5,max=2)) #braccio superiore
  loc2=cbind(runif(750,min=0.95,max=6), runif(750,min=-0.95,max=0.95),runif(750,min=-2,max=-0.5)) # braccio inferiore
  loc3=cbind(runif(750,min=-0.95,max=0.95), runif(750,min=-0.95,max=0.95),runif(750,min=-2,max=2)) #circonferenza
  loc=rbind(loc1,loc2,loc3)
  nloc=dim(loc)[1]

  #check if any point falls outside and remove it
  #checkout = smooth.FEM.basis(observations = func_evaluation2,locations=loc,
  #                            FEMbasis = FEMbasis, lambda = lambda,
  #                            CPP_CODE = TRUE, GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG)

  out=c(44,52,66,70,114,129,135,153,176,193,200,212,256,275,281,283,288,320,346,361,367,371,381,396,404,411,423,
        464,480,490,501,506,529,541,543,546,547,549,555,582,589,601,603,618,627,637,673,686,705,718,729,741,747,
        752,765,766,770,781,783,787,789,790,802,809,817,826,865,885,903,920,935,940,945,988,989,1001,1004,1026,1042,
        1054,1069,1076,1111,1118,1122,1132,1200,1209,1212,1213,1214,1223,1251,1273,1306,1311,1341,1358,1364,1367,
        1375,1387,1388,1396,1425,1445,1453,1457,1494,1503,1504,1508,1509,1511,1514,1517,1521,1522,1524,1525,1530,1531,
        1533,1539,1543,1548,1549,1551,1553,1556,1557,1558,1560,1565,1566,1571,1576,1578,1583,1584,1586,1587,1588,1589,
        1596,1611,1625,1631,1634,1638,1639,1642,1646,1650,1651,1657,1658,1662,1670,1672,1674,1676,1685,1695,1700,1703,
        1704,1709,1714,1715,1718,1722,1724,1726,1729,1730,1736,1737,1738,1739,1741,1745,1747,1748,1751,1762,1765,1768,
        1770,1771,1773,1774,1776,1779,1781,1782,1784,1787,1788,1790,1792,1794,1799,1800,1801,1812,1814,1818,1821,1822,
        1823,1824,1827,1829,1831,1836,1839,1844,1846,1849,1850,1853,1855,1858,1859,1860,1861,1864,1865,1866,1868,1871,
        1875,1877,1878,1889,1891,1894,1895,1901,1905,1907,1912,1914,1922,1925,1928,1930,1934,1936,1938,1939,1943,1945,
        1946,1950,1953,1955,1958,1962,1964,1968,1970,1973,1976,1977,1979,1980,1987,1989,1991,1992,1999,2000,2001,2002,
        2004,2005,2008,2009,2012,2014,2015,2021,2032,2042,2043,2044,2045,2047,2051,2052,2055,2057,2064,2070,2078,2080,
        2082,2088,2089,2093,2096,2100,2101,2103,2104,2108,2109,2113,2128,2130,2132,2133,2140,2143,2147,2149,2152,2154,
        2156,2157,2163,2164,2166,2174,2177,2179,2182,2184,2185,2192,2193,2199,2202,2203,2204,2211,2215,2218,2227,2232,
        2234,2237,2239,2240,2244,2246)

  #update locations and field
  loc=loc[-out,]
  nloc=dim(loc)[1]
}

# Exact test function - locations different from nodes
func_evaluation2=numeric(nloc)
for (i in 0:(nloc-1)){
  func_evaluation2[i+1] = a1* sin(2*pi*loc[i+1,1]) +  a2* sin(2*pi*loc[i+1,2]) +  a3*sin(2*pi*loc[i+1,3]) +1
}
ran2=range(func_evaluation2)

# check if all point shave been correctly removed
# checkin =smooth.FEM.basis(observations = func_evaluation2,locations=loc,
#                           FEMbasis = FEMbasis, lambda = lambda,
#                           GCV = GCVFLAG, GCVmethod = GCVMETHODFLAG)


# Perturbed data - locations different from nodes
data2=func_evaluation2+rnorm(nloc,mean=0,sd=0.05*(ran2[2]-ran2[1]))

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 =smooth.FEM(observations = func_evaluation,
                        FEMbasis = FEMbasis,
                        lambda = lambda,
                        GCV=GCVFLAG,
                        GCVmethod = GCVMETHODFLAG)
output_CPP1$fit.FEM$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points1)))
length(which(is.na(points2)))
diff

# above results should be all 0

### 2) Tree mesh: through smooth.FEM ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM$FEMbasis$mesh
data = func_evaluation
output_CPP_tree<-smooth.FEM(observations=data,
                            FEMbasis=FEMbasis,
                            lambda=lambda,
                            GCV=GCVFLAG,
                            GCVmethod = GCVMETHODFLAG)

#reuse the tree information
points2=eval.FEM(output_CPP_tree$fit.FEM,
                 locations = mesh$nodes)

len = length(points1)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0


### 3) Tree mesh: through create.FEM.basis ####

#create tree mesh information in advance (default is FALSE):
#check tree mesh information under FEMbasis$mesh
#check tree information between "output_CPP_tree" and "FEMbasis2" are the same
FEMbasis2 = create.FEM.basis(mesh, saveTree=TRUE)
dim(FEMbasis2$mesh$node_box)
dim(output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM$FEMbasis$mesh$node_box
compare2=FEMbasis2$mesh$node_box


diff = 0
for (i in 1:row) {
  for (j in 1:col) {
    diff = diff + abs(compare1[i,j]-compare2[i,j])
  }

}

diff
# above results should be all 0


#####4) Tree vs naive: locations diff from mesh nodes   ######
output_CPP1 = smooth.FEM(locations=loc,
                         observations = data2,
                         FEMbasis = FEMbasis,
                         lambda = lambda,
                         GCV=GCVFLAG,
                         GCVmethod = GCVMETHODFLAG)

#tree search (default)
points1=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes)

#naive search
points2=eval.FEM(output_CPP1$fit.FEM,
                 locations = mesh$nodes,
                 search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0


### 5) Barycenter : through smooth.FEM ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM(locations=loc,
                             observations = data2,
                             FEMbasis=FEMbasis,
                             lambda = lambda,
                             GCV=GCVFLAG,
                             GCVmethod = GCVMETHODFLAG)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM(output_CPP_bary$fit.FEM,
                 bary.locations = output_CPP_bary$bary.locations)

#C++ requres to calcuate the barycenter
points3=eval.FEM(output_CPP_bary$fit.FEM,
                 locations = loc)

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points3[i]-points2[i])
}

length(which(is.na(points2)))
length(which(is.na(points3)))
diff
# above results should be all 0

### 6) Location comparison in locations vs bary.locations ####
# If locations in 'bary.locations' and 'locations' are the same, stops
loc2=loc
loc2[1,]=c(0,0,0)
output_CPP2_bary = smooth.FEM(locations=loc2,
                              bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM(output_CPP_bary$fit.FEM,
                locations=loc2,
                bary.locations = output_CPP_bary$bary.locations)
# give error of 'Locations are not same as the one in barycenter information.'

### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM(bary.locations =output_CPP_bary$bary.locations,
                              observations = data2,
                              FEMbasis = FEMbasis,
                              lambda = lambda,
                              GCV=GCVFLAG,
                              GCVmethod = GCVMETHODFLAG)
