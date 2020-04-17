#### 2D TEST  ####
rm(list=ls())
library(fdaPDE)

data(simpleMesh2Ddata) # Load the mesh
mesh=create.mesh.2D(nodes=nodes, triangles = triangles, order=2)

FEMbasis = create.FEM.basis(mesh) # Create the FEM basis object
time_mesh=seq(0,2,length.out = 5)
set.seed(5847947)

f<-function(loc) {pi*sin(loc[,1])*loc[,3]} # Exact solution

# Exact data - Nodes locations
SpaceTimePoints=cbind(rep(mesh$nodes[,1],length(time_mesh)),
                      rep(mesh$nodes[,2],length(time_mesh)),
                      rep(time_mesh,each=nrow(mesh$nodes)))
data_exact = f(SpaceTimePoints)

# Perturbed data - Nodes locations
data = data_exact + rnorm(n = length(data_exact), sd = 0.01)
tdata=matrix(data,nrow(mesh$nodes),length(time_mesh))


# Locations different from nodes
xobs=runif(min=-0.5,max=0.5,n=80)
yobs=runif(min=-0.5,max=0.5,n=80)

# Exact data - locations different from nodes
SpaceTimePoints=cbind(rep(xobs,5),
                      rep(yobs,5),
                      rep(time_mesh,each=80))
data_exact = f(SpaceTimePoints)

# Perturbed data - locations different from nodes
data = data_exact + rnorm(n = length(data_exact), sd = 0.01)
tdata2=matrix(data,length(xobs),length(time_mesh))
loc = cbind(xobs,yobs)

lambdaS = 1e-3
lambdaT = 1e-3

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 = smooth.FEM.time(observations = tdata,
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis,
                              lambdaS = lambdaS,
                              lambdaT = lambdaT,
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
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

### 2) Tree mesh: through smooth.FEM.time ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM.time$FEMbasis$mesh
output_CPP_tree<-smooth.FEM.time(observations = tdata,
                                 time_mesh = time_mesh,
                                 FEMbasis = FEMbasis,
                                 lambdaS = lambdaS,
                                 lambdaT = lambdaT,
                                 GCV=FALSE)

#reuse the tree information
points2=eval.FEM.time(output_CPP_tree$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

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
dim(output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box
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
output_CPP1 = smooth.FEM.time(locations = loc,
                              observations = tdata2, 
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis, 
                              lambdaS = lambdaS,
                              lambdaT = lambdaT, 
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
                      search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0


### 5) Barycenter : through smooth.FEM.time ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM.time(locations = loc,
                                  observations = tdata2, 
                                  time_mesh = time_mesh,
                                  FEMbasis = FEMbasis, 
                                  lambdaS = lambdaS,
                                  lambdaT = lambdaT, 
                                  GCV=FALSE)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      bary.locations = output_CPP_bary$bary.locations,
                      locations = loc,
                      time.instants=time_mesh)

#C++ requres to calcuate the barycenter
points3=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      locations = loc,
                      time.instants=time_mesh)


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
output_CPP2_bary = smooth.FEM.time(locations = loc2,
                                   bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                     locations=loc2,
                     bary.locations = output_CPP_bary$bary.locations,
                     time.instants=time_mesh)
# give error of 'Locations are not same as the one in barycenter information.'

### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM.time(bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)
#####.##########################################################
#### 2.5D TEST ####
library(fdaPDE)
rm(list=ls())

data(hub25Ddata)
mesh <- create.mesh.2.5D(nodes = nodes,triangles = triangles)

FEMbasis <- create.FEM.basis(mesh)


lambdaS=1
lambdaT=1

# Locations at nodes
nodesLocations=mesh$nodes

# Exact data - Locations at nodes
nnodes = mesh$nnodes
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
time_mesh = 0:4

locations = cbind(rep(time_mesh,each=nnodes),
                  rep(nodesLocations[,1],length(time_mesh)),
                  rep(nodesLocations[,2],length(time_mesh)),
                  rep(nodesLocations[,3],length(time_mesh)))

func = function(x)
{
  (a1*sin(2*pi*x[,2])+a2*sin(2*pi*x[,3])+a3*sin(2*pi*x[,4])+1)*cos(x[,1])
}

func_evaluation = func(locations)
ran = range(func_evaluation)
data = func_evaluation +rnorm(nnodes,mean=0,sd=0.05*(ran[2]-ran[1]))
tdata = matrix(data,mesh$nnodes,length(time_mesh))

# Locations different from nodes
x <- seq(-3,3,by=0.8)
y <- seq(-3,3,by=0.8)
z <- seq(-3,3,by=0.8)
grid = expand.grid(x=x, y=y, z=z)
loc = points.projection.2.5D(mesh, grid)


# Create observations for locations
locations = cbind(rep(time_mesh,each=dim(loc)[1]),
                  rep(loc[,1],length(time_mesh)),
                  rep(loc[,2],length(time_mesh)),
                  rep(loc[,3],length(time_mesh)))


func_evaluation = func(locations)
ran = range(func_evaluation)
data = func_evaluation +rnorm(dim(loc)[1],mean=0,sd=0.05*(ran[2]-ran[1]))
tdata2 = matrix(data,dim(loc)[1],length(time_mesh))

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 = smooth.FEM.time(observations = tdata,
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis,
                              lambdaS = lambdaS,
                              lambdaT = lambdaT,
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
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

### 2) Tree mesh: through smooth.FEM.time ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM.time$FEMbasis$mesh
output_CPP_tree<-smooth.FEM.time(observations = tdata,
                                 time_mesh = time_mesh,
                                 FEMbasis = FEMbasis,
                                 lambdaS = lambdaS,
                                 lambdaT = lambdaT,
                                 GCV=FALSE)

#reuse the tree information
points2=eval.FEM.time(output_CPP_tree$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

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
dim(output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box)
dim(FEMbasis2$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box
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
output_CPP1 = smooth.FEM.time(locations = loc,
                              observations = tdata2, 
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis, 
                              lambdaS = lambdaS,
                              lambdaT = lambdaT, 
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
                      search="naive")


len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0

### 5) Barycenter : through smooth.FEM.time ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM.time(locations = loc,
                                  observations = tdata2, 
                                  time_mesh = time_mesh,
                                  FEMbasis = FEMbasis, 
                                  lambdaS = lambdaS,
                                  lambdaT = lambdaT, 
                                  GCV=FALSE)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      bary.locations = output_CPP_bary$bary.locations,
                      locations = loc,
                      time.instants=time_mesh)

#C++ requres to calcuate the barycenter
points3=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      locations = loc,
                      time.instants=time_mesh)

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
output_CPP2_bary = smooth.FEM.time(locations = loc2,
                                   bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                     locations=loc2,
                     bary.locations = output_CPP_bary$bary.locations,
                     time.instants=time_mesh)
# give error of 'Locations are not same as the one in barycenter information.'

### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM.time(bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)

#####.##########################################################
#### 3D TEST ####
library(fdaPDE)
rm(list=ls())

# Build mesh: Sphere
data(sphere3Ddata)
mesh<-create.mesh.3D(sphere3Ddata$nodes, sphere3Ddata$tetrahedrons)

FEMbasis <- create.FEM.basis(mesh)
nodesLocations=mesh$nodes
nnodes = mesh$nnodes
time_mesh = seq(0,1,length.out = 5)

Locations = cbind(rep(time_mesh,each=nnodes),
                  rep(nodesLocations[,1],length(time_mesh)),
                  rep(nodesLocations[,2],length(time_mesh)),
                  rep(nodesLocations[,3],length(time_mesh)))

# Exact test function
a1 = rnorm(1,mean = 1, sd = 1)
a2 = rnorm(1,mean = 1, sd = 1)
a3 = rnorm(1,mean = 1, sd = 1)
a4 = rnorm(1,mean = 1, sd = 1)

func = function(x)
{
  a1*sin(2*pi*(x[,1]*x[,2]))+a2*cos(2*pi*x[,2])+a3*cos(2*pi*x[,3])+a4*sin(2*pi*x[,4])
}

func_evaluation = func(Locations)
ran=range(func_evaluation)

# Generate locations
nloc = 1000
loc=matrix(data=runif(3*nloc, min=-1,max=1),nrow=nloc,ncol=3,byrow=T)

ind=NULL
for(row in 1:nloc){
  normvec = (loc[row,1]^2+loc[row,2]^2+loc[row,3]^2)
  if(normvec>0.975)   # check points outside the sphere and remove them
    ind = c(ind,row)
}

loc=loc[-ind,]
nloc=dim(loc)[1]
loc = cbind(rep(time_mesh,each=nloc),
            rep(loc[,1],length(time_mesh)),
            rep(loc[,2],length(time_mesh)),
            rep(loc[,3],length(time_mesh)))


# Exact test function - locations different from nodes
func_evaluation2=func(loc)

lambdaS=1
lambdaT=1
# GCVFLAG=T
# GCVMETHODFLAG='Stochastic'

# beta_exact= c(0.7,2.0)
ran = range(func_evaluation)
data = func_evaluation +rnorm(nrow(Locations),mean=0,sd=0.05*(ran[2]-ran[1]))

ran = range(func_evaluation2)
data_noloc = func_evaluation2 +rnorm(nrow(loc),mean=0,sd=0.05*(ran[2]-ran[1]))

tdata = matrix(data,nnodes,length(time_mesh))
tdata2 = matrix(data_noloc,nloc,length(time_mesh))
loc = loc[1:nloc,2:4]

##### 1) Tree vs naive search comparison : mesh nodes location  #####
output_CPP1 = smooth.FEM.time(observations = tdata,
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis,
                              lambdaS = lambdaS,
                              lambdaT = lambdaT,
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
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

### 2) Tree mesh: through smooth.FEM.time ####

#create tree mesh information (default):
#check tree mesh information under output_CPP_tree$fit.FEM.time$FEMbasis$mesh
output_CPP_tree<-smooth.FEM.time(observations = tdata,
                                 time_mesh = time_mesh,
                                 FEMbasis = FEMbasis,
                                 lambdaS = lambdaS,
                                 lambdaT = lambdaT,
                                 GCV=FALSE)

#reuse the tree information
points2=eval.FEM.time(output_CPP_tree$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

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
dim(output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box)

row =dim(FEMbasis2$mesh$node_box)[1]
col = dim(FEMbasis2$mesh$node_box)[2]
compare1=output_CPP_tree$fit.FEM.time$FEMbasis$mesh$node_box
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
output_CPP1 = smooth.FEM.time(locations = loc,
                              observations = tdata2, 
                              time_mesh = time_mesh,
                              FEMbasis = FEMbasis, 
                              lambdaS = lambdaS,
                              lambdaT = lambdaT, 
                              GCV=FALSE)
output_CPP1$fit.FEM.time$FEMbasis = FEMbasis #for now, don't use tree info

#tree search (default)
points1=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh)

#naive search
points2=eval.FEM.time(output_CPP1$fit.FEM.time,
                      locations = mesh$nodes,
                      time.instants=time_mesh,
                      search="naive")

len = length(points2)

diff = 0
for (i in 1:len) {
  diff = diff + abs(points1[i]-points2[i])
}

length(which(is.na(points2)))
diff
# above results should be all 0


### 5) Barycenter : through smooth.FEM.time ####
#create bary.locations information (default):
#bary.locations is the barycenter of the location inside a certain element (given "locations" option)
#check bary.locations information under output$bary.locations
output_CPP_bary = smooth.FEM.time(locations = loc,
                                  observations = tdata2, 
                                  time_mesh = time_mesh,
                                  FEMbasis = FEMbasis, 
                                  lambdaS = lambdaS,
                                  lambdaT = lambdaT, 
                                  GCV=FALSE)

#use bary.locations information directly (C++ doesn't calculate the barycenter)
points2=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      bary.locations = output_CPP_bary$bary.locations,
                      locations = loc,
                      time.instants=time_mesh)

#C++ requres to calcuate the barycenter
points3=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                      locations = loc,
                      time.instants=time_mesh)

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
loc2=loc
loc2[1,]=c(0,0,0)
output_CPP2_bary = smooth.FEM.time(locations = loc2,
                                   bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)
# give error of 'Locations are not same as the one in barycenter information.'

points=eval.FEM.time(output_CPP_bary$fit.FEM.time,
                     locations=loc2,
                     bary.locations = output_CPP_bary$bary.locations,
                     time.instants=time_mesh)
# give error of 'Locations are not same as the one in barycenter information.'

### 7) If locations is null but bary.locations is not null, use the locations in bary.locations  ####
#check that there is bary.locations information under output_CPP3_bary$bary.locations
output_CPP3_bary = smooth.FEM.time(bary.locations =output_CPP_bary$bary.locations,
                                   observations = tdata2, 
                                   time_mesh = time_mesh,
                                   FEMbasis = FEMbasis, 
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT, 
                                   GCV=FALSE)
