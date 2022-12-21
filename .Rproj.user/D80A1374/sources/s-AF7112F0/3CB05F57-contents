# devtools::document()
# in NAMESPACE:
# useDynLib(pepit)
# importFrom(Rcpp, evalCpp)
# in rstudio:
# clean and rebuilt

# 2pdzA:B
library(pepit)
library(bio3d)
pdb=read.pdb("2pdz")
pdb=trim.pdb(pdb, string="protein")
pdb=trim.pdb(pdb, string="noh")

pdb=type_atoms(pdb)

target=access_atoms(pdb, chains="A", add="calpha", surface=FALSE)
bs=get_binding_sites(target, target_chains="A", ligand_chains="B", add="calpha", add.ligand=TRUE, surface=10)
target=surface_atoms(target,minacc=10)


bs=bio3d::trim.pdb(bs,chain="A")

X=target$atom[,c("x","y","z")]
X=as.matrix(X)
Y=bs$atom[,c("x","y","z")]
Y=as.matrix(Y)
XProp=target$atom$elesy
YProp=bs$atom$elesy
result=cliques(X, XProp, Y, YProp, deltadist=1.0, mindist=0.0, maxdist=25.0, types=c("a"))#si maxdist=15, pas tout

mapping=get_mapping(result, nrow(X))
I=mapping[[1]]$I
J=mapping[[1]]$J
length(I)
XProp[I]
YProp[J]
bio3d::rmsd(as.vector(t(X[I,])),as.vector(t(Y[J,])),fit=TRUE)

result=extend_cliques(X, XProp, Y, YProp, result, deltadist=1.0)
result=remove_redundant_clusters (clusters)
mapping=get_mapping(result, nrow(X))
I=mapping$cliques[[1]]$I
J=mapping$cliques[[1]]$J
bio3d::rmsd(as.vector(t(X[I,])),as.vector(t(Y[J,])),fit=TRUE)

pdb2=read.pdb("1g3j")
pdb2=bio3d::trim.pdb(pdb2,chain="A")
pdb2=type_atoms(pdb2)
target2=access_atoms(pdb2, chains="A", add="calpha", surface=FALSE)
target2=surface_atoms(target2,minacc=10)
X=target2$atom[,c("x","y","z")]
X=as.matrix(X)
XProp=target2$atom$elesy
result=cliques(X, XProp, Y, YProp, deltadist=1.0, mindist=0.0, maxdist=15.0, types=c("a"))#si maxdist=15, pas tout
mapping=get_mapping(result, nrow(X))
I=mapping[[1]]$I
J=mapping[[1]]$J
length(I)
XProp[I]
YProp[J]
bio3d::rmsd(as.vector(t(X[I,])),as.vector(t(Y[J,])),fit=TRUE)

