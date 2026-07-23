# from Propedia list to Propedia table file
L = readLines("Propedia3567.txt")
L = unlist(strsplit(L,split=", "))
L = strsplit(L,split="_")
cat("target tchain pchain\n", file="Propedia3567.dat")
for (i in 1:length(L)) {
    line = L[[i]]
    cat(line[1], line[3], line[2], "\n", append=TRUE,  file="Propedia3567.dat")
}

# D = read.table("Propedia3567.dat", header = TRUE)
# head(D)