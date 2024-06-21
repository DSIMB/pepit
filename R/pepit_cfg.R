#
# user defined parameters
#
cat("pepit config\n")
RESIDUES = ""
CONTACT = 6.0
PRECISION = 1.0
ADD = "calpha"
ACC = 10
POSE = TRUE
NBCLIQUES = 5
NBHITS = 20
SCORE = mapping_dist_sum2
MAXDELTADIST = 25
MINDIST = 0
MINCLIQUE = 4
BCMIN = 0
INTERCLIQUE = 4
NBNEI = 4
MAXEDGES = 3000000
TYPES = c("a","b","c","o","n","A","C","O","N")
MINSCORE  =  10
MAXCLASHES  =  10
radius = 1.4
PVALUE = FALSE
AAGAP = 10
PCCODE = c(A = 1, C = 3, D = 5, E = 5, F = 1, G = 2, H = 6, I = 1, K = 6, L = 1, M = 1, N = 4, P = 1, Q = 4, R = 6, S = 4, T = 4, V = 1, W = 2, X = 0, Y = 4)
HSECUTOFF = 13
PROTEIN = TRUE
## utile ?
EMPTY_GRAPH = -3
NO_CLIQUE = -2
TOO_BIG = -1
