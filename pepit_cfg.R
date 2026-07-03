#
# user defined parameters
#
message("Warning : default pepit config")
# user defined
POSE = TRUE
PEPSEQ = TRUE
PEPTRIM = TRUE
MINPEPLEN = 4
MAXPEPLEN = 50
MAXCLASHES  =  100
NBCLIQUES = 5
MINCLIQUESIZE = 6
MINALEN = 6 
NBHITS = 20
ALIGNFILE = FALSE
#
# for BuildBSBank.R
#
ADD = "residue"
CONTACT = 6.0
HSECUTOFF = 13
#
# algorithm parameters
#
TYPES = c("a","b","c","o","n","A","C","O","N")
PRECISION = 1.0
MAXDELTADIST = 25
MINDIST = 0
BCMIN = 0
INTERCLIQUE = 4
NBNEI = 4
EXTEND = FALSE
MAXEDGES = 3000000
PVALUE = FALSE
MODE = 1 # vertex: all feature's comparison
RESIDUES = ""
PROTEIN = TRUE
HSERADIUS = 13
SCORE = mapping_dist_sum2
radius = 1.4
AAGAP = 10
PCCODE = c(A = 1, C = 3, D = 5, E = 5, F = 1, G = 2, H = 6, I = 1, K = 6, L = 1, M = 1, N = 4, P = 1, Q = 4, R = 6, S = 4, T = 4, V = 1, W = 2, X = 0, Y = 4)
MINSCORE  =  0 # alen>=10
EXTEND = FALSE
SCORE = mapping_dist_sum2
