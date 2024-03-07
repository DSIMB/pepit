library(bio3d)

# In the case of Gly, which obviously lacks a Catom,
# we construct a pseudo-Catom by rotating the N atom over - 120° along the Ca-C axis
# Rw = cos(theta)v + (1-cos(theta))u%*%t(u)v + sin(theta) u%vec%v
# u = axis;norm(u)=1
# pseudo Cb calculation: 
# center at origin
#        n_v = n_v - ca_v
#        c_v = c_v - ca_v
#        # rotation around c-ca over -120 deg
#        rot = rotaxis(-pi * 120.0 / 180.0, c_v)
#        cb_at_origin_v = n_v.left_multiply(rot)
#        # move back to ca position
#        cb_v = cb_at_origin_v + ca_v
rotm120 = function(v, u){
	theta = -2*pi/3
	#ctheta = cos(theta)
	#stheta = sin(theta)
	ctheta = -0.5
	stheta = -sqrt(3)/2
	vec = cbind(u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3], u[1]*v[2]-u[2]*v[1])
	ctheta*v + (1-ctheta)*sum(u*v)*u + stheta*vec
}


HSEA = function(pdb, cutoff=13) {
     XCA = as.matrix(pdb$atom[pdb$calpha, c("x","y","z")])
     nres = unique(pdb$atom$resno)
     n = nrow(XCA)
     hseau = c()
     hsead = c()
     for (i in 2:(n-1)) {
         cacb = XCA[i+1,] - 2*XCA[i,] + XCA[i-1,]
    	 sel = atom.select(pdb, resno=nres[i], elety="CA")
    	 selA = atom.select(pdb, string="calpha")
    	 bs = binding.site(a=pdb, a.inds=selA, b.inds=sel, byres=FALSE, cutoff=cutoff)
	 ind = bs$inds$atom
    	 #bs = binding.site(a=pdb,  b.inds=sel, cutoff=cutoff)
	 #ind = intersect(bs$inds$atom, selA$atom)
	 Xnei = as.matrix(pdb$atom[ind, c("x","y","z")])
    	 Xnei = -t(t(Xnei) - XCA[i,])
	 up = Xnei%*%cacb
    	 hseau = c(hseau, sum(up>0)) #HSEAU
	 hsead = c(hsead, sum(up<0)) #HSEAD
    }
    seq = pdb$atom$resid[pdb$calpha]
    names(hseau) = seq[2:(n-1)]
    names(hsead) = seq[2:(n-1)]
    list(hseau=hseau, hsead=hsead)
}


HSEB = function(pdb, cutoff=13) {
     XCA = as.matrix(pdb$atom[pdb$calpha, c("x","y","z")])
     nres = pdb$atom$resno[pdb$calpha]
     n = nrow(XCA)
     hsebu = c()
     hsebd = c()
     seq = pdb$atom$resid[pdb$calpha]
     for (i in 1:n) {
     	 res = trim(pdb, resno=nres[i])
    	 if (res$atom$resid[1]=="GLY") {#pseudo-CB
 	   xn = res$atom[atom.select(res, elety="N")$atom, c("x", "y", "z")]
	   xc = res$atom[atom.select(res, elety="C")$atom, c("x", "y", "z")]
	   can = xn - XCA[i,]
	   cac = xc - XCA[i,]
	   cac = cac/sqrt(sum(cac**2))
	   cacb = rotm120(can, cac)
    	 }
	 if (res$atom$resid[1]!="GLY") {
       	   xcb = res$atom[atom.select(res, elety="CB")$atom, c("x", "y", "z")]
       	   xcb = t(xcb)
    	   cacb = xcb - XCA[i,]
   	 }
    	 sel = atom.select(pdb, resno=nres[i], elety="CA")
    	 selA = atom.select(pdb, string="calpha")
    	 bs = binding.site(a=pdb, a.inds=selA, b.inds=sel, byres=FALSE, cutoff=cutoff)
    	 Xnei = as.matrix(pdb$atom[bs$inds$atom, c("x","y","z")])
	 up = down = 0
	 for (j in 1:nrow(Xnei)) {
	     ps = sum((Xnei[j,]-XCA[i,])*cacb)
	     up = up + ifelse(ps<0,1,0)
	     down = down + ifelse(ps>=0,1,0)
	 }
   	 hsebu = c(hsebu, up)
    	 hsebd = c(hsebd, down)
     }
     names(hsebu) = seq
     names(hsebd) = seq
     list(hsebu=hsebu, hsebd=hsebd)
}

CN = function(pdb, cutoff=13) {
   pdb = trim(pdb, string="calpha")
   n = nrow(pdb$atom)
   CN = c()
   for (i in 1:n) {
       sel = atom.select(pdb, eleno=pdb$atom$eleno[i])
       bs = binding.site(a=pdb, b.inds=sel, cutoff=cutoff, byres=FALSE)
       CN = c(CN, length(bs$inds$atom)-1)
   }
   CN
}

