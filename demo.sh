#
# Creation de la banque de BS PepPro2024-2
#
cd ~/3D/PepPro
mkdir  PepPro2024-2
Rscript ~/pepit/inst/BuiltBSBank.R PepPro_v2_holo.dat PepPro2024-2
#
#
#
cd ~/3D/PepPro
Rscript ~/pepit/inst/pepit.R ./PepPro2024-2/1dkd_bA\:P.pdb A PepPro2024-2/1dkd_bA\:P.dat test # target = pdb
#
#
#
cd 3D
Rscript ~/pepit/inst/getBSResi.R PepPro2024-2/1dkd_bA\:P.dat
# A229,A230,A231,A234,A237,A238,A241,A257,A260,A261,A262,A264,A265,A268,A270,A271

mkdir 1dkd
Rscript ~/pepit/inst/pepit.R 1dkd A ~/3D/Propedia/Propedia1305 1dkd1/1dkd "A229,A230,A231,A234,A237,A238,A241,A257,A260,A261,A262,A264,A265,A268,A270,A271"


#/home/alma/guyon/3D/Propedia/Propedia1305/1a2xA:B.dat ./1dkd.cif 1 174 15 0.607 0.01392758 0.2798144 162.0418 0.9404762 7 
cd 1dkd1
Rscript ~/pepit/inst/getSequence.R 1dkd-1.pdb
# EKRNRAITARR
Rscript ~/pepit/inst/getSequence.R ~/3D/Propedia/Propedia1305/1a2xA:B.pdb
# EEKRNRAITARRQHLKSVMLQIAATELEKEE
grep ^100 1dkd.score
# 100 /home/alma/guyon/3D/Propedia/Propedia1305/2l35A:B.dat ./1dkd.cif 1 193 17 1.745 0.01578459 0.4628265 155.2431 0.8759259 1 
Rscript ~/pepit/inst/getSequence.R 1dkd-100.pdb
# DLVLTVLI
Rscript ~/pepit/inst/getSequence.R ~/3D/Propedia/Propedia1305/2l35A:B.pdb
# CSTVSPGVLAGIVVGDLVLTVLIALAVYFLGR
grep "^2 " 1dkd.score
# 2 /home/alma/guyon/3D/Propedia/Propedia1305/1akhB:A.dat ./1dkd.cif 1 88 18 1.601 0.01671309 0.3832396 199.8304 0.8986486 6
Rscript ~/pepit/inst/getSequence.R 1dkd-2.pdb
# FLEEVFRRKQSLNSKEKEEVAKK
Rscript ~/pepit/inst/getSequence.R ~/3D/Propedia/Propedia1305/1akhB:A.pdb
# ISPQARAFLEEVFRRKQSLNSKEKEEVAKKCGITPLQVRVWFINKRMRS
