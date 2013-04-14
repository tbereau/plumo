# Zun-Jing Wang 
# 2011 August    

# WALP23 peptide
# GLY TRP TRP LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU TRP TRP ALA
# NO ionizable residue.
# The atom types in this peptide are respectively 8:N, 9:CA, 10: Pro-N (no H-bond), 11:C, 12:O, 15:term-N, 16:term-C
# 20-39: CB (Side chain beads)
#
############
# WARNING: ionizable residues are parametrized in their charged form.
#   *** END CAPS ADDED ***
###########
#
#
# CAP _40_
# GLY 20
# TRP 39
# TRP 39
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# ALA 21
# LEU 34
# TRP 39
# TRP 39
# ALA 21
# CAP _40_
source configs/peptide2popc.tcl

set resilist [list NAP GLY TRP TRP LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU ALA LEU TRP TRP ALA CAP]
#puts "resilist: $resilist"
set beadlist {}
foreach resname $resilist {
	foreach partinfo_this_resi $partlist_per_res3letter {
	     set this_resi_name [lindex $partinfo_this_resi 0]
	     #puts "this_resi_name: $this_resi_name"
	     #puts "resname: $resname"
	     if { $resname == $this_resi_name } {
		     set beadlist_this_resi [lindex $partinfo_this_resi 1]
		     foreach thispartnum $beadlist_this_resi {
			lappend beadlist $thispartnum
	    	     }
	     }
	}
}
#puts "beadlist: $beadlist"
unset partlist_per_res3letter
# Don't specify topology here
set bondlist [list ]
set angllist [list ]
set dihelist [list ]
lappend respartlist $beadlist
lappend respartlist $bondlist
lappend respartlist $angllist
lappend respartlist $dihelist
lappend respartlist $resilist
unset beadlist
unset bondlist
unset angllist
unset dihelist
unset resilist 

#moltypeid
lappend molpeptidepartlist "1" 
#molspec
lappend molpeptidepartlist "PROT" 
#
lappend molpeptidepartlist $respartlist 
lappend molpeptidepartlist $resparttypelist 
lappend molpeptidepartlist $respartcharmmbeadlist 
#puts "molpeptidepartlist: $molpeptidepartlist"
#exit
unset respartlist
unset resparttypelist
unset respartcharmmbeadlist
::mmsg::send [namespace current] "molpeptidepartlist set"


