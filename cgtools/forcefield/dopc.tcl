# Zun-Jing Wang 2010 June 23 
#
# dopc.tcl

# --- Specify mole type: DOPC lipid -----#
# The atom types in this lipid are respectively 0:CH, 1:PH, 2:GL, 3:E1, 4:E2, 5:AS, 6:AD, 7:AE. 
# BOND:
# The bonding interactions are respectively 0:CHPH, 1:PHGL, 2:GLE1, 3:GLE2, 4:E1AS, 5:E2AS, 6:ASAS, 7:ASAD, 8:ASAE. 
# The bending interactions are respectively 9:E1GLE2 10:CHPHGL, 11:PHGLE1,  12:PHGLE2, 13:GLE1AS, 14:GLE2AS, 15:E1ASAS, 16:E2ASAS, 17:ASASAS, 18:ASASAD, 19:ASADAS, 20:ASASAE. 
# list of beads:{beadid}
set beadlist [list 0 1 2 3 5 5 6 5 5 7 4 5 5 6 5 5 7 ]
# list of bonds:{bondid {bead1# bead2# in the bead lists}}
set bondlist [list { 0 { 0 1 } } { 1 { 1 2 } } { 2 { 2 3 } }  { 4 { 3 4 } } { 6 { 4 5 } } { 7 { 5 6 }} { 7 { 6 7 } } { 6 { 7  8 } } { 8 { 8 9 }} { 3 { 2 10 } } { 5 { 10 11 } } { 6 { 11 12 } } { 7 { 12 13 } } { 7 { 13 14 } } { 6 { 14 15 } } { 8 { 15 16 } } { 9 { 3 10 } } { 10 { 0 2 } } { 11 { 1 3 } } { 13 { 2 4 } } { 15 { 3 5 } } { 18 { 4 6 } } { 19 { 5 7 } } { 18 { 6 8 } } { 20 { 7 9 } }  { 12 { 1 10 } } { 14 { 2 11 } } { 16 { 10 12 } } { 18 { 11 13 } } { 19 { 12 14 } } { 18 { 13 15 } } { 20 { 14 16 } } ]
# list of angls:{bondid {bead1# bead2# bead3# in the bead lists}}
set angllist [list ]
# list of dihes:{bondid {bead1# bead2# bead3# bead4# in the bead lists}} 
set dihelist [list ]
lappend dopclist $beadlist
lappend dopclist $bondlist
lappend dopclist $angllist
lappend dopclist $dihelist
unset beadlist
unset bondlist
unset angllist
unset dihelist

# list of beads names for each of types:{beadid beadname beadmass beadcharge}
set beadtypelist [list { 0 CH 87.166 1.1} { 1 PH  94.970 -1.2} { 2 GL 41.073 0.36} { 3 E1 58.036 -0.13} { 4 E2 58.036 -0.13 } { 5 AS 42.081 0.0} { 6 AD 26.038 0.0} { 7 AE 29.062 0.0} ]
# list of bond names for each of types:{bondid {bondnames}} 
set bondtypelist [list { 0 CH PH } { 1 PH GL } { 2 GL E1 } { 3 GL E2 } { 4 E1 AS } { 5 E2 AS } { 6 AS AS } { 7 AS AD } { 8 AS AE } { 9 E1 GL E2 } { 10 CH PH GL } { 11 PH GL E1 } { 12 PH GL E2 } { 13 GL E1 AS } { 14 GL E2 AS } { 15  E1 AS AS } { 16 E2 AS AS } { 17 AS AS AS } { 18 AS AS AD } { 19 AS AD AS } {20 AS AS AE } ]
# list of angl names for each of types:{bondid {anglnames}}  
set angltypelist [list ]
# list of dihe names for each of types:{bondid {dihenames}} 
set dihetypelist [list ]
# list of nonbonded interaction names for each of types:{{typeid1 typeid2} {nonbnames}} 
#Notice rename E1 E2 to ES, total 28 types of nonbonded interactions
set nonbtypelist [list { { 0 0 } { CH CH } } { { 0 1 } { CH PH } } { { 0 2 } { CH GL } } { { 0 3 } { CH ES } } { { 0 4 } { CH ES } } { { 0 5 } { CH AS } } { { 0 6 } { CH AD } } { { 0 7 } { CH AE } } { { 1 1 } { PH PH } } { { 1 2 } { PH GL } } { { 1 3 } { PH ES } } { { 1 4 } { PH ES } } { { 1 5 } { PH AS } } { { 1 6 } { PH AD } } { { 1 7 } { PH AE } } { { 2 2 } { GL GL } } { { 2 3 } { GL ES } } { { 2 4 } { GL ES } } { { 2 5 } { GL AS } } { { 2 6 } { GL AD } }  { { 2 7 } { GL AE } } { { 3 3 } { ES ES } } { { 3 4 } { ES ES } } { { 3 5 } { ES AS } } { { 3 6 } { ES AD } } { { 3 7 } { ES AE } } { { 4 4 } { ES ES } } { { 4 5 } { ES AS } } { { 4 6 } { ES AD } } { { 4 7 } { ES AE } } { { 5 5 } { AS AS } } { { 5 6 } { AS AD } } { { 5 7 } { AS AE } } { { 6 6 } { AD AD } } { { 6 7 } { AD AE } } { { 7 7 } { AE AE } } ]
lappend dopctypelist $beadtypelist
lappend dopctypelist $bondtypelist
lappend dopctypelist $angltypelist
lappend dopctypelist $dihetypelist
lappend dopctypelist $nonbtypelist
unset beadtypelist
unset bondtypelist
unset angltypelist
unset dihetypelist
unset nonbtypelist

# list of charmm beads names for each of types:{beadid cahrmmbeadname}
set charmmbeadlist [list { 0 CH } { 1 PH } { 2 GL } { 3 ES1 } { 4 AS11 } { 5 AS12 } { 6 AD13 } { 7 AS14 } { 8 AS15 } { 9 AE16 } { 10 ES2 } { 11 AS21 } { 12 AS22 } { 13 AD23 } { 14 AS24 } { 15 AS25 } { 16 AE26 } ]
lappend dopccharmmbeadlist $charmmbeadlist
unset charmmbeadlist

#moltypeid
lappend moldopclist "0" 
#molspec
lappend moldopclist "DOPC" 
#
lappend moldopclist $dopclist 
lappend moldopclist $dopctypelist 
lappend moldopclist $dopccharmmbeadlist 
unset dopclist
unset dopctypelist
unset dopccharmmbeadlist
lappend moltypelists $moldopclist
::mmsg::send [namespace current] "moldopclist is set"


##########################################################
###############  Potentials and Forces ####################
##########################################################

##########################################################
# ---- Potentials and Forces of LIPIDS only ----#
set partbondnamelist [lindex $moldopclist 3]
set beadnamelist [lindex $partbondnamelist 0 ]
set bondnamelist [lindex $partbondnamelist 1 ]
set anglnamelist [lindex $partbondnamelist 2 ]
set dihenamelist [lindex $partbondnamelist 3 ]
set nonbnamelist [lindex $partbondnamelist 4 ]
unset partbondnamelist

# BOND Potentials
#----------------------------------------------------------#
set nbond [llength  $bondnamelist]
for { set p 0 } { $p < $nbond } { incr p } {
	set bondname [lindex $bondnamelist $p]
	set bondfeature [llength $bondname]
	if {$bondfeature == 3} {
		set bondid [lindex $bondname 0] 
		set bn0 [lindex $bondname 1]
		set bn1 [lindex $bondname 2]
		set overlapnamenow overlap_bond.$bn0$bn1\.coff
		lappend overlapnames $overlapnamenow 
		require_feature OVERLAPPED
		lappend bonded_parms [list $bondid overlapped bond $::cgtools::overlapdir/$overlapnamenow]
	}
	if {$bondfeature == 4} {
		set bondid [lindex $bondname 0] 
		set bn0 [lindex $bondname 1]
		set bn1 [lindex $bondname 2]
		set bn2 [lindex $bondname 3]
		set overlapnamenow overlap_bend.$bn0$bn1$bn2\.coff
		lappend overlapnames $overlapnamenow 
		require_feature OVERLAPPED
		lappend bonded_parms [list $bondid overlapped bond $::cgtools::overlapdir/$overlapnamenow]
	}
}
unset nbond bondname bondid bn0 bn1 overlapnamenow
unset bondnamelist
#----------------------------------------------------------#

# ANGL Potentials
#----------------------------------------------------------#
set nangl [llength  $anglnamelist]
for { set p 0 } { $p < $nangl } { incr p } {
	set anglname [lindex $anglnamelist $p]
	set bondid [lindex $anglname 0] 
	set an0 [lindex $anglname 1]
	set an1 [lindex $anglname 2]
	set an2 [lindex $anglname 3]
	set overlapnamenow overlap_cosangl.$an0$an1$an2\.coff
	lappend overlapnames $overlapnamenow 
	require_feature OVERLAPPED
	lappend bonded_parms [list $bondid overlapped angle $::cgtools::overlapdir/$overlapnamenow]
}


if {$nangl > 0} {
	unset nangl anglname bondid an0 an1 an2 overlapnamenow
	unset anglnamelist
}
#----------------------------------------------------------#

# DIHE Potentials
#----------------------------------------------------------#
set ndihe [llength  $dihenamelist]
for { set p 0 } { $p < $ndihe } { incr p } {
	set dihename [lindex $dihenamelist $p]
	set bondid [lindex $dihename 0] 
	set dn0 [lindex $dihename 1]
	set dn1 [lindex $dihename 2]
	set dn2 [lindex $dihename 3]
	set dn3 [lindex $dihename 4]
	set overlapnamenow overlap_dihe.$dn0$dn1$dn2$dn3\.coff
	lappend overlapnames $overlapnamenow 
	require_feature OVERLAPPED
	lappend bonded_parms [list $bondid overlapped dihedral $::cgtools::overlapdir/$overlapnamenow]
}
if {$ndihe > 0} {
        unset ndihe dihename bondid dn0 dn1 dn2 dn3 overlapnamenow
        unset dihenamelist
}
#----------------------------------------------------------#


# Non Bonded Potentials 
#----------------------------------------------------------#
# Define the interactions between lipid beads
# NO Non-bonded interactions for this one molecular system

set nnonb [llength $nonbnamelist]
for { set i 0 } { $i < $nnonb } { incr i } {
	set nonbinfo [lindex $nonbnamelist $i]
	set btype [lindex $nonbinfo 0]
	set typei [lindex $btype 0]
	set typej [lindex $btype 1]
	set tabname [lindex $nonbinfo 1]
	set tabnamei [lindex $tabname 0]
	set tabnamej [lindex $tabname 1]
	set tablenamenow pot_force_nonbond.$tabnamei$tabnamej\.tab
	lappend tablenames $tablenamenow 
	require_feature TABULATED
	lappend nb_interactions [list $typei $typej tabulated $::cgtools::tabledir/$tablenamenow]
}
unset nnonb nonbinfo btype typei typej tabname tabnamei tabnamej tablenamenow 
unset beadnamelist  
unset nonbnamelist  
::mmsg::send [namespace current] "overlapnames, tablenames AND bonded_parms, nb_interactions are set"


