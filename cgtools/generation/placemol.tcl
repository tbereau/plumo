#
# Routines for placing Positions, Bonds, Angles, Dihedrals of a molecule 
#
# Author: Zun-Jing Wang
# Sep.26 - Oct 1 2008 Done 


namespace eval cgtools::generation {}

# ::cgtools::generation::placemol-- 
# general routine for placing molecules
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
proc ::cgtools::generation::placemol { mol partlist args } { 

    set options {
        {changepos.arg "0" "if change positions the particles"  }
        {move.arg   "{ 0. 0. 0. }" "movement vector of the center of mass of mol"  }
        {rotate.arg  "{ { 1. 0. 0.} { 0. 1. 0.} { 0. 0. 1.} }" "rotation matrix for the mol " }
    }
    set usage "Usage: create_bilayer topo boxl \[changepos:move:rotate]"
    array set params [::cmdline::getoptions args $options $usage]
    
    # Retrieve the molecule information for this molecule type  
    set typekey [::cgtools::utils::matchtype [lindex $mol 0]]

    # Place the molecule depending on type
    switch [lindex $typekey 1] {
        "POPC" {
            if {$params(changepos)==0} {
                place_lipid $mol $partlist
            }
            if {$params(changepos)!=0} {
                set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
                place_lipid $mol $partlist_new 
            }
        }
        "DOPC" {
            if {$params(changepos)==0} {
                place_lipid $mol $partlist
            }
            if {$params(changepos)!=0} {
                set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
                place_lipid $mol $partlist_new 
            }
        }
        "DPPC" {
            if {$params(changepos)==0} {
                place_lipid $mol $partlist
            }
            if {$params(changepos)!=0} {
                set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
                place_lipid $mol $partlist_new 
            }
        }
        "PART" {
            if {$params(changepos)==0} {
                place_part $mol $partlist
            }
            if {$params(changepos)!=0} {
                set partlist_new [place_posmoverotation $partlist $params(move) $params(rotate)]
                place_part $mol $partlist_new 
            }
        }
        "PROT" {
            place_protein $mol $partlist
        }
        "default" {
            ::mmsg::err [namespace current] "couldn't place molecule of type [lindex $typekey 1], possibilities are: \n lipid \n hollowsphere \n sphericalconstraint"
        }
    }

    return
}

# ::cgtools::generation::place_lipid-- 
# Place positions of particles and bonds between particles of a lipid molecule 
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
# Note that both linkbond and bendbond will have previously been set by set_bonded_interactions
proc ::cgtools::generation::place_lipid { mol partlist } {

    set moltype [lindex $mol 0]
    set typeinfo [::cgtools::utils::matchtype $moltype]
    
    set partbondlists [lindex $typeinfo 2]
    set partbondtypelists [lindex $typeinfo 3]

    #particle positions
    set beadlists [lindex $partbondlists 0]
    set nbeads [llength $beadlists]

    set beadtypelists [lindex $partbondtypelists 0]

    for { set b 0 } { $b < $nbeads } {incr b } {

        #current positions of particles
        set partnum [lindex $mol [ expr $b + 1] ]
        set parttype [lindex $beadlists $b]
        set parttypeinfo [lindex $beadtypelists $parttype]
        set partmass [lindex $parttypeinfo 2]

        set curpart [lindex $partlist $b]
        set curpos [lindex $curpart 1]
        set posx [lindex $curpos 0]
        set posy [lindex $curpos 1]
        set posz [lindex $curpos 2]
        
        part $partnum pos $posx $posy $posz type $parttype mass $partmass 
    }

    #bonds
    set bondlists [lindex $partbondlists 1]
    set nbonds [llength $bondlists]
    for { set b 0 } { $b < $nbonds } {incr b } {
        set curbond [lindex $bondlists $b ]
        # bond type i.e. bondid 
        set btype [lindex $curbond 0]
        # index of the particles inside mol 
        set partlists_inmol [lindex $curbond 1]
        set part1_inmol [lindex $partlists_inmol 0]
        set part2_inmol [lindex $partlists_inmol 1]
        # pid [0:npart-1] to link the bond 
        set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
        set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
        
        part $partnum2 bond $btype $partnum1
        #puts "$btype $partnum1 $partnum2"
    }
    
    #angl 
    set angllists [lindex $partbondlists 2]
    set nangls [llength $angllists]
    for { set b 0 } { $b < $nangls } {incr b } {
        set curbond [lindex $angllists $b ]
        # bond type i.e. bondid 
        set btype [lindex $curbond 0]
        # index of the particles inside mol 
        set partlists_inmol [lindex $curbond 1]
        set part1_inmol [lindex $partlists_inmol 0]
        set part2_inmol [lindex $partlists_inmol 1]
        set part3_inmol [lindex $partlists_inmol 2]
        # pid [0:npart-1] to link the angle 
        set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
        set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
        set partnum3 [lindex $mol [ expr $part3_inmol + 1] ]
        
        part $partnum2 bond $btype $partnum1 $partnum3
        #puts "$btype $partnum1 $partnum2 $partnum3"
    }
    
    #dihe 
    set dihelists [lindex $partbondlists 3]
    set ndihes [llength $dihelists]
    for { set b 0 } { $b < $ndihes } {incr b } {
        set curbond [lindex $dihelists $b ]
        # bond type i.e. bondid 
        set btype [lindex $curbond 0]
        # index of the particles inside mol 
        set partlists_inmol [lindex $curbond 1]
        set part1_inmol [lindex $partlists_inmol 0]
        set part2_inmol [lindex $partlists_inmol 1]
        set part3_inmol [lindex $partlists_inmol 2]
        set part4_inmol [lindex $partlists_inmol 3]
        # pid [0:npart-1] to link the dihedral
        set partnum1 [lindex $mol [ expr $part1_inmol + 1] ]
        set partnum2 [lindex $mol [ expr $part2_inmol + 1] ]
        set partnum3 [lindex $mol [ expr $part3_inmol + 1] ]
        set partnum4 [lindex $mol [ expr $part4_inmol + 1] ]
        
        part $partnum2 bond $btype $partnum1 $partnum3 $partnum4
        #puts "$btype $partnum1 $partnum2 $partnum3 $partnum4"
    }
}

proc ::cgtools::generation::place_posmoverotation {partlist move rotate} {
    #::mmsg::send [namespace current] "entering place_posmoverotation"
    
    set npart [llength $partlist]

    set partlist_new 0  
    unset partlist_new
    for { set i 0 } { $i < $npart } {incr i } {
        set pos_old [lindex [lindex $partlist $i] 1]
        set pos_move [::cgtools::utils::add_vecs $pos_old $move]
        set pos_rotate [::cgtools::utils::matrix_vec_multiply $rotate $pos_move]

        set ninfo [llength [lindex $partlist $i]]
        set poslist_new 0   
        unset poslist_new
        for { set j 0 } { $j < $ninfo } {incr j } {
            if { $j!=1 } {
                lappend poslist_new [lindex [lindex $partlist $i] $j]
            }
            if { $j==1 } {
                lappend poslist_new $pos_rotate
            }
        }
        lappend partlist_new $poslist_new
    }

    return $partlist_new
}

# ::cgtools::generation::place_part-- 
# Place positions of particles for sigle-part mol
# Arguments:
# mol : list of atoms in the molecule
# partlist : list of particle information of the atoms 
# Note that both linkbond and bendbond will have previously been set by set_bonded_interactions
proc ::cgtools::generation::place_part { mol partlist } {
    set moltype [lindex $mol 0]
    set typeinfo [::cgtools::utils::matchtype $moltype]
    #puts "typeinfo= $typeinfo"
    
    set partbondlists [lindex $typeinfo 2]
    set partbondtypelists [lindex $typeinfo 3]

    #particle positions
    set beadlists [lindex $partbondlists 0]
    set nbeads [llength $beadlists]
    set beadtypelists [lindex $partbondtypelists 0]
    set itype_begin [lindex [lindex $beadtypelists 0] 0]
    #puts "itype_begin= $itype_begin"
    #puts "beadlists= $beadlists"
    #puts "nbeads= $nbeads"
    #puts "beadtypelists= $beadtypelists"

    for { set b 0 } { $b < $nbeads } {incr b } {

        #current positions of particles
        set partnum [lindex $mol [ expr $b + 1] ]
        set parttype [lindex $beadlists $b]
        #puts "parttype= $parttype"
        set parttypeinfo [lindex $beadtypelists [expr $parttype - $itype_begin]]
        set partmass [lindex $parttypeinfo 2]

        set curpart [lindex $partlist $b]
        set curpos [lindex $curpart 1]
        set posx [lindex $curpos 0]
        set posy [lindex $curpos 1]
        set posz [lindex $curpos 2]
        
        part $partnum pos $posx $posy $posz type $parttype mass $partmass
        #puts "part $partnum pos $posx $posy $posz type $parttype mass $partmass"
        #exit
    }
}

# Place particles that belong to a protein
proc ::cgtools::generation::place_protein { mol partlist } {
    variable count_proteins
    set this [namespace current]

    incr count_proteins

    set moltype [lindex $mol 0]
    set typeinfo [::cgtools::utils::matchtype $moltype]
    
    set partbondlists [lindex $typeinfo 2]
    set partbondtypelists [lindex $typeinfo 3]

    #particle positions
    set beadlists [lindex $partbondlists 0]
    set nbeads [llength $beadlists]

    set beadtypelists [lindex $partbondtypelists 0]
    if { [expr [llength $mol]-1] != [llength $partlist] } {
        ::mmsg::err $this "Incompatibility between input file and\
            topology for protein (ID: [lindex $mol 0])."
    }
    set resilists [lindex $partbondlists 4]
    for { set b 0 } { $b < $nbeads } {incr b } {

        #current positions of particles
        set partnum [lindex $mol [ expr $b + 1] ]
        set parttype [lindex $beadlists $b]

        # Loop over all bead types to find the correct one
        set parttypeindex -1
        for { set i 0 } { $i < [llength $beadtypelists] } { incr i } {
            if { [lindex [lindex $beadtypelists $i] 0] == $parttype } {
                set parttypeindex $i
                break
            }
        }
        if { $parttypeindex == -1} {
            ::mmsg::err $this "Can't find particle type $parttype"
        }

        set parttypeinfo [lindex $beadtypelists $parttypeindex]
        set partmass [lindex $parttypeinfo 2]

        set curpart [lindex $partlist $b]
        set curpos [lindex $curpart 1]
        set posx [lindex $curpos 0]
        set posy [lindex $curpos 1]
        set posz [lindex $curpos 2]

        set mts ""
        variable ::cgtools::multitimestep
        if { $multitimestep > 0} {
            set mts " smaller_timestep $multitimestep"
        }
        set part_cmd "part $partnum pos $posx $posy $posz type $parttype \
            mass $partmass fix 1 1 1 mol $count_proteins"
        append part_cmd $mts
        eval $part_cmd
    }
    ::mmsg::send $this "Assigned positions to protein (ID: [lindex $mol 0])"

    #bonds
    set index_res 0 
    for { set b 0 } { $b < $nbeads } {incr b 5 } {
        set resname_now [lindex $resilists $index_res]
    #puts "resname_now: $resname_now"
    set partnum_N [lindex $mol [expr $b+1]]
    set partnum_CA [lindex $mol [expr $b+2]]
    set partnum_CB [lindex $mol [expr $b+3]]
    set partnum_C [lindex $mol [expr $b+4]]
    set partnum_O [lindex $mol [expr $b+5]]

    if { $b > 0 } {
        # angle _C-N-CA
        part $partnum_N bond 108 [expr $partnum_N-2] $partnum_CA 
            ## set auto_exclusion --- membrane: 1, peptide: 2  --- _C-N-CA add exclusion between _C and CA
        part [expr $partnum_N-2] exclude $partnum_CA
        # dihedral phi
        part $partnum_N bond 111 [expr $partnum_N-2] $partnum_CA $partnum_C
    }


    # bonds
    part $partnum_N bond 100 $partnum_CA
        if {$resname_now == "GLY"} {
            part $partnum_CA bond 115 $partnum_CB
        } else {
            part $partnum_CA bond 116 $partnum_CB
        } 
    part $partnum_CA bond 102 $partnum_C
    part $partnum_C bond 117 $partnum_O
    ### bond 120: virtual bond (because auto_exclusions=1...) --- deleted and replaced with exclude
    ##part $partnum_N bond 120 $partnum_CB
    ##part $partnum_N bond 120 $partnum_C
    ##part $partnum_CB bond 120 $partnum_C
    # angles
    part $partnum_CA bond 104 $partnum_N $partnum_CB
        ## set auto_exclusion --- membrane: 1, peptide: 2  --- N-CA-CB add exclusion between N and CB
    part $partnum_N exclude $partnum_CB
    part $partnum_CA bond 105 $partnum_N $partnum_C 
        ## set auto_exclusion --- membrane: 1, peptide: 2  --- N-CA-C add exclusion between N and C
    part $partnum_N exclude $partnum_C
    part $partnum_CA bond 106 $partnum_CB $partnum_C 
        ## set auto_exclusion --- membrane: 1, peptide: 2  --- CB-CA-C add exclusion between CB and C
    part $partnum_CB exclude $partnum_C
    part $partnum_C  bond 118 $partnum_CA $partnum_O 
        ## set auto_exclusion --- membrane: 1, peptide: 2  --- CA-C=O add exclusion between CA and O 
    part $partnum_CA exclude $partnum_O
    # dihedrals
    part $partnum_CA bond 112 $partnum_N $partnum_C $partnum_CB

    if { $b < [expr $nbeads-5]} {
        # bond
        part $partnum_C bond 103 [expr $partnum_C+2]
        ### bond 120: virtual bond (because auto_exclusions=1...)--- deleted and replaced with exclude      
        ##part $partnum_CA bond 120 [expr $partnum_O+1]
        ##part $partnum_CB bond 120 [expr $partnum_O+1]
        ##part $partnum_C  bond 120 [expr $partnum_O+2]
        ##part $partnum_C  bond 120 [expr $partnum_O+3]
        # angle
        part $partnum_C bond 107 $partnum_CA [expr $partnum_C+2] 
            ## set auto_exclusion --- membrane: 1, peptide: 2  --- CA-C-N+ add exclusion between CA and N+ 
        part $partnum_CA exclude [expr $partnum_C+2]
            ## NOT SURE IF NEEDED set auto_exclusion --- membrane: 1, peptide: 2  --- O=C-N+ no angle, but add exclusion between O and N+
        part $partnum_O exclude [expr $partnum_C+2]
        # dihedral psi CA*-N-CB-C
        part $partnum_CA bond 111 $partnum_N $partnum_C [expr $partnum_C+2] 
        # dihedral omega C*-CA-O-N
            if {$resname_now == "PRO"} {
            part $partnum_C bond 114 $partnum_CA [expr $partnum_C+2] [expr $partnum_C+3]
            } else {
            part $partnum_C bond 110 $partnum_CA [expr $partnum_C+2] [expr $partnum_C+3]
            }
        # improper dihedral for C*-CA-O-N+
        part $partnum_C bond 119 $partnum_CA $partnum_O [expr $partnum_O+1]
    }
    incr index_res
    }
    ::mmsg::send $this "Assigned bonds, angles, and dihedrals to protein (ID: [lindex $mol 0])."

    # Now turn on Hbond interaction (after bonded partners are defined)
    ::cgtools::forcefield::source_hbond_ff
}
