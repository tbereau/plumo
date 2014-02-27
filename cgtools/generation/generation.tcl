# cgtools::generation --
#
# Generate molecular topology, parameter and positions
# 
# Author: Zun-Jing Wang
# Sep. 25 2008

package require ::mmsg 
package require ::cgtools::utils

package provide ::cgtools::generation 1.0.0

# Create the empty namespace to which we shall add our routines
namespace eval ::cgtools::generation {
  # Global routine for setting up the system
  namespace export generate_system

  # Global variables for setup routines
  variable topology
  variable boxl
  variable notopo
  # Bookkeeps coordinates read in each file. Useful when reading 
  # the same file for different types of molecules (e.g., proteins AND lipids).
  variable coordinput []

  variable trappedmols
  variable trappedmolsupdated
  variable userfixedparts
  variable interbead

  ### temporary
  variable count_proteins 0
  #####

  # Read in all the routines
  # Files with separate namespaces corresponding to geometries
  source [file join [file dirname [info script]] creat_fromfile.tcl]
  source [file join [file dirname [info script]] creat_random.tcl]
  source [file join [file dirname [info script]] creat_fromvectors.tcl]


  # Helper function files
  source [file join [file dirname [info script]] topologies.tcl]
  source [file join [file dirname [info script]] placemol.tcl]
  source [file join [file dirname [info script]] placeparticles.tcl]
}

# ::cgtools::generation::generate_system -- 
#
# A large scale wrapper routine for all of the system setup commands
# in this module.
# 
# This routine should allow flexible use of multiple geometrical
# structures and topologies in a simulation.  A complete topology and
# geometry specification is provided for each object (eg. singlemol,
  # bilayer etc) in the system and these are grouped into the list
# <system_specs>.  generate_system loops through each of these objects
# setting each up independently and then combining topologies at the
# end.
#
# Arguments:
# 
# system_specs: A list of specified system objects.  Each object is
# itself a list of the form < geometry n_lipidslist beads_per_mol >
# where geometry may be any allowable geometry (eg singlemol, bilayer 
  # etc). <n_lipidslist> is a list containing the number of molecules of
# each of the molecule types in this object and beads_per_mol
# specifies the number of particles in each of the molecule types.
#
# 
proc ::cgtools::generation::generate_system { system_specs iboxl } {
  ::mmsg::send [namespace current] "setting up system "

  # The molecule types spec should be globally accessible
  variable boxl
  variable topology 
  variable notopo

  set boxl $iboxl

  set topolist 0
  unset topolist

  set topologieslist 0
  unset topologieslist

  # Starting value for particle ids
  set startp 0
  # Current particle id
  set currpid $startp


  foreach spec $system_specs {
    variable ::cgtools::moltypelists
    #puts "::cgtools::generation::generate_system::spec = $spec"
    # Flags for tracking input settings
    set geometryset 0
    set n_molslistset 0
    set n_molslist 0
    set sequenceset 0
    set notopo 0
    foreach item $spec {
      switch [lindex $item 0] {
        "geometry" {
          set geometry [lindex $item 1]
          set geometryreadfile [lindex $item 2]
          set geometryset 1
        }
        "n_molslist" {
          set n_molslist [lindex $item 1]
          set n_molslistset 1
        }
        "sequence" {
          set sequence [lindex $item 1]
          set sequenceset 1
        }
        "default" {
          ::mmsg::warn [namespace current] "unknown item [lindex $item 0] in system spec. allowed values are: \n geometry \n n_molslist  "
        }
      }
    }


    if { !$geometryset } {
      mmsg::err [namespace current] "geometry not specified"
    }
    if { !$n_molslistset } {
      mmsg::err [namespace current] "n_molslist not specified for [lindex $geometry 0]"
    }
    if { [lindex [lindex $n_molslist 0] 0] == "1" } {
      if { !$sequenceset } {
        mmsg::err [namespace current] "Missing sequence for peptide"
      } else {
        # Generate peptide topology
        lappend moltypelists [gen_peptide_topol $sequence]
        ::cgtools::utils::initmoltypeskey $moltypelists
      }
    }

    # Generate a topology from a list of the number and size of
    # each molecule
    foreach mtp $n_molslist {
      set thismoltypeid [lindex $mtp 0]
      set nmols [lindex $mtp 1]
      set tpspec [::cgtools::utils::matchtype [lindex $mtp 0]]
      set nbeads_mol [llength [lindex [lindex $tpspec 2] 0]]
      # Create the topology for this lipid type
      set topo [create_simple_topo $nmols $nbeads_mol -moltype  $thismoltypeid -startpart $currpid ]		

      # Just in case zero molecules were specified we need
      # to check if topo was actually created at all by the
      # last command
      if { $topo == 0 } {
        ::mmsg::err [namespace current] "no topo created for molecule type $thismoltypeid"
      } else {
        lappend topolist $topo
        set currpid [expr [::cgtools::utils::maxpartid $topo ] + 1]
      }

    }

    # Join all of the previously made topologies
    set first 1
    foreach topo $topolist {
      if { $first } { 
        set topology $topo 
        set first 0
      } else {
        set topology [join_topos $topology $topo]
      }
    }
    unset topolist

    # Now wrap the topology onto a specified geometry and perform any
    # other geometry specific tasks

    # Shuffle the topology
    set topology [shuffle_topo $topology ]
    # Now run the creation command for the specified geometry
    set createprefix "create_"
    set namespaceprefix "::cgtools::generation::"
    # Construct the name of the create command
    set command $geometry
    set geometry [lindex [split $geometry " "] 0]
    set createcommand "$namespaceprefix$geometry\:\:$createprefix$command "
    ::mmsg::debug [namespace current] "executing $command"
    #eval $createcommand -readfile $geometryreadfile
    #puts "::cgtools::generation::generate_system::geometryreadfile = $geometryreadfile"
    if { [catch  {eval $createcommand -readfile $geometryreadfile} errm ] } {
      mmsg::err [namespace current] "couldn't execute creation command for $command \n $errm"
    }

    if {!$notopo} {
      lappend topologieslist $topology
    }
    #puts "topology: $topology"

  }

  # Join all of the previously made topologies
  set first 1
  foreach topo $topologieslist {
    if { $first } { 
      set topology $topo 
      set first 0
    } else {
      set topology [join_topos $topology $topo]
    }

  }
  
  # HACKED IN. CLEAN UP EVENTUALLY.
  # If we find two PARTs (ID==2) in the second part of the system and variables 
  # $::k_umb and ::$d_umb, we turn on umbrella  between the COM of the bilayer 
  # and each particle automatically.
  if { [llength $topologieslist] == 2 } {
    puts "step1"
    if { [llength [lindex $topologieslist 1]] == 2 } {
      puts "step2"
      if { [lindex [lindex [lindex $topologieslist 1] 0] 0] == 2 &&
	[lindex [lindex [lindex $topologieslist 1] 1] 0] == 2 &&
	[info exists ::k_umb] == 1 && [info exists ::z0_umb_1] == 1 &&
	[info exists ::z0_umb_2] == 1 } {
	# Turn on umbrella between COM of bilayer and each PART.
	set nlipids [llength [lindex $topologieslist 0]]
	set nlipid  [expr [llength [lindex [lindex $topologieslist 0] 0]]-1]
	puts "$nlipids $nlipid"
	set bead1 [expr $nlipids*$nlipid]
	set bead2 [expr $nlipids*$nlipid+1]
	
	analyze set chains 1 1 [expr $nlipids*$nlipid]
	set memcomz [::cgtools::utils::compute_membrane_comz $topology]
	set comzpart [setmd n_part]
	part $comzpart pos 1 1 $memcomz virtual 1 molecule 0 type 9
	
	inter 80 umbrella $::k_umb 2 $::z0_umb_1
	inter 91 umbrella $::k_umb 2 $::z0_umb_2
	part $bead1 bond 80 $comzpart
	part $bead2 bond 81 $comzpart
	exit
      }
    }
  }

  #puts "topology= $topology"
  set topology [sort_topo $topology]
  return $topology

}

proc ::cgtools::generation::gen_peptide_topol { sequence } {
  variable ::cgtools::forcefield::partlist_per_res3letter
  variable ::cgtools::forcefield::resparttypelist
  variable ::cgtools::forcefield::respartcharmmbeadlist

  set beadlist {}
  foreach resname $sequence {
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
  lappend respartlist $sequence
  unset beadlist
  unset bondlist
  unset angllist
  unset dihelist

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
  return $molpeptidepartlist
}

proc ::cgtools::generation::get_trappedmols {  } {
  variable trappedmols
  variable topology

  variable trappedmolsupdated

  if { [catch { set dum $trappedmols } ] } {
   ::mmsg::warn [namespace current] "no trappedmols defined"
   return -1
   } else {
     if { !$trappedmolsupdated } {
       set didntfindmol 1
       for { set j 0 } { $j < [llength $trappedmols] } { incr j } {
        # Update trappedmols
        set fmol [lindex $trappedmols $j]
        for { set i 0 } { $i < [llength $topology] } { incr i } {
          set mol [lindex $topology $i]
          if { [lindex $fmol 0] == [lindex $mol 1] } {
           lset trappedmols $j 0 $i
           set didntfindmol 0
           break			
         }
       }
       if { $didntfindmol } {
        ::mmsg::err [namespace current] "could not get_trappedmols unable to find the corresponding particles"
      }

    }
  }
  set trappedmolsupdated 1
  return $trappedmols
}
}

proc ::cgtools::generation::get_userfixedparts {  } {
  variable userfixedparts

  if { [catch { set dum $userfixedparts } ] } {
   ::mmsg::warn [namespace current] "no user fixed particles defined"
   return -1
   } else {
     return $userfixedparts
   }
 }

 proc ::cgtools::generation::get_interbead {  } {
  variable interbead
  return $interbead
}
