# Zun-Jing Wang 
# 2011 August    

# Tristan Bereau
# Spring 2013


# Check that all interactions are compiled in
require_feature LENNARD_JONES
require_feature LENNARD_JONES_GENERIC

# The atom types in this peptide are respectively 8:N, 9:CA, 10: Pro-N (no H-bond), 11:C, 12:O, 15:term-N, 16:term-C
# 20-39: CB (Side chain beads)
#
############
# WARNING: ionizable residues are parametrized in their charged form.
#   *** END CAPS ADDED ***
###########
#

# list of beads names for each of types:{beadid beadname beadmass beadcharge}
# Type    Atom
#  8       N
#  9       Ca
#  10      N-proline
#  11      C
#  12      O
#  13-14   Empty
#  15      term-N
#  16      term-C
#  17-19   Empty
#  20-43   Cb
##    20   gly Cb
##    21   ala Cb
##    22   pro Cb
##    23   glu Cb
##    24   gln Cb
##    25   asp Cb
##    26   asn Cb
##    27   ser Cb
##    28   his Cb
##    29   lys Cb
##    30   arg Cb
##    31   thr Cb
##    32   val Cb
##    33   ile Cb
##    34   leu Cb
##    35   met Cb
##    36   phe Cb
##    37   tyr Cb
##    38   cys Cb
##    39   trp Cb
##    40   end Cb

set beadtypelist [list \
  { 8  N  14.0 0} \
  { 9 CA  12.0 0} \
  {10 NP  14.0 0} \
  {11  C  12.0 0} \
  {12  O  16.0 0} \
  {13 XX   1.0 0} \
  {14 XX   1.0 0} \
  {15  N  14.0 0} \
  {16  C  12.0 0} \
  {17 XX   1.0 0} \
  {18 XX   1.0 0} \
  {19 XX   1.0 0} \
  {20 CB  75.0 0 GLY} \
  {21 CB  89.0 0 ALA} \
  {22 CB 115.0 0 PRO} \
  {23 CB 147.0 0 GLU} \
  {24 CB 146.0 0 GLN} \
  {25 CB 133.0 0 ASP} \
  {26 CB 132.0 0 ASN} \
  {27 CB 105.0 0 SER} \
  {28 CB 155.0 0 HIS} \
  {29 CB 146.0 0 LYS} \
  {30 CB 174.0 0 ARG} \
  {31 CB 119.0 0 THR} \
  {32 CB 117.0 0 VAL} \
  {33 CB 131.0 0 ILE} \
  {34 CB 131.0 0 LEU} \
  {35 CB 149.0 0 MET} \
  {36 CB 165.0 0 PHE} \
  {37 CB 181.0 0 TYR} \
  {38 CB 121.0 0 CYS} \
  {39 CB 204.0 0 TRP} \
  {40 CB  75.0 0 CAP}]
##Notice, bond,angle,dihedral,hydrogen bonds are all set in cgtools/generation/placemol.tcl "proc ::cgtools::generation::place_protein"
set bondtypelist [list ]
set angltypelist [list ]
set dihetypelist [list ]
set nonbtypelist [list ]
set partlist_per_res3letter [list \
  {GLY { 8 9 20 11 12}}\
  {ALA { 8 9 21 11 12}}\
  {PRO {10 9 22 11 12}}\
  {GLU { 8 9 23 11 12}}\
  {GLN { 8 9 24 11 12}}\
  {ASP { 8 9 25 11 12}}\
  {ASN { 8 9 26 11 12}}\
  {SER { 8 9 27 11 12}}\
  {HIS { 8 9 28 11 12}}\
  {LYS { 8 9 29 11 12}}\
  {ARG { 8 9 30 11 12}}\
  {THR { 8 9 31 11 12}}\
  {VAL { 8 9 32 11 12}}\
  {ILE { 8 9 33 11 12}}\
  {LEU { 8 9 34 11 12}}\
  {MET { 8 9 35 11 12}}\
  {PHE { 8 9 36 11 12}}\
  {TYR { 8 9 37 11 12}}\
  {CYS { 8 9 38 11 12}}\
  {TRP { 8 9 39 11 12}}\
  {NAP {15 9 40 11 12}}\
  {CAP { 8 9 40 16 12}}]

lappend resparttypelist $beadtypelist
lappend resparttypelist $bondtypelist
lappend resparttypelist $angltypelist
lappend resparttypelist $dihetypelist
lappend resparttypelist $nonbtypelist
lappend resparttypelist $partlist_per_res3letter
unset beadtypelist
unset bondtypelist
unset angltypelist
unset dihetypelist
unset nonbtypelist
#puts "partlist_per_res3letter: $partlist_per_res3letter"

# list of charmm beads names for each of types:{beadid cahrmmbeadname}
set charmmbeadlist [list {8 N } {9 CA } {10 N } {11 C } {12 O } \
      {13 XX } {14 XX }  {15 N } {16 C }\
      {17 XX } {18 XX } {19 XX}\
      {20 CB } {21 CB } {22 CB } {23 CB }\
      {24 CB } {25 CB } {26 CB } {27 CB }\
      {28 CB } {29 CB } {30 CB } {31 CB }\
      {32 CB } {33 CB } {34 CB } {35 CB }\
      {36 CB } {37 CB } {38 CB } {39 CB }\
      {40 CB }]
lappend respartcharmmbeadlist $charmmbeadlist
unset charmmbeadlist

# Source interaction parameters if needed
if { [llength [inter]] == 0 } {
  source [file join [file dirname [info script]] peptide_parameters.tcl    ]
  source [file join [file dirname [info script]] peptide_sc_parameters.tcl ]

  # Check for LJGEN_SOFTCORE feature if it's turned on in the simulation
  if { ![info exists peptideb::softcore_flag] } {
    set peptideb::softcore_flag 0
  }
  if { $peptideb::softcore_flag != 0 } {
    require_feature LJGEN_SOFTCORE
  }
}

# set parameters
# Bonded interactions

# Pairs
# bond 0 is between N and Ca
lappend bonded_parms [list 100 harmonic $peptideb::k_bond $peptideb::bond_NCa]
# bonds 15 for GLY and 16 for non-GLY - between Ca and Cb
lappend bonded_parms [list 115 harmonic $peptideb::k_bond $peptideb::bondCaCb_GLY]
lappend bonded_parms [list 116 harmonic $peptideb::k_bond $peptideb::bondCaCb_XXX]
# bond 2 is between Ca and C
lappend bonded_parms [list 102 harmonic $peptideb::k_bond $peptideb::bond_CaC]
# bond 3 is between C and N
lappend bonded_parms [list 103 harmonic $peptideb::k_bond $peptideb::bond_CN]
# bond 17 is between C and O
lappend bonded_parms [list 117 harmonic $peptideb::k_bond $peptideb::bond_CO]

# Angles

# Angle  N - Ca - Cb
lappend bonded_parms [list 104 angle $peptideb::k_angle $peptideb::angleR_NCaCb]
# Angle  N - Ca - C
lappend bonded_parms [list 105 angle $peptideb::k_angle $peptideb::angleR_NCaC ]
# Angle Cb - Ca - C
lappend bonded_parms [list 106 angle $peptideb::k_angle $peptideb::angleR_CbCaC]
# Angle Ca - C - N
lappend bonded_parms [list 107 angle $peptideb::k_angle $peptideb::angleR_CaCN ]
# Angle C - N - Ca
lappend bonded_parms [list 108 angle $peptideb::k_angle $peptideb::angleR_CNCa ]
# Angle Ca - C - O
lappend bonded_parms [list 118 angle $peptideb::k_angle $peptideb::angleR_CaCO ]

# Dihedrals

# psi and phi now include the dipolar interaction, that biases
# the potential for beta-sheets rather than alpha helices.
# Dihedral \omega
lappend bonded_parms [list 110 dihedral 1 $peptideb::k_dih_omega $peptideb::pi]
# Dihedral \omega for proline
lappend bonded_parms [list 114 dihedral 2 $peptideb::k_dih_w_pro  0.]
# Dipolar interactions
lappend bonded_parms [list 111 dihedral 1 $peptideb::k_dipolar    0.]
# Improper dihedral - for chirality
lappend bonded_parms [list 112 dihedral 1 $peptideb::k_imp_dih   $peptideb::angleR_Im_Cb]
# Dihedral for Oxygen
lappend bonded_parms [list 119 dihedral 1 $peptideb::k_dih_omega $peptideb::pi]

# Subtract Lennard-Jones for bonded partners    
lappend bonded_parms [list 113 subt_lj 0 $peptideb::max_length]

# Virtual bond: deleted
lappend bonded_parms [list 120 virtual_bond]

# Nonbonded interactions
# Type    Atom
#  8       N
#  10      N-proline
#  9       Ca
#  11      C
#  12      O
#  15      term-N
#  16      term-C
#  20-43   Cb

# interaction N  and N
lappend nb_interactions \
    [list 8 8 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_N]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_N ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction N and Ca
lappend nb_interactions \
    [list 8 9 lennard-jones $peptideb::lj_eps $peptideb::lb_NCa  [expr $peptideb::cut_factor*$peptideb::lb_NCa ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction N and C
lappend nb_interactions \
    [list 8 11 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction N and Pro-N
lappend nb_interactions \
    [list 8 10 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_N]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_N ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction Ca and Ca
lappend nb_interactions \
    [list 9 9 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_Ca] [expr $peptideb::cut_factor*2.*$peptideb::rvdw_Ca] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction Ca and C
lappend nb_interactions \
    [list 9 11 lennard-jones $peptideb::lj_eps $peptideb::lb_CaC  [expr $peptideb::cut_factor*$peptideb::lb_CaC ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction C  and C
lappend nb_interactions \
    [list 11 11 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_C]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_C ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction Pro-N and Pro-N
lappend nb_interactions \
    [list 10 10 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_N]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_N ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction Pro-N and Ca
lappend nb_interactions \
    [list 10 9 lennard-jones $peptideb::lj_eps $peptideb::lb_NCa  [expr $peptideb::cut_factor*$peptideb::lb_NCa ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction Pro-N and C
lappend nb_interactions \
    [list 10 11 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-N and N
lappend nb_interactions \
    [list 15 8 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_N]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_N ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-N and Ca
lappend nb_interactions \
    [list 15 9 lennard-jones $peptideb::lj_eps $peptideb::lb_NCa  [expr $peptideb::cut_factor*$peptideb::lb_NCa ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-N and C
lappend nb_interactions \
    [list 15 11 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-N and Pro-N
lappend nb_interactions \
    [list 15 10 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_N]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_N ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-N and term-C
lappend nb_interactions \
    [list 15 16 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-C and N
lappend nb_interactions \
    [list 8 16 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-C and Ca
lappend nb_interactions \
    [list 9 16 lennard-jones $peptideb::lj_eps $peptideb::lb_CaC  [expr $peptideb::cut_factor*$peptideb::lb_CaC ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-C and C
lappend nb_interactions \
    [list 11 16 lennard-jones $peptideb::lj_eps [expr 2.*$peptideb::rvdw_C]  [expr $peptideb::cut_factor*2.*$peptideb::rvdw_C ] \
   $peptideb::lj_shift $peptideb::ljoffset]
# interaction term-C and Pro-N
lappend nb_interactions \
    [list 10 16 lennard-jones $peptideb::lj_eps $peptideb::lb_NC   [expr $peptideb::cut_factor*$peptideb::lb_NC  ] \
   $peptideb::lj_shift $peptideb::ljoffset]

### there are no nonbonded interactions between O and "X": N, CA, CB, C (4-beads model)


# Side chains... part ID from 20 to 39.
set index_sc 20
foreach name $peptideb::3letter_list {
    # No interaction with Gly
    if {$name != "GLY"} {
  
  # interaction N and Sc
  set sigma [expr $peptideb::rvdw_N+$peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 8 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor *$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  # interaction Ca and Sc
  set sigma [expr $peptideb::rvdw_Ca+$peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 9 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor*$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  # interaction C and Sc
  set sigma [expr $peptideb::rvdw_C+$peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 11 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor*$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  # interaction Pro-N and Sc
  set sigma [expr $peptideb::rvdw_N+$peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 10 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor*$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  # interaction term-N and Sc
  set sigma [expr $peptideb::rvdw_N + $peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 15 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor*$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  # interaction term-C and Sc
  set sigma [expr $peptideb::rvdw_C + $peptideb::rvdw_XXX]
  lappend nb_interactions \
      [list 16 $index_sc lennard-jones $peptideb::lj_eps $sigma [expr $peptideb::cut_factor*$sigma] $peptideb::lj_shift $peptideb::ljoffset]
  

  # interaction Sc and Sc - Matrix of coefficients
  set index2_sc 20
  foreach partner $peptideb::3letter_list {
      if {$partner != "GLY"} {

    set sigma [expr $peptideb::rvdw_XXX+$peptideb::rvdw_XXX]
    # epsilon is the geometric mean of the two side chain energies
    set epsilon [expr sqrt([set peptideb::hp_$name] * [set peptideb::hp_$partner])]
    # Multiply epsilon by the side chain coupling (for HREMD...)
    set epsilon [expr $epsilon * $peptideb::hp_coupling]
    # WCA shift
    set eps_rel [expr .25*(1-$epsilon/$peptideb::lj_hp)]
    # WCA - repulsive part (use lj-gen because we can't have lennard-jones twice!)
    lappend nb_interactions \
        [list $index_sc $index2_sc lj-gen $peptideb::lj_hp $sigma [expr $peptideb::cut_factor*$sigma] \
       $eps_rel $peptideb::ljoffset 12 6 1.0 1.0]
    # LJ - attractive part
    lappend nb_interactions \
        [list $index_sc $index2_sc lennard-jones $epsilon $sigma $peptideb::ljhp_cut 0.00 \
       $peptideb::ljoffset 0. [expr $peptideb::cut_factor*$sigma]]
      }
      incr index2_sc
  }
    }
    incr index_sc
}
unset index_sc
unset index2_sc
unset sigma
unset epsilon


# Softcore potentials (these initial values will not affect the original
# potentials)
set lambda    1.0
set lambdaLJ  1.0
set lambdaWCA 1.0
set delta     0.0
if { $peptideb::softcore_flag != 0 } {
  # Soft-core potential is turned on
  set lambda    $peptideb::lambda_coupling
  set lambdaLJ  [::cgtools::utils::max [expr 2*($lambda-0.5)] 0.0]
  set lambdaWCA [::cgtools::utils::min [expr 2* $lambda     ] 1.0]
  set delta     $peptideb::softcore_delta  
}
set sigmaScale 1.0
if { $cgtools::hremd == 1 } {
  # HREMD: scale Hbond strength and sigma of peptide-lipid interactions
  set sigmaScale $peptideb::lambda_coupling
}

### Peptide-lipid interaction ###
# interaction N, Ca, C' with lipid bead types
# lipid bead types: {0-7}
# assume lipid bead has radius = 3
for { set lipid_i 0} { $lipid_i < 8 } { incr lipid_i } {
    set wca_eps   $peptideb::lj_eps
    set wca_sig   [expr ($peptideb::rvdw_Ca + 3.) * $sigmaScale]
    set wca_cut   [expr $wca_sig * sqrt(pow(2,1/3.) \
                            -(1-$lambdaWCA)*$delta)]
    set wca_shift 0.25
    set wca_off   0.0
    set wca_cap   0.0
    set wca_soft  ""
    set wca_command "$lipid_i 9 lj-gen \
      $wca_eps $wca_sig \
      $wca_cut $wca_shift $wca_off \
      12 6 1.0 1.0"
    if { $peptideb::softcore_flag != 0 } {
      set wca_soft " 1.0 $lambdaWCA $delta"
      append wca_command $wca_soft

    }
    lappend nb_interactions $wca_command
    
    set wca_sig [expr ($peptideb::rvdw_N + 3.) * $sigmaScale]
    set wca_cut   [expr $wca_sig * sqrt(pow(2,1/3.) \
                            -(1-$lambdaWCA)*$delta)]
    set wca_command "$lipid_i 8 lj-gen \
      $wca_eps $wca_sig \
      $wca_cut $wca_shift $wca_off \
      12 6 1.0 1.0"
    if { $peptideb::softcore_flag != 0 } {
      set wca_soft " 1.0 $lambdaWCA $delta"
      append wca_command $wca_soft

    }
    lappend nb_interactions $wca_command

    set wca_sig [expr ($peptideb::rvdw_C + 3.) * $sigmaScale]
    set wca_cut [expr $wca_sig * sqrt(pow(2,1/3.) \
                            -(1-$lambdaWCA)*$delta)]
     set wca_command "$lipid_i 11 lj-gen \
      $wca_eps $wca_sig \
      $wca_cut $wca_shift $wca_off \
      12 6 1.0 1.0"
    if { $peptideb::softcore_flag != 0 } {
      set wca_soft " 1.0 $lambdaWCA $delta"
      append wca_command $wca_soft

    }
    lappend nb_interactions $wca_command
}

# input amino acid number
# output interaction with all lipid bead types. Format:
# {
#   { CH_type CH_eps CH_sigma } 
#   { PH_type PH_eps PH_sigma } 
#   { GL_type GL_eps GL_sigma } 
#   { ES_type ES_eps ES_sigma } 
#   { AS_type AS_eps AS_sigma }
#   { AE_type AE_eps AE_sigma }
# }
############
# WARNING: ionizable residues are parametrized in their charged form.
###########
proc aa_lipid_inter { aa_num } {
    switch -glob $aa_num \
  21 { # ala
      return [list { "lj" 1.0 5.96 } { "wca" 1.0 9.54 } \
      { "lj" 4.5 5.26 } { "wca" 1.0 4.61 } \
      { "lj" 0.85 4.21 } { "lj" 0.85 4.21 }]
  } 30 { # arg
      return [list { "lj" 1.2 6.66 } { "wca" 1.0 6.66 } \
          { "lj" 6.5 4.17 } { "lj" 2.0 4.52 } \
          { "wca" 1.0 4.17 } { "wca" 1.0 13.71 }]
  } 26 { # asn
      return [list { "lj" 1.1 6.21 } { "wca" 1.0 7.45 } \
          { "lj" 7.5 3.86 } { "wca" 1.0 3.01 } \
          { "wca" 1.0 3.86 } { "wca" 1.0 9.92 }]
  } 25 { # asp
      return [list { "lj" 1.2 6.18 } { "wca" 1.0 9.27 } \
          { "lj" 4.0 5.48 } { "wca" 1.0 2.39 } \
          { "no" 0.0 0.0 } { "wca" 1.0 13.7 }]
  } 38 { # cys
      return [list { "lj" 1.0 6.16 } { "wca" 1.0 8.32 } \
          { "lj" 7.0 3.82 } { "no" 0.0 0.0 } \
          { "wca" 1.0 1.09 } { "wca" 1.0 1.37 }]
  } 24 { # gln
      return [list { "lj" 1.2 6.45 } { "wca" 1.0 7.74 } \
          { "lj" 7.5 4.03 } { "wca" 1.0 2.50 } \
          { "wca" 1.0 4.60 } { "wca" 1.0 8.05 }]
  } 23 { # glu
      return [list { "lj" 1.2 6.41 } { "wca" 1.0 9.94 } \
            { "wca" 1.0 5.71 } { "no" 0.0 0.0 } \
            { "wca" 1.0 5.71 } { "wca" 1.0 10.28 }]
  } 20 { # gly
      return [list { "lj" 1.0 5.63 } { "wca" 1.0 8.73 } \
            { "lj" 4.5 4.44 } { "wca" 1.0 4.34 } \
            { "lj" 0.7 3.94 } { "lj" 0.7 3.94 }]
  } 28 { # his
      return [list { "no" 0.0 0.0 } { "no" 0.0 0.0 } \
            { "no" 0.0 0.0 } { "no" 0.0 0.0 } \
            { "no" 0.0 0.0 } { "no" 0.0 0.0 }]
  } 33 { # ile
      return [list { "lj" 1.2 6.61 } { "wca" 1.0 9.92 } \
          { "lj" 5.5 5.91 } { "wca" 1.0 6.73 } \
          { "lj" 1.45 3.55 } { "lj" 1.5 3.55 }]
  } 34 { # leu
      return [list { "lj" 1.2 6.61 } { "wca" 1.0 10.25 } \
          { "lj" 5.5 5.91 } { "wca" 1.0 6.73 } \
          { "lj" 1.2 3.55 } { "lj" 1.2 3.55 }]
  } 29 { # lys
      return [list { "lj" 1.2 6.63 } { "wca" 1.0 6.63 } \
          { "lj" 6.5 4.15 } { "wca" 1.0 2.57 } \
          { "wca" 1.0 4.74 } { "wca" 1.0 10.67 }]
  } 35 { # met
      return [list { "lj" 1.0 6.59 } { "wca" 1.0 8.96 } \
          { "lj" 7.5 4.36 } { "no" 0.0 0.0 } \
          { "wca" 1.0 1.47 } { "wca" 1.0 1.18 }]
  } 36 { # phe
      return [list { "lj" 1.2 6.77 } { "wca" 1.0 9.82 } \
          { "lj" 7.5 5.46 } { "wca" 1.0 6.37 } \
          { "lj" 1.1 3.64 } { "lj" 1.2 3.64 }]
  } 22 { # pro 
      return [list { "lj" 1.0 6.20 } { "wca" 1.0 8.99 } \
          { "lj" 7.5 4.95 } { "wca" 1.0 0.60 } \
          { "wca" 1.0 4.68 } { "wca" 1.0 8.25 }]
  } 27 { # ser
      return [list { "lj" 1.0 5.97 } { "wca" 1.0 8.24 } \
          { "lj" 5.5 4.22 } { "wca" 1.0 1.73 } \
          { "wca" 1.0 3.16 } { "wca" 1.0 10.01 }]
  } 31 { # thr
      return [list { "lj" 1.0 6.23 } { "wca" 1.0 9.28 } \
          { "lj" 7.5 4.98 } { "wca" 1.0 0.60 } \
          { "wca" 1.0 4.15 } { "wca" 1.0 7.74 }]
  } 39 { # trp
      return [list { "lj" 1.2 6.99 } { "wca" 1.0 7.69 } \
          { "lj" 7.0 4.15 } { "no" 0.0 0.0 } \
          { "wca" 1.0 1.89 } { "wca" 1.0 3.15 }]
  } 37 { # tyr
      return [list { "lj" 1.2 6.79 } { "wca" 1.0 8.15 } \
          { "lj" 6.5 4.45 } { "wca" 1.0 1.98 } \
          { "wca" 1.0 3.65 } { "wca" 1.0 4.87 }]
  } 32 { # val
      return [list { "lj" 1.2 6.42 } { "wca" 1.0 10.08 } \
          { "lj" 7.0 5.15 } { "wca" 1.0 1.87 } \
          { "lj" 1.0 3.43 } { "lj" 1.0 3.43 }]
  } 40 { # end
      return [list { "wca" 1.0 4.0 } { "wca" 1.0 4.0 } \
          { "wca" 1.0 4.0 } { "wca" 1.0 4.0 } \
          { "wca" 1.0 4.0 } { "wca" 1.0 4.0 }]
  } default {
      ::mmsg::err [namespace current] "No such residue $aa_name defined."
  }
}


# lipid bead types CH:0; PH:1; GL:2; ES:3 4; AS:5; AD:6; AE:7
set lipid_beads [list {0} {1} {2} {3 4} {5 6} {7}]
set lj_offset 0.0
set lj_cutoff 15.0


# loop over all amino acids
for { set cb_type 20 } { $cb_type < 40 } { incr cb_type } {
  # read in lipid-peptide parameters for side chain cb_type
  set interaction_aa_lipid [aa_lipid_inter $cb_type]
 
  # loop over 6 different lipid bead types (as far as aa-lipid interaction is concerned)
  for { set bead_type 0 } { $bead_type < [llength $lipid_beads] } { incr bead_type } {
    # read interaction
    set inter_details [lindex $interaction_aa_lipid $bead_type]
    set inter_type [lindex $inter_details 0]
    set inter_eps [lindex $inter_details 1]
    set inter_sig [lindex $inter_details 2] 
    
    # loop over all types of beads (e.g., AS -> AS & AD)
    foreach type [lindex $lipid_beads $bead_type] {     
      if { $inter_type == "lj" } {
        # Separate potential in two parts: repulsive (WCA-like) and attractive (LJ)

        # 1) repulsive WCA-like. Vary only if lambda < 0.5
        set wca_eps   $inter_eps ;# Will need to incorporate lambda coupling
        set wca_sig   [expr $inter_sig * $sigmaScale]
        set wca_cut   [expr $wca_sig * sqrt(pow(2,1/3.) \
                             -(1-$lambdaWCA)*$delta)]
        set wca_shift [expr 0.25*(1-$lambdaLJ)]
        set wca_off   0.0
        set wca_cap   0.0
        set wca_soft  ""
        set wca_command "$type $cb_type lj-gen \
          $wca_eps $wca_sig $wca_cut $wca_shift $wca_off \
          12 6 1.0 1.0"
        if { $peptideb::softcore_flag != 0 } {
          set wca_soft " 1.0 $lambdaWCA $delta"
          append wca_command $wca_soft
        }
        lappend nb_interactions $wca_command

        # 2) attractive LJ-like. Vary only between 0.5 < lambda < 1.
        #    Interaction is turned off at lambda <= 0.5.
        #    Lambda coupling only enters in the epsilon parameter.
        set lj_eps   [expr $lambdaLJ * $inter_eps]
        set lj_sig   [expr $inter_sig * $sigmaScale]
        set lj_cut   $lj_cutoff
        set lj_shift [calc_lj_shift $lj_sig $lj_cut]
        set lj_off   0.0
        set lj_cap   0.0
        set lj_min   $wca_cut
        lappend nb_interactions [list $type $cb_type lennard-jones \
          $lj_eps $lj_sig $lj_cut $lj_shift $lj_off $lj_cap $lj_min]

      } elseif { $inter_type == "wca" } {
        set wca_eps   $inter_eps
        set wca_sig   [expr $inter_sig * $sigmaScale]
        set wca_cut   [expr $wca_sig * sqrt(pow(2,1/3.) \
                                -(1-$lambdaWCA)*$delta)]
        set wca_shift 0.25
        set wca_off   0.0
        set wca_cap   0.0
        set wca_soft  ""
        set wca_command "$type $cb_type lj-gen \
          $wca_eps $wca_sig $wca_cut $wca_shift $wca_off \
          12 6 1.0 1.0"
        if { $peptideb::softcore_flag != 0 } {
          set wca_soft " 1.0 $lambdaWCA $delta"
          append wca_command $wca_soft

        }
        lappend nb_interactions $wca_command
      } elseif { $inter_type == "no" } {
        # No interaction
      } else {
        puts "Unknown interaction $inter_type"
        exit
      }
    }
  }
}


