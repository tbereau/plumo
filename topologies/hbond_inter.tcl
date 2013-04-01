# Hydrogen bonding
# Careful with the bonded partners !  
# if {$::cgtools::peptideb::HB_bilayer_dz > 0.} {
#     inter 8 11 lj-angle $::cgtools::peptideb::ljangle_eps \
# 				 $::cgtools::peptideb::hbond_NC $::cgtools::peptideb::ljangle_cut 1 -1 2 \
# 				 -2 0 $::cgtools::peptideb::HB_bilayer_z0 $::cgtools::peptideb::HB_bilayer_dz \
# 				 $::cgtools::peptideb::HB_bilayer_kappa $::cgtools::peptideb::ljangle_eps_bilayer
# } else {
    inter 8 11 lj-angle $::cgtools::peptideb::ljangle_eps \
				 $::cgtools::peptideb::hbond_NC $::cgtools::peptideb::ljangle_cut 1 -1 2 -2
# }
