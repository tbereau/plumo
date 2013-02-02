# ::cgtools::analysis::analyze_rdf_intermol --
#
# Calculate the rdf_intermol of the system.  
# Author: Zun-Jing Wang 
# Dec. 2008

namespace eval ::cgtools::analysis {}

namespace eval ::cgtools::analysis::rdf_intermol {
    variable rdf_intermol_infolist
    variable rlist 
    variable avg_rdf_intermollist_list
    variable rdf_intermol_cnt
    variable rdf_intermol_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf_intermol
#newnew
    variable rdf_intermol_output_infolist 
    variable rdf_intermol_merge_coeflist 

    namespace export setup_rdf_intermol
    namespace export analyze_rdf_intermol
    namespace export printav_rdf_intermol
    namespace export resetav_rdf_intermol
}

proc ::cgtools::analysis::rdf_intermol::resetav_rdf_intermol { } {
    variable avg_rdf_intermollist_list
    variable rdf_intermol_cnt
    variable n_bin 
    variable n_rdf_intermol

    set rdf_intermol_cnt 0
    set avg_rdf_intermollist_list [::cgtools::utils::init_matrix $n_rdf_intermol $n_bin] 
}

proc ::cgtools::analysis::rdf_intermol::printav_rdf_intermol { } {
   global ::cgtools::analysis::outputdir
   variable rdf_intermol_infolist
   variable rlist 
   variable avg_rdf_intermollist_list
   variable rdf_intermol_cnt
   variable rdf_intermol_cntfile
   variable range_min 
   variable range_max 
   variable n_bin 
   variable n_rdf_intermol
#newnew
   variable rdf_intermol_output_infolist 
   variable rdf_intermol_merge_coeflist 

   mmsg::send [namespace current] "printing rdf_intermol at rdf_intermol_cnt = $rdf_intermol_cnt"

   if {  $rdf_intermol_cnt > 0 } {

     # If <$outputdir/rdfcgintermol> doesn't exist then create it
     catch { exec mkdir $outputdir/rdfcgintermol }

#newnew
     set n_rdf_intermoloutput [llength $rdf_intermol_output_infolist]
     mmsg::send [namespace current] "n_rdf_intermoloutput = $n_rdf_intermoloutput"
     for { set i_rdf_intermoloutput 0 } { $i_rdf_intermoloutput  < $n_rdf_intermoloutput } { incr i_rdf_intermoloutput } {

        set current_rdf_intermollist [lindex $rdf_intermol_output_infolist $i_rdf_intermoloutput]
        set nrdf_intermol_plus [llength $current_rdf_intermollist]

        set i_rdf_intermol [lindex $current_rdf_intermollist 0]

        set rdf_intermol_info_now [lindex $rdf_intermol_infolist $i_rdf_intermol]
        set typei [lindex $rdf_intermol_info_now 0]
        set typej [lindex $rdf_intermol_info_now 1]
        set rdf_intermolfile [lindex $rdf_intermol_info_now 2]

        set nbin_rlist [llength $rlist]

        set avg_rdf_intermollist [lindex $avg_rdf_intermollist_list $i_rdf_intermol]
        set rscale_rdf_intermol [expr 1.0/($rdf_intermol_cnt*1.0)]
        set avg_rdf_intermollist [::cgtools::utils::scalevec $avg_rdf_intermollist $rscale_rdf_intermol]

        if { $nrdf_intermol_plus > 1 } {
           #mmsg::send [namespace current] "merge rdf_intermol $rdf_intermolfile, # of files is  $nrdf_intermol_plus"
           if { $nrdf_intermol_plus == 2 } {
                set coeffnow [lindex $rdf_intermol_merge_coeflist 0]
                if { [llength $coeffnow] != $nrdf_intermol_plus } {
                       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
                }
           }
           if { $nrdf_intermol_plus == 3 } {
                set coeffnow [lindex $rdf_intermol_merge_coeflist 1]
                if { [llength $coeffnow] != $nrdf_intermol_plus } {
                       mmsg::err [namespace current] "coeffnow is wrong: $coeffnow "
                } 
           }
           if { $nrdf_intermol_plus != 3 && $nrdf_intermol_plus != 2 } {
                mmsg::err [namespace current] "nrdf_intermol_plus is wrong: nrdf_intermol_plus = $nrdf_intermol_plus "
           }

           set f_rdf_intermol [open "$rdf_intermolfile.0" w]
           puts $f_rdf_intermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
           foreach r $rlist rdf_intermol $avg_rdf_intermollist {
                puts $f_rdf_intermol "$r $rdf_intermol"
           }
           close $f_rdf_intermol
            
           set avg_rdf_intermollist [::cgtools::utils::scalevec $avg_rdf_intermollist [lindex $coeffnow 0] ]


           for { set irdf_intermol_index 1 } { $irdf_intermol_index  < $nrdf_intermol_plus } { incr irdf_intermol_index} {

                set i_rdf_intermol [lindex $current_rdf_intermollist $irdf_intermol_index]
                set avg_rdf_intermolpluslist [lindex $avg_rdf_intermollist_list $i_rdf_intermol]

                set rscale_rdf_intermol [expr 1.0/($rdf_intermol_cnt*1.0)]
                set avg_rdf_intermolpluslist [::cgtools::utils::scalevec $avg_rdf_intermolpluslist $rscale_rdf_intermol]

           	set f_rdf_intermol [open "$rdf_intermolfile.$irdf_intermol_index" w]
           	puts $f_rdf_intermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
           	foreach r $rlist rdf_intermol $avg_rdf_intermolpluslist {
                	puts $f_rdf_intermol "$r $rdf_intermol"
           	}
           	close $f_rdf_intermol

                set coeffvalue [lindex $coeffnow $irdf_intermol_index]
                set avg_rdf_intermolpluslist [::cgtools::utils::scalevec $avg_rdf_intermolpluslist $coeffvalue ]

                set avg_rdf_intermollist [::cgtools::utils::add2vec $avg_rdf_intermollist $avg_rdf_intermolpluslist]
           }
        }

        set f_rdf_intermol [open "$rdf_intermolfile" w]
        puts $f_rdf_intermol "\# $nbin_rlist     [lindex $rlist 0]       [lindex $rlist [expr $nbin_rlist - 1] ]"
        foreach r $rlist rdf_intermol $avg_rdf_intermollist {
                puts $f_rdf_intermol "$r $rdf_intermol"
        }
        close $f_rdf_intermol
        #after 100
     } 
     #end for { set i_rdf_intermoloutput 0 }

     # backup <$outputdir/rdfcgintermol> 
     incr rdf_intermol_cntfile

     # checkwith if the filedir exists already 
     set newdir 0
     while { $newdir == 0 } { 
    	if { [file isdirectory  $outputdir/rdfcgintermol.$rdf_intermol_cntfile ] } {
		incr rdf_intermol_cntfile
    	} else {
		set newdir 1
	}
     }

     puts "rdf_intermol_cntfile = $rdf_intermol_cntfile"
     catch { exec cp -r $outputdir/rdfcgintermol  $outputdir/rdfcgintermol.$rdf_intermol_cntfile }

   }
}

#newnew
proc ::cgtools::analysis::rdf_intermol::setup_rdf_intermol { rdf_intermolcglist rdf_intermolcgoutputlist args } {
    variable rdf_intermol_infolist
    variable avg_rdf_intermollist_list
    variable rdf_intermol_cnt
    variable rdf_intermol_cntfile
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf_intermol
#newnew
    variable rdf_intermol_output_infolist 
    variable rdf_intermol_merge_coeflist 
    variable ::cgtools::analysis::topology

    set options {
	{rmin.arg "0."  "minimum distange in rdf_intermol computation"}
	{rmax.arg "15.01" "maximum distange in rdf_intermol computation"}
	{nbin.arg "1501" "number of bins in rdf_intermol computation "}
    }
    set usage "Usage: ::cgtools::analysis::rdf_intermol::setup_rdf_intermol rdf_intermolcglist rdf_intermolcgoutputlist \[rmin:rmax:nbin]"
    array set params [::cmdline::getoptions args $options $usage]

    set rdf_intermol_infolist $rdf_intermolcglist
    set n_rdf_intermol [llength $rdf_intermol_infolist]
#newnew
    set rdf_intermol_output_infolist $rdf_intermolcgoutputlist

    #set rdf_intermol_merge_coeflist which is used in printav_rdf
    set nmols [llength $topology]
    mmsg::send [namespace current] "nmols = $nmols"

    set coeflist_now 0.5
    lappend coeflist_now 0.5
    lappend rdf_intermol_merge_coeflist $coeflist_now

    set coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr ($nmols * 1.0)/(2.*$nmols-1.) ]
    lappend coeflist_now [expr 0.5*($nmols*1.0-1.0)/(2.*$nmols-1.) ]
    lappend rdf_intermol_merge_coeflist $coeflist_now
    mmsg::send [namespace current] "$rdf_intermol_merge_coeflist"


    set rdf_intermol_cnt 0
    set rdf_intermol_cntfile 0
    set range_min $params(rmin)
    set range_max $params(rmax)
    set n_bin $params(nbin)

    set avg_rdf_intermollist_list [::cgtools::utils::init_matrix $n_rdf_intermol $n_bin] 
}

proc ::cgtools::analysis::rdf_intermol::analyze_rdf_intermol { } {
    variable rdf_intermol_infolist
    variable rlist 
    variable avg_rdf_intermollist_list
    variable rdf_intermol_cnt
    variable range_min 
    variable range_max 
    variable n_bin 
    variable n_rdf_intermol

    mmsg::send [namespace current] "analyzing rdf_intermol"

    set rdf_intermollist_list_now "" 
    for { set i_rdf_intermol 0 } { $i_rdf_intermol < $n_rdf_intermol } { incr i_rdf_intermol } {
	set rdf_intermol_info_now [lindex $rdf_intermol_infolist $i_rdf_intermol]
	set typei [lindex $rdf_intermol_info_now 0]
	set typej [lindex $rdf_intermol_info_now 1]
	set rdf_intermolfile [lindex $rdf_intermol_info_now 2]

    	set rdf_intermol [analyze rdf-intermol $typei $typej $range_min $range_max $n_bin]
    	set rlist ""
    	set rdf_intermollist ""
    	foreach value [lindex $rdf_intermol 1] {
   		lappend rlist [lindex $value 0]
# time factor 0.5 to the current rdf_intermol, since the boxsize is 2 times big of the AA simulations
    		lappend rdf_intermollist [expr [lindex $value 1]*0.5] 
	}
	lappend rdf_intermollist_list_now $rdf_intermollist
    }
    #mmsg::send [namespace current] "$rdf_intermollist_list_now"
    set avg_rdf_intermollist_list [::cgtools::utils::add_matrixs $avg_rdf_intermollist_list $rdf_intermollist_list_now] 

    incr rdf_intermol_cnt

    mmsg::send [namespace current] "analyzing rdf_intermol done"
}
