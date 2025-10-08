
    set PROBEDATA [dict create]
    set PROBETYPES [dict create core "Core probes"]
    dict set PROBETYPES polar "Polar probes"
    dict set PROBETYPES hydrophobe "Hydrophobes"
    dict set PROBETYPES negative "Negatively charged"
    dict set PROBETYPES positive "Positively charged"
    dict set PROBETYPES ring5 "5-membered rings"
    dict set PROBETYPES ring6 "6-membered rings"
    set PROBETOPPAR [dict create PBDA "probe2.top probe.prm"]
    dict set PROBETOPPAR CGenFF "top_all36_cgenff.rtf par_all36_cgenff.prm"
    set PACKAGEPATH /Users/carlosventura/Desktop/prody_drugui/ProDy/prody/drugui/DruGUI-script
    
    foreach {key} [dict keys $PROBETYPES] {
    dict set PROBEDATA $key ""
    }

    set inf [open [file join $PACKAGEPATH "probeV2.dat"] r]
    foreach line [split [read -nonewline $inf] \n] {
        if {[llength $line] < 3 || [string equal -length 1 $line "#"]} {continue}
        set resi [lindex $line 0]
        set key [lindex $line 1]
        if {![dict exists $PROBEDATA $resi]} {
            dict set PROBEDATA $resi default 0
            dict set PROBEDATA $resi alias ""
            dict set PROBEDATA $resi atomnames ""
            dict set PROBEDATA $resi charge 0
            dict set PROBEDATA $resi source ""
        }
        if {$key == "default" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA defaults $resi
        }
        if {$key == "type" && [lindex $line 2] > 0} {
            dict lappend PROBEDATA [lindex $line 2] $resi
        }
        dict set PROBEDATA $resi $key [lrange $line 2 end]
        if {$key == "alias" && [expr [llength [dict get $PROBEDATA $resi $key]] % 2] != 0} {
            error "Problem in aliases of $resi!"
        }
     }
    close $inf

    set ipronames "C2 H21 C1 H11 H12 H13 C3 C4 C5 C6 H31 H32 H33 OH2 HO2"
    dict for {src toppar} $PROBETOPPAR {
        set inf [open [file join $PACKAGEPATH [lindex $toppar 0]] r]
        puts "$inf"
        set resi "____"
        set prev "____"
        foreach line [split [read -nonewline $inf] \n] {
            if {[dict exists $PROBEDATA $resi] && [string range $line 0 3] == "ATOM"} {
                dict set PROBEDATA $resi atomnames "[dict get $PROBEDATA $resi atomnames] [lindex $line 1]"
            } elseif {[string range $line 0 3] == "RESI"} {
              set prev $resi
              set resi [string trim [lindex $line 1]]
              if {[dict exists $PROBEDATA $resi]} {
                dict set PROBEDATA $resi charge [lindex $line 2]
                dict set PROBEDATA $resi source $src
              }
              if {[dict exists $PROBEDATA $prev]} {
                set thisnames [dict get $PROBEDATA $prev atomnames]
                set remove 0
                set aliases [dict get $PROBEDATA $prev alias]
                foreach {name alias} $aliases {
                    if {[lsearch $ipronames $name] == -1} {
                        error "Invalid $prev alias, $name is not found in IPRO!"
                        #set remove 1
                    } elseif {[lsearch $thisnames $alias] == -1} {
                      error "Invalid $prev alias, $alias is not found in $prev!"
                    #set remove 1
                    }
                    #if {$remove} {
                    #  set PROBEDATA [dict remove $PROBEDATA $prev]
                    #  break
                    #}
                }
            }
        }
     }
    close $inf
    }

    set probe_data [open "probedata.txt" w]
    puts $probe_data $PROBEDATA
    close $probe_data

    set probetypes [open "probetypes.txt" w]
    puts $probetypes $PROBETYPES
    close $probetypes

    set probetoppar [open "probetoppar.txt" w]
    puts $probetoppar $PROBETOPPAR
    close $probetoppar
    exit
    