# NMWiz: Normal Mode Visualization, Animation, and Plotting
# 
# University of Illinois Open Source License
# Copyright 2010-2011 Ahmet Bakan
# All rights reserved.
# 
# Developed by:		
#       Ahmet Bakan
# 			http://www.pitt.edu/~ahb12/
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the Software), to deal with 
# the Software without restriction, including without limitation the rights to 
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
# of the Software, and to permit persons to whom the Software is furnished to 
# do so, subject to the following conditions:
# 
# Redistributions of source code must retain the above copyright notice, 
# this list of conditions and the following disclaimers.
# 
# Redistributions in binary form must reproduce the above copyright notice, 
# this list of conditions and the following disclaimers in the documentation 
# and/or other materials provided with the distribution.
# 
# Neither the names of Theoretical and Computational Biophysics Group, 
# University of Illinois at Urbana-Champaign, nor the names of its contributors 
# may be used to endorse or promote products derived from this Software without 
# specific prior written permission.
# 
# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
# THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS WITH THE SOFTWARE.

package require exectool 1.2
package provide nmwiz 0.8

# 2-D plotting tool
#
# $Id: multiplot.tcl,v 1.31 2009/07/20 21:32:45 saam Exp $
#
# Author:
# Jan Saam
# Beckman Institute
# University of Illinois
# saam@ks.uiuc.edu
namespace eval ::nmwiz_MultiPlot:: {
   proc initialize {} {
      variable plotlist {}
      variable plotcount -1
      variable parent
      variable verbose 0
   }
   initialize
}
proc sign x {expr {($x>0) - ($x<0)}}
proc ::nmwiz_MultiPlot::init_plot {args} {
   variable parent
   variable verbose

   set  parent [lindex $args 0]
   incr ::nmwiz_MultiPlot::plotcount
   set ns "::nmwiz_MultiPlot::Plot${::nmwiz_MultiPlot::plotcount}"

   if {$verbose} {
     if {[namespace exists $ns]} {
       vmdcon -info "Reinitializing namespace $ns."
     } else {
       vmdcon -info "Creating namespace $ns"
     }
   }

   namespace eval $ns {
      # Default values
      variable nsets 0
      variable datasets
      array unset datasets
      lappend args -set 0
      variable curset 0
      variable title {}
      variable titlefontsize 10
      variable ticfontsize   8
      variable labelfontsize 10
      variable titlefont  "Helvetica $titlefontsize"
      variable ticfont    "Helvetica $ticfontsize"
      variable labelfont  "Helvetica $labelfontsize"
      variable infoFont   {Courier 9}
      variable postscript "multiplot.ps"
      variable printstats 0
      variable replot 0
      variable canh   700;   # canvas height
      variable canw   1000;  # canvas width
      variable resize 0
      variable ticlen 10; # length of the major tic marks
      variable rim    8;  # extra space around the plot
      variable xlabeloffset
      variable xlabeltext ""
      variable ylabeloffset
      variable ylabeltext ""
      variable lines     1;        # connect data points with lines?
      variable marker    "none";   # display data points [circle|point|square|none]
      variable radius    2;        # radius of circles and points , size of squares
      variable linewidth 1;        # width of the line connecting data points
      variable fillcolor Skyblue2; # fill color of data point markers   
      variable linecolor black;    # color of lines connecting data points
      variable dashed    {{}};     # Draw dashed lines (uses the same format as -dash for Tk canvas)
      variable legend    {{}};     # legend string for current dataset
      variable colorlist {black red green blue magenta orange OliveDrab2 cyan maroon gold2 yellow gray60 SkyBlue2 orchid3 ForestGreen PeachPuff LightSlateBlue}

      variable predefRange 0
      variable givenXmin auto
      variable givenYmin auto
      variable givenXmax auto
      variable givenYmax auto

      variable xmin   0
      variable xmin_y 0
      variable ymin   0
      variable ymin_x 0
      variable xmax   0
      variable xmax_y 0
      variable ymax   0
      variable ymax_x 0
      variable spanx  0
      variable spany  0
      variable anglescalex 0
      variable anglescaley 0
      variable xmajortics {}
      variable ymajortics {}
      variable xminortics {}
      variable yminortics {}

      variable hline    {}
      variable vline    {}
      variable xplotmin {}
      variable yplotmin {}
      variable xplotmax {}
      variable yplotmax {}
      variable scalex {}
      variable scalex {}
      variable dimx 0
      variable dimy 0
      variable minorticx 5
      variable minorticy 5

      variable objectlist {}; # other drawn objects like circles, text, lines,...
      variable redraw 0;      # redraw all objects

      variable w ${::nmwiz_MultiPlot::parent}.plotwindow${::nmwiz_MultiPlot::plotcount}
      variable istoplevel 1;  # set to 0 for embedded widgets
      variable c
      variable namespace ::Plothandle${::nmwiz_MultiPlot::plotcount}

      if {${::nmwiz_MultiPlot::parent} != ""} {
        set istoplevel 0
        set canh 350
        set canw 500
      } else {
        set istoplevel 1
      } 

      catch {destroy $w}
      if {$istoplevel} {
        toplevel $w
        wm title $w "MultiPlot"
        wm iconname $w "MultiPlot"
        wm protocol $w WM_DELETE_WINDOW "[namespace current]::plothandle quit"
        wm withdraw $w
      } else {
        frame $w -bd 0
        pack $w -side top -fill x -fill y
      }

      frame $w.menubar -relief raised -bd 2
      menubutton $w.menubar.file -text "File" -underline 0 \
         -menu $w.menubar.file.menu
      $w.menubar.file config -width 3 

      menu $w.menubar.file.menu -tearoff 0
      
      $w.menubar.file.menu add command -label "Export to PostScript" -command "[namespace current]::savedialog "
      $w.menubar.file.menu add command -label "Export to Xmgrace" -command "[namespace current]::xmgracedialog "
      $w.menubar.file.menu add command -label "Print plothandle in Console" -command "vmdcon [namespace current]::plothandle"
      if {$istoplevel} {
        $w.menubar.file.menu add command -label "Quit" -command "[namespace current]::plothandle quit"
      }
      pack $w.menubar.file -side left
      pack $w.menubar -anchor w -fill x

      if {![winfo exists $w.f.cf]} {
         variable canw
         variable canh
         variable c $w.f.cf
         frame $w.f 
         pack $w.f -fill x -fill y 

         canvas $c -relief flat -borderwidth 0 -width $canw -height $canh 
         scrollbar $w.f.y -orient vertical   -command [namespace code {$c yview}]
         scrollbar $w.f.x -orient horizontal -command [namespace code {$c xview}]
         $c configure  -yscrollcommand [namespace code {$w.f.y set}] -xscrollcommand [namespace code {$w.f.x set}]
         $c configure  -scrollregion   "0 0 $canw $canh"
         grid $c $w.f.y 
         grid $w.f.x    
         grid rowconfigure    $w.f 0 -weight 1
         grid columnconfigure $w.f 0 -weight 1
         grid configure $w.f.y  -sticky ns
         grid configure $w.f.x  -sticky we
      }

      # Create a plothandle procedure that provides some commands to control the plot.
      # It's full name will be returned when you invoke multiplot.
      proc plothandle { command args } {
         variable w
         switch $command {
            namespace { return [namespace current] }
            replot    { variable replot 1; plot_update; return }
            add       { 
               set newX [lindex $args 0]
               set newY [lindex $args 1]

               set lenX [llength $newX]
               set lenY [llength $newY]
               if {!$lenX} { error "X data vector is empty!" }
               if {!$lenY} { error "Y data vector is empty!" }

               # Check, if we have several coordinate sets:
               if {[llength [join $newX]]>$lenX || [llength [join $newY]]>$lenY} {
                  if {$lenX != $lenY} {
                     error "Different number of datasets for x and y ($lenX!=$lenY)"
                  }
                  foreach x $newX y $newY {
                     eval add_data [list $x] [list $y] [lrange $args 2 end]
                  }
               } else {
                  eval add_data [list $newX] [list $newY] [lrange $args 2 end]
               }
               plot_update 
            }
            draw {
               # Register the new object
               variable objectlist
               lappend objectlist $args

               # Make sure that the plot geometry was calculated already and draw
               variable xplotmin
               if {![llength $xplotmin]} {
                  variable redraw 1;
                  variable resize 1;
                  plot_update; # implicitely draws all objects
               } else {
                  draw_object $args
               }
            }
             undraw {
                 undraw_object $args
             }
            configure { 
               variable datasets
               variable nsets
               variable resize

               variable curset {} 
               set pos [lsearch $args "-set"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  variable curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
               }
               if {![llength $curset]} { set curset 0 }

               set havedata 0
               set pos [lsearch $args "-x"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  if {$nsets==0} {
                     lappend datasets(X) {}
                     lappend datasets(Y) {}
                     lappend datasets(xmin) {}
                     lappend datasets(xmax) {}
                     lappend datasets(xmin_y) {}
                     lappend datasets(xmax_y) {}
                     lappend datasets(ymin) {}
                     lappend datasets(ymax) {}
                     lappend datasets(ymin_x) {}
                     lappend datasets(ymax_x) {}
                     incr nsets
                  }
                  lset datasets(X) $curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
                  variable resize 1
                  set havedata 1
               }

               set pos [lsearch $args "-y"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  if {$nsets==0} {
                     lappend datasets(X) {}
                     lappend datasets(Y) {}
                     lappend datasets(xmin) {}
                     lappend datasets(xmax) {}
                     lappend datasets(xmin_y) {}
                     lappend datasets(xmax_y) {}
                     lappend datasets(ymin) {}
                     lappend datasets(ymax) {}
                     lappend datasets(ymin_x) {}
                     lappend datasets(ymax_x) {}
                     incr nsets
                  }
                  lset datasets(Y) $curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
                  variable resize 1
                  set havedata 1
               }

               plot_scan_options $args; 

               if {$resize && $havedata} {
                  init_dataset
               }

               plot_update 
            }
            nsets      {
               variable datasets;
               return [llength $datasets(Y)]
            }
            xdata      {
               variable datasets;
               return $datasets(X)
            }
            ydata      {
               variable datasets;
               return $datasets(Y)
            }
            data      { 
               variable datasets;
               return [list $datasets(X) $datasets(Y)]
            }
            getpath   {
               variable w;
               return $w
            }
            export    {
              variable datasets
              variable title
              variable legend
              variable nsets
              if { [llength $args] < 2} {
                vmdcon -err "Incorrect export syntax"
                return
              }
              set progname [lindex $args 0]
              set filename [lindex $args 1] 

              switch $progname {
                grace - 
                xmgr  -
                xmgrace {
                  vmdcon -info "Exporting plot in xmgrace format as filename $filename"
                  set fd [open $filename "w"]
                  puts $fd "@type xy"
                  puts $fd "@title "$title""
                  set ylen [llength $datasets(Y)]
                  for {set s 0} {$s < $nsets} {incr s} {
                    if {[lindex $legend $s] != ""} {
                      puts $fd "@s$s legend "[lindex $legend $s]""
                    }
                  }
                  for {set s 0} {$s < $nsets} {incr s} {
                    set xlen [llength [lindex $datasets(X) $s]]
                    for {set i 0} {$i < $xlen} {incr i} {
                      puts $fd "[lindex $datasets(X) $s $i] [lindex $datasets(Y) $s $i]"
                    }
                    puts $fd "&"
                  }                    
                  close $fd
                }
              }

              return
            }
            quit   { 
               destroy $w;
               namespace delete [namespace current]
               return
            }
         }
      }

      proc init_dataset {} {
         variable datasets
         variable curset
         set minx [lindex $datasets(X) $curset 0]
         set maxx [lindex $datasets(X) $curset end]
         set miny [lindex $datasets(Y) $curset 0]
         set maxy [lindex $datasets(Y) $curset end]
         set minx_y [lindex $datasets(Y) $curset 0]
         set maxx_y [lindex $datasets(Y) $curset end]
         set miny_x [lindex $datasets(X) $curset 0]
         set maxy_x [lindex $datasets(X) $curset end]
         foreach x [lindex $datasets(X) $curset] y [lindex $datasets(Y) $curset] {
            if {$x<$minx} {
               set minx   $x
               set minx_y $y
            }
            if {$x>$maxx} {
               set maxx   $x
               set maxx_y $y
            }
            if {$y<$miny} {
               set miny   $y
               set miny_x $x
            }
            if {$y>$maxy} {
               set maxy   $y
               set maxy_x $x
            }
         }
         lset datasets(xmin)   $curset $minx
         lset datasets(xmin_y) $curset $minx_y
         lset datasets(xmax)   $curset $maxx
         lset datasets(xmax_y) $curset $maxx_y
         lset datasets(ymin)   $curset $miny
         lset datasets(ymin_x) $curset $miny_x
         lset datasets(ymax)   $curset $maxy
         lset datasets(ymax_x) $curset $maxy_x
      }

      proc plot_scan_options { arg } {
         set drawlines {}
         set points 0
         variable printstats
         variable anglescalex
         variable anglescaley

         # Scan for single options
         set argnum 0
         set arglist $arg
         foreach i $arg {
            if {$i=="-lines"}  then {
               set drawlines 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-nolines"}  then {
               set drawlines 0
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-stats"}  then {
               set printstats 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-plot"}  then {
               variable replot 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-nostats"}  then {
               set printstats 0
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-xanglescale"}  then {
               set anglescalex 1
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-yanglescale"}  then {
               set anglescaley 1
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-autoscale"} then {
               variable predefRange 0
               variable givenXmin auto
               variable givenXmax auto
               variable givenYmin auto
               variable givenYmax auto
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1
               continue
            }
            incr argnum
         }

         # must search for the dataset option first
         variable nsets
         variable curset 
         foreach {i j} $arglist {
            if {$i=="-set"}       then { 
               if {$j>=$nsets} {
                  error "Dataset $j doesn't exist"
               }
               variable curset $j;
            }
         }

         #variable curset
         if {[llength $drawlines]} {
            variable lines
            if {![llength $curset]} {
               for {set s 0} {$s<$nsets} {incr s} {
                  lset lines $s $drawlines
               }
            } else {
               lset lines $curset $drawlines
            }
         }

         # Scan for options with one argument
         variable hline
         variable vline
         variable datasets
         foreach {i j} $arglist {
#           if {$i=="-x"}          then { 
#              if {![llength [array get datasets X]]} {
#                 lappend datasets(X) $curset $j;
#              } else {
#              lset datasets(X) $curset $j
#              }
#              variable resize 1
#           }
#           if {$i=="-y"}          then { 
#              if {![llength [array get datasets Y]]} {
#                 lappend datasets(Y) $curset $j;
#              } else {
#                 lset datasets(Y) $curset $j
#              }
#              variable resize 1
#           }
            if {$i=="-title"}      then { variable title $j; variable resize 1 }
            if {$i=="-xlabel"}     then { variable xlabeltext $j; variable resize 1 }
            if {$i=="-ylabel"}     then { variable ylabeltext $j; variable resize 1 }
            if {$i=="-xmajortics"} then { variable xmajortics $j; variable resize 1 }
            if {$i=="-ymajortics"} then { variable ymajortics $j; variable resize 1 }
            if {$i=="-xminortics"} then { variable xminortics $j; variable resize 1 }
            if {$i=="-yminortics"} then { variable yminortics $j; variable resize 1 }
            if {$i=="-xsize"}      then { variable canw $j; variable resize 1 }
            if {$i=="-ysize"}      then { variable canh $j; variable resize 1 }

            if {$i=="-xmin"}       then { variable givenXmin $j; variable predefRange 1; variable resize 1 }
            if {$i=="-xmax"}       then { variable givenXmax $j; variable predefRange 1; variable resize 1 }
            if {$i=="-ymin"}       then { variable givenYmin $j; variable predefRange 1; variable resize 1 }
            if {$i=="-ymax"}       then { variable givenYmax $j; variable predefRange 1; variable resize 1 }

            if {$i=="-hline"}      then { lappend hline $j; variable resize 1 }
            if {$i=="-vline"}      then { lappend vline $j; variable resize 1 }
            if {$i=="-radius"}     then { 
               variable radius 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset radius $s $j
                  }
               } else {
                  lset radius $curset $j
               }
            }
            if {$i=="-dash"}     then { 
               variable dashed
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset dashed $s $j
                  }
               } else {
                  lset dashed $curset $j
               }
            }
            if {[string match "-fill*" $i]} then { 
               variable fillcolor 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset fillcolor $s $j;
                  }
               } else {
                  lset fillcolor $curset $j
               }
            }
            if {$i=="-linewidth"} then { 
               variable linewidth 
               variable datasets
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset linewidth $s $j
                  }
               } else {
                  lset linewidth $curset $j
               }
            }
            if {$i=="-linecolor"} then {
               variable linecolor 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset linecolor $s $j;
                  }
               } else {
                  lset linecolor $curset $j
               }
            }
            if {[string match "-mark*" $i]} then { 
               variable marker
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset marker $s $j
                  }
               } else {
                  lset marker $curset $j
               }
            }
            if {$i=="-legend"}      then { 
               variable legend
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset legend $s $j
                  }
               } else {
                  lset legend $curset $j
               }
            }
            if {$i=="-nmwizns"}     then {
               variable nmwizns $j
            }
         }
      }

      proc undraw_object {args} {
          variable c
          
          $c delete $args
      }
      
      proc add_data {x y args} {
         if {[llength $x] != [llength $y]} {
            error "Different number of x and y coordinates ([llength $x]!=[llength $y])"
         }
         variable datasets
         variable nsets
         variable curset $nsets
         variable lines 
         variable linewidth 
         variable linecolor
         variable marker
         variable fillcolor 
         variable dashed
         variable radius
         variable legend
         variable colorlist

         lappend datasets(X) $x
         lappend datasets(Y) $y
         lappend datasets(xmin)   {}
         lappend datasets(xmax)   {}
         lappend datasets(xmin_y) {}
         lappend datasets(xmax_y) {}
         lappend datasets(ymin)   {}
         lappend datasets(ymax)   {}
         lappend datasets(ymin_x) {}
         lappend datasets(ymax_x) {}
         lappend lines     1
         lappend linewidth 1
         lappend linecolor [lindex $colorlist [expr {(1+$nsets)%[llength $colorlist]}]]
         lappend marker    "none"
         lappend fillcolor [lindex $colorlist [expr {(1+$nsets)%[llength $colorlist]}]]
         lappend radius    2
         lappend dashed    {}
         lappend legend    {}

         # Evaluate the command line options
         lappend args -set $nsets
         incr nsets
         plot_scan_options $args

         #variable replot 1
         init_dataset
      }

      proc plot_update {} {
         variable datasets
         set lenx [llength [lindex [array get datasets X] 1 0]]
         set leny [llength [lindex [array get datasets Y] 1 0]]
         if {!$leny} {
            vmdcon -warn "multiplot: Data vector empty, ignoring plot!"
            variable replot 0; return
         }
         if {$lenx && $lenx!=$leny} {
            vmdcon -warn "multiplot: Different size of X and Y data, ignoring plot!"
            variable replot 0; return
         }

         variable replot 
         variable redraw

         if {!$replot && !$redraw} { return }

         # Use index number if no X-coordinate was specified
         set j 0
         foreach X $datasets(X) Y $datasets(Y) {
            if {![llength $X]} {
               set x {}
               for {set i 0} {$i<[llength $Y]} {incr i} {
                  lappend x $i
               }
               lset datasets(X) $j $x
               init_dataset
            }
            incr j
         }

         variable w
         variable c
         variable resize
         variable istoplevel

         # Display some statistics in an info frame
         variable printstats
         if {![winfo exists $w.info] && $printstats} { draw_infobox }
         if {[winfo exists $w.info] && !$printstats} { destroy $w.info; pack $w.cf }

         if {[winfo exists $c] && $resize} {
            variable canw
            variable canh
            $c configure -width $canw -height $canh
         }

         calculate_range
         calculate_ticspacing

         if {$resize} {
            # Clear the canvas
            $c addtag all all
            $c delete all
            calculate_labelsize
            calculate_plot_geometry
            draw_periphery
            redraw_objects 
            variable redraw 0
            variable resize 0
         }

         if {$replot} {
            if {$istoplevel} {
              wm deiconify $w
            }
            plot_data
            variable replot 0
         }
      }

      proc plot_data {} {
         # Plot the values
         variable c
         variable datasets
         variable marker
         variable radius
         variable dashed
         variable fillcolor
         variable lines
         variable linewidth
         variable linecolor
         variable xmin
         variable ymin
         variable xmax
	       variable ymax
	       variable xplotmin
         variable xplotmax
         variable yplotmin
         variable yplotmax
         variable scalex
         variable scaley
         variable legend
         variable legendheight
         $c delete legend
         $c delete point
         $c delete lines
         variable nsets
         set ds 0
         foreach X $datasets(X) Y $datasets(Y) dsxmin $datasets(xmin) dsxmin_y $datasets(xmin_y) {
            set stride 1
            #set len [llength $dset(X)]
            #if {$len>[expr 2*($xplotmax-$xplotmin)]} { 
            #   set stride [expr int([llength $dset(X)]/[expr $xplotmax-$xplotmin])]
            #   puts "Using stride $stride for data set $ds ([llength $dset(X)]/[expr $xplotmax-$xplotmin])"
            #}
            set fc   [lindex $fillcolor $ds]
            set lc   [lindex $linecolor $ds]
            set rad  [lindex $radius $ds]
            set dash [lindex $dashed $ds]
            set leg  [lindex $legend $ds]

            if {[lindex $lines $ds]} {
               set i 0
               foreach cx $X cy $Y {
                  set cxf [format "%10g" $cx]
                  set cyf [format "%10g" $cy]

                  incr i
                  if {[expr $i%$stride]} { continue }
                  set x [expr {$xplotmin + ($scalex*($cx-$xmin))}]
                  set y [expr {$yplotmin + ($scaley*($cy-$ymin))}]

                  set outofBounds 1
                  if { $cxf < $xmin } { set outofBounds [expr $outofBounds * 2] }
                  if { $cxf > $xmax } { set outofBounds [expr $outofBounds * 3] }
                  if { $cyf < $ymin } { set outofBounds [expr $outofBounds * 5] }
                  if { $cyf > $ymax } { set outofBounds [expr $outofBounds * 7] }

                  if { $i == 1 } { 
                     set oldcx $cx
                     set oldcy $cy
                     set oldx $x
                     set oldy $y
                     set oldoutofBounds $outofBounds
                     continue
                  }

                  if { $outofBounds == 1 } {
                     if { $oldoutofBounds == 1 } {
                        set item [$c create line $oldx $oldy $x $y -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     } else {
                        set xyinters [calcIntersect $oldcx $oldcy $cx $cy $oldoutofBounds]
                        set xinterplot [expr {$xplotmin + ($scalex*([lindex $xyinters 0] -$xmin))}]
                        set yinterplot [expr {$yplotmin + ($scaley*([lindex $xyinters 1] -$ymin))}]
                        set item [$c create line $xinterplot $yinterplot $x $y -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     }
                  } else {
                     if { $oldoutofBounds == 1 } {
                        set xyinters [calcIntersect $oldcx $oldcy $cx $cy $outofBounds]
                        set xinterplot [expr {$xplotmin + ($scalex*([lindex $xyinters 0] -$xmin))}]
                        set yinterplot [expr {$yplotmin + ($scaley*([lindex $xyinters 1] -$ymin))}]
                        set item [$c create line $oldx $oldy $xinterplot $yinterplot -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     } else {
                        if { ($outofBounds % 2 != 0 || $oldoutofBounds % 2 != 0) && ($outofBounds % 3 != 0 || $oldoutofBounds % 3 != 0) && \
                              ($outofBounds % 5 != 0 || $oldoutofBounds % 5 != 0) && ($outofBounds % 7 != 0 || $oldoutofBounds % 7 != 0) } {
                           set xyinters1 [calcIntersect $oldcx $oldcy $cx $cy $oldoutofBounds]
                           if {$xmin <= [lindex $xyinters1 0] && [lindex $xyinters1 0] <= $xmax && $ymin <= [lindex $xyinters1 1] && [lindex $xyinters1 1] <= $ymax} {
                              set xinterplot1 [expr {$xplotmin + ($scalex*([lindex $xyinters1 0] -$xmin))}]
                              set yinterplot1 [expr {$yplotmin + ($scaley*([lindex $xyinters1 1] -$ymin))}]
                              set xyinters2 [calcIntersect $oldcx $oldcy $cx $cy $outofBounds]
                              if {$xmin <= [lindex $xyinters2 0] && [lindex $xyinters2 0] <= $xmax && $ymin <= [lindex $xyinters2 1] && [lindex $xyinters2 1] <= $ymax} {               
                                 set xinterplot2 [expr {$xplotmin + ($scalex*([lindex $xyinters2 0] -$xmin))}]
                                 set yinterplot2 [expr {$yplotmin + ($scaley*([lindex $xyinters2 1] -$ymin))}]
                                 set item [$c create line $xinterplot1 $yinterplot1 $xinterplot2 $yinterplot2 -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                                 $c addtag lines withtag $item
                              }
                           }
                        }
                     }
                  }
                  set oldcx $cx
                  set oldcy $cy
                  set oldx $x
                  set oldy $y
                  set oldoutofBounds $outofBounds
               }
            }

            if {[lindex $marker $ds]!="none"} {
               set i 0
               foreach cx $X cy $Y {
                  set cxf [format "%10g" $cx]
                  set cyf [format "%10g" $cy]
                  if { $cxf >= $xmin && $cxf <= $xmax && $cyf >= $ymin && $cyf <= $ymax } {
                     incr i
                     if {[expr $i%$stride]} { continue }
                     set x [expr {$xplotmin + ($scalex*($cx-$xmin))}]
                     set y [expr {$yplotmin + ($scaley*($cy-$ymin))}]
                     if {[string match "point*" [lindex $marker $ds]]} {
                        set item [$c create oval [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 0 -fill $fc -tags "$nsets"]
                                  
                        $c addtag point withtag $item
                     } elseif {[lindex $marker $ds]=="circle"} {
                        set item [$c create oval [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 1 -outline $lc \
                                     -fill $fc -tags "$nsets"]
                        $c addtag point withtag $item
                     } elseif {[lindex $marker $ds]=="square"} {
                        set item [$c create rectangle [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 1 -outline $lc \
                                     -fill $fc -tags "$nsets"]
                        $c addtag point withtag $item
                     }
                  }
               }
            }

            # Draw the legend
            if {[llength $leg]} {
               variable ticfont
               variable ticfontsize
               set ylegpos [expr $yplotmax+2*$ticfontsize+$ds*2.2*$ticfontsize]
               set xlegpos [expr $xplotmin+30]
               set item [$c create line $xlegpos $ylegpos [expr $xlegpos+30] $ylegpos \
                  -width [lindex $linewidth $ds] -fill $lc -dash $dash]
               $c addtag legend withtag $item
               set item [$c create text [expr $xlegpos+30+$ticfontsize] $ylegpos -text $leg \
                            -font $ticfont -anchor w]
               $c addtag legend withtag $item
               if {[lindex $marker $ds]=="points"} {
                  set item [$c create oval [expr {$xlegpos-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+$rad}] [expr {$ylegpos+$rad}] -width 1 -fill $fc]
                  $c addtag legend withtag $item
                  set item [$c create oval [expr {$xlegpos+30-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+30+$rad}] [expr {$ylegpos+$rad}] -width 1 -fill $fc]
                  $c addtag legend withtag $item
               } elseif {[lindex $marker $ds]=="circle"} {
                  set item [$c create oval [expr {$xlegpos-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+$rad}] [expr {$ylegpos+$rad}] -width 1 -outline $lc \
                               -fill $fc]
                  $c addtag legend withtag $item
                  set item [$c create oval [expr {$xlegpos+30-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+30+$rad}] [expr {$ylegpos+$rad}] -width 1 -outline $lc \
                               -fill $fc]
                  $c addtag legend withtag $item
               }
            }

            incr ds
         }
         $c bind point <Any-ButtonPress> [namespace code {
            #$c itemconfig current -fill red;
            print_datapoint %x %y [$c itemcget current -fill]
         }]
         #$c bind point <Any-Leave> "$c itemconfig current -fill $fc"
         $c bind legend <1> "[namespace current]::grab_legend $c %x %y"
         $c bind legend <B1-Motion> "[namespace current]::move_legend $c %x %y"
      }


      proc calcIntersect { oldcx oldcy cx cy outofBounds } {
 
         variable xmin
         variable ymin
         variable xmax
         variable ymax
  
         set slope [expr {($cy - $oldcy)*1.0/($cx - $oldcx)}]
         if { $outofBounds % 2 == 0 } {
            set xinter $xmin
            set yinter [expr {$slope*($xmin - $oldcx) + $oldcy}]
         } elseif { $outofBounds % 3 == 0 } {
            set xinter $xmax
            set yinter [expr {$slope*($xmax - $oldcx) + $oldcy}]
         }
         if { $outofBounds % 5 == 0 } {
            if { $outofBounds % 2 == 0 || $outofBounds % 3 == 0 } {
               if { $yinter < $ymin } {
                  set yinter $ymin
                  set xinter [expr {($ymin - $oldcy)/$slope + $oldcx}]
               }
            } else {
               set yinter $ymin
               set xinter [expr {($ymin - $oldcy)/$slope + $oldcx}]
            }
         } elseif { $outofBounds % 7 == 0 } {
            if { $outofBounds % 2 == 0 || $outofBounds % 3 == 0 } {
               if { $yinter > $ymax } {
                  set yinter $ymax
                  set xinter [expr {($ymax - $oldcy)/$slope + $oldcx}]
               }
            } else {
               set yinter $ymax
               set xinter [expr {($ymax - $oldcy)/$slope + $oldcx}]
            }
         }
         return [list $xinter $yinter]
      }

      # Transforms coordinates from plot coordinates to canvas coords
      proc world2canvascoor {wx wy} {
         variable xplotmin
         variable yplotmin
         variable scalex
         variable scaley
         variable xmin
         variable ymin
         set x [expr {$xplotmin + ($scalex*($wx-$xmin))}]
         set y [expr {$yplotmin + ($scaley*($wy-$ymin))}]
         return [list $x $y]
      }                    
      
      proc redraw_objects {} {
        variable objectlist
        foreach object $objectlist {
           draw_object $object
        }
      }
      
      proc draw_object {object} {
        variable c
        set oname [lindex $object 0]
        set optpos [lsearch -regexp $object {^-[^[:digit:]]}]
        set options {}
        if {$optpos<0} {
          set optpos end
        } else {
          set options [lrange $object $optpos end]
          incr optpos -1
        }
        set coords [join [lrange $object 1 $optpos]]
        foreach {wx wy} $coords {
          lappend plotcoords [world2canvascoor $wx $wy]
        }
        if {$oname=="circle" || $oname=="square"} {
          set rad 1.0
          set pos [lsearch $options "-radius"]
          if {$pos>=0} {
             if {$pos+1<[llength $options]} {
                set rad [lindex $options [expr {$pos+1}]]
             }
             set options [lreplace $options $pos [expr {$pos+1}]]
          }
          foreach {x y} [join $plotcoords] {break}
          set plotcoords  [list [expr {$x-$rad}] [expr {$y-$rad}] [expr {$x+$rad}] [expr {$y+$rad}]]
          if {$oname=="circle"} { 
             set oname "oval"
          } else { set oname "rectangle" }
        }
        
        set evalstr "$c create $oname [join $plotcoords] $options"
        set item [eval $evalstr]
        $c addtag objects withtag $item
      }
      
      # grab_legend --
      # This procedure is invoked when the mouse is pressed over one of the
      # legend items.  It sets up state to allow the legend to be dragged.
      #
      # Arguments:
      # w -             The canvas window.
      # x, y -  The coordinates of the mouse press.

      proc grab_legend {w x y} {
         variable legendpos
         #$w dtag selected
         #$w addtag selected withtag current
         $w raise legend
         set legendpos(lastX) $x
         set legendpos(lastY) $y
      }

      # move_legend --
      # This procedure is invoked during mouse motion events.  It drags the
      # legend.
      #
      # Arguments:
      # w -             The canvas window.
      # x, y -  The coordinates of the mouse.

      proc move_legend {w x y} {
         variable legendpos
         $w move legend [expr {$x-$legendpos(lastX)}] [expr {$y-$legendpos(lastY)}]
         set legendpos(lastX) $x
         set legendpos(lastY) $y
      }

      proc calculate_range {} {
         # Get min/max values

         variable predefRange
         variable givenXmin
         variable givenYmin
         variable givenXmax
         variable givenYmax

         variable datasets
         set lxmin {}
         set lxmax {}
         set lymin {}
         set lymax {}
         foreach dsxmin $datasets(xmin) dsxmax $datasets(xmax) \
                 dsymin $datasets(ymin) dsymax $datasets(ymax) \
                 dsxmin_y $datasets(xmin_y) dsxmax_y $datasets(xmax_y) \
                 dsymin_x $datasets(ymin_x) dsymax_x $datasets(ymax_x) {
               lappend lxmin [list $dsxmin $dsxmin_y]
               lappend lymin [list $dsymin $dsymin_x]
               lappend lxmax [list $dsxmax $dsxmax_y]
               lappend lymax [list $dsymax $dsymax_x]
            }

         if { $predefRange } {
            if { $givenXmin == "auto" || $givenXmin == "Auto" } {
               set givenXmin [lindex [lsort -real -index 0 $lxmin] 0 0]
            }
            if { $givenXmax == "auto" || $givenXmax == "Auto" } {
               set givenXmax [lindex [lsort -real -index 0 $lxmax] end 0]
            }
            if { $givenYmin == "auto" || $givenYmin == "Auto" } {
               set givenYmin [lindex [lsort -real -index 0 $lymin] 0 0]
            }
            if { $givenYmax == "auto" || $givenYmax == "Auto" } {
               set givenYmax [lindex [lsort -real -index 0 $lymax] end 0]
            }
            if { $givenXmin < $givenXmax && $givenYmin < $givenYmax } {
               set tmpxmin $givenXmin
               set tmpymin $givenYmin
               set tmpxmax $givenXmax
               set tmpymax $givenYmax
            } else {
               variable predefRange 0
               set givenXmin auto
               set givenXmax auto
               set givenYmin auto
               set givenYmax auto
            }
         } 

         if { !$predefRange } {
            set tmpxmin [lindex [lsort -real -index 0 $lxmin] 0 0]
            set tmpymin [lindex [lsort -real -index 0 $lymin] 0 0]
            set tmpxmax [lindex [lsort -real -index 0 $lxmax] end 0]
            set tmpymax [lindex [lsort -real -index 0 $lymax] end 0]
         }

         variable xmin     
         variable ymin     
         variable xmax     
         variable ymax     
         if {$tmpxmin<$xmin || $tmpxmax>$xmax || $tmpymin<$ymin || $tmpymax>$ymax} {
            variable resize 1
         }

         variable xmin     [format "%10g" $tmpxmin]
         variable ymin     [format "%10g" $tmpymin]
         variable xmax     [format "%10g" $tmpxmax]
         variable ymax     [format "%10g" $tmpymax]
         variable spanx    [expr $xmax-$xmin]
         variable spany    [expr $ymax-$ymin]

         if { $predefRange } {
            variable xmin_y $ymin 
            variable ymin_x $xmin
            variable xmax_y $ymax
            variable ymax_x $xmax
         } else {
            variable xmin_y   [format "%10g" [lindex [lsort -real -index 0 $lxmin] 0 1]]
            variable ymin_x   [format "%10g" [lindex [lsort -real -index 0 $lymin] 0 1]]
            variable xmax_y   [format "%10g" [lindex [lsort -real -index 0 $lxmax] end 1]]
            variable ymax_x   [format "%10g" [lindex [lsort -real -index 0 $lymax] end 1]]
         }

         # Order of magnitude of value range
         if {$spanx==0.0} { variable spanx 1 } 
         if {$spany==0.0} { variable spany 1 } 
         variable dimx [expr 0.5*pow(10,floor(log10($spanx)))]
         variable dimy [expr 0.5*pow(10,floor(log10($spany)))]
      }
         
      proc calculate_ticspacing {} {
         variable spanx
         variable spany
         variable dimx
         variable dimy

         # Total number of tics between two major tics
         variable minorticx 5
         if {[expr $spanx/$dimx]>5} {
            variable minorticx 2
         }
         
         variable minorticy 5
         if {[expr $spany/$dimy]>5} {
            variable minorticy 2
         }

         variable anglescalex
         variable anglescaley
         if {$anglescalex} {
            set dimx 90
            set minorticx 3
         }
         if {$anglescaley} {
            set dimy 90
            set minorticy 3
         }

         variable xmajortics
         variable ymajortics
         variable xminortics
         variable yminortics
         if {[llength $xmajortics]} { set dimx $xmajortics }
         if {[llength $ymajortics]} { set dimy $ymajortics }
         if {[llength $xminortics]} { set minorticx $xminortics }
         if {[llength $yminortics]} { set minorticy $yminortics }
#        set i 0
#        while {1} {
#           variable loticx [expr $i*$minorticx]
#           if {$loticx<$xmin} { return [expr $i*$minorticx]}
#           incr i
#        }
         variable xmin
         variable ymin
         variable xmax
         variable ymax
         if {${::nmwiz_MultiPlot::verbose}} {
           vmdcon -info "dimx=$dimx xmin=$xmin xmax=$xmax ceil=[expr ceil($xmin/$dimx*$minorticx)]"
           vmdcon -info "dimy=$dimy ymin=$ymin ymax=$ymax ceil=[expr ceil($ymin/$dimy*$minorticy)]"
         }
         variable loticx [expr $dimx*ceil($xmin/$dimx*$minorticx)/$minorticx]
         variable loticy [expr $dimy*ceil($ymin/$dimy*$minorticy)/$minorticy]
         variable hiticx [expr $dimx*floor($xmax/$dimx*$minorticx)/$minorticx]
         variable hiticy [expr $dimy*floor($ymax/$dimy*$minorticy)/$minorticy]
      }

      proc calculate_labelsize {} {
         # Measure y-axis label size
         variable c
         variable labelfont
         variable xlabeltext
         variable ylabeltext
         if {[llength $ylabeltext]} {
            set item [$c create text 0 0 -text $ylabeltext -font $labelfont -anchor nw]
            set bbox [$c bbox $item]
            variable ylabelheight [expr [lindex $bbox 3]-[lindex $bbox 1]]
            variable ylabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0] + $ylabelheight]
            $c delete $item
         } else {
            variable ylabelwidth 0.0
         }

         # Measure x-axis label height
         if {[llength $xlabeltext]} {
            set item [$c create text 0 0 -text $xlabeltext -font $labelfont -anchor nw]
            set bbox [$c bbox $item]
            $c delete $item
            variable xlabelheight  [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         } else {
            variable xlabelheight 0.0
         }

         
         ## Measure x-axis ticlabel size
         variable loticx
         variable hiticx
         variable ticfont
         # Compare smallest and biggest tics
         set absxmax [lindex [lsort -real [list [expr abs($loticx)] [expr abs($hiticx)]]] end]
         set item [$c create text 0 0 -text [format "-%g" $absxmax] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable ticlabelheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         variable xticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]

         ## Measure y-axis ticlabel size
         variable dimx
         variable loticy
         variable hiticy
         variable ticfont
         # Compare smallest and biggest tics
         set absymax [lindex [lsort -real [list [expr abs($loticy)] [expr abs($hiticy)]]] end]
         set item [$c create text 0 0 -text [format "-%g" $absymax] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable ticlabelheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         variable yticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]
         # Check if the neighboring ticlabel is wider since it could involve more decimal places
         set item [$c create text 0 0 -text [format "-%g" [expr $absymax+$dimx]] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         if {[expr 1.5*[lindex $bbox 2]-[lindex $bbox 0]]>$yticlabelwidth} {
            variable yticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]
         }

         # Measure title height
         variable title
         variable titlefont
         if {![llength $title]} { set title [namespace current] }
         set item [$c create text 0 0 -text $title -font $titlefont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable titleheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
      } 
      
      proc calculate_plot_geometry {} {
         # Compute legend height
         variable legend         
         variable ticfontsize
         variable legendheight 0.0
         foreach legitem $legend {
            if {[llength $legitem]} {
               set legendheight [expr $legendheight+1.8*$ticfontsize]
            }
         }
         
         ## Plot geometry
         variable rim
         variable canh
         variable canw
         variable ticlen
         variable xlabelheight
         variable ylabelwidth
         variable xticlabelwidth
         variable yticlabelwidth
         variable ticlabelheight
         variable titleheight
         variable xplotmin [expr $rim+$ylabelwidth+$yticlabelwidth+$ticlen]
         variable yplotmin [expr $canh-($rim+$xlabelheight+$ticlabelheight+$ticlen)]
         variable xplotmax [expr $canw-$rim-0.5*$xticlabelwidth]
         variable yplotmax [expr $rim+$titleheight]

         # Scaling factor to convert world coordinates into plot coordinates
         variable spanx
         variable spany
         variable scalex [expr ($xplotmax-$xplotmin)/(1.0*$spanx)]      
         variable scaley [expr ($yplotmax-$yplotmin)/(1.0*$spany)]      

         variable dimx
         variable scalex
         if {[expr $xticlabelwidth]>[expr $dimx*$scalex*0.7]} {
            set dimx [expr 2.0*$dimx]
            calculate_ticspacing
            calculate_labelsize
            calculate_plot_geometry
         }
      }

      proc draw_periphery {} {
         # Draw title
         variable c
         variable rim
         variable canw
         variable canh
         variable title
         variable titlefont
         $c create text [expr $canw/2] $rim -anchor n -text $title -font $titlefont -fill brown

         # Draw bounding box
         variable xplotmin
         variable yplotmin
         variable xplotmax
         variable yplotmax
         $c create line $xplotmin $yplotmin $xplotmax $yplotmin -width 2
         $c create line $xplotmin $yplotmin $xplotmin $yplotmax -width 2
         $c create line $xplotmax $yplotmin $xplotmax $yplotmax -width 2
         $c create line $xplotmin $yplotmax $xplotmax $yplotmax -width 2

         # Scaling factor to convert plot coordinates into canvas coordinates
         variable spanx
         variable spany
         variable scalex        
         variable scaley
         variable xmin  
         variable ymin  
         variable xmax  
         variable ymax  

         # x-axis, y=0
         if {$ymin<0 && $ymax>0} {
            set zero [expr $yplotmin-($scaley*$ymin)]
            $c create line $xplotmin $zero $xplotmax $zero -width 1 -dash -
         }
         # y-axis, x=0
         if {$xmin<0 && $xmax>0} {
            set zero [expr $xplotmin-($scalex*$xmin)]
            $c create line $zero $yplotmin $zero $yplotmax -width 1 -dash -
         }

         # x-label
         variable ticlen
         variable labelfont
         variable xlabeltext
         variable xlabelheight
         if {[llength $xlabeltext]} {
            set labelposx [expr $xplotmin+($xplotmax-$xplotmin)*0.5]
            #set labelposy [expr $yplotmin+$ticlen+$ticlabelheight+0.2*$xlabelheight]
            set labelposy [expr $canh-$rim-0.2*$xlabelheight]
            $c create text $labelposx $labelposy -text $xlabeltext -font $labelfont -anchor s
         }

         # y-label
         variable ylabeltext
         if {[llength $ylabeltext]} {
            set labelposy [expr $yplotmin+($yplotmax-$yplotmin)*0.5]
            $c create text $rim $labelposy -text $ylabeltext -font $labelfont -anchor w
         }

         # Draw x-tics
         variable ticfont
         set i 0 
         set ticval $xmin
         variable dimx
         variable hiticx
         variable loticx
         variable minorticx
         variable ticlabelheight
         set firstmajor [expr abs(int($loticx-int($loticx/$minorticx)*$minorticx))]
         #set firstmajor [expr $loticx%$minorticx]
         set ticlabelposy [expr $yplotmin+$ticlen+0.2*$ticlabelheight]
         while {$ticval<$hiticx} {
            set ticval [expr $loticx+$i*$dimx/$minorticx]
            set x [expr $xplotmin + ($ticval-$xmin)*$scalex]
            if {![expr ($i-$firstmajor)%$minorticx]} {
               $c create line $x $yplotmin $x [expr $yplotmin+$ticlen] -width 2
               $c create text $x $ticlabelposy -text [format "%g" $ticval] -anchor n -font $ticfont
            } else {
               $c create line $x $yplotmin $x [expr $yplotmin+0.5*$ticlen] -width 2
            }
            incr i
         }

         # Draw y-tics
         set i 0
         variable dimy
         variable hiticy
         variable loticy
         variable minorticy
         set firstmajor [expr abs(int($loticy-int($loticx/$minorticy)*$minorticy))]
         set ticlabelposx [expr $xplotmin-$ticlen-0.2*$ticlabelheight]
         set ticval $ymin
         while {$ticval<$hiticy} {
            set ticval [expr $loticy+$i*$dimy/$minorticy]
            set y [expr $yplotmin + ($ticval-$ymin)*$scaley]
            if {![expr ($i-$firstmajor)%$minorticy]} {
               $c create line $xplotmin $y [expr $xplotmin-$ticlen] $y -width 2
               $c create text $ticlabelposx $y -text [format "%g" $ticval] -anchor e -font $ticfont
            } else {
               $c create line $xplotmin $y [expr $xplotmin-0.5*$ticlen] $y -width 2
            }
            incr i
         }

         # Draw user specified horizontal lines
         variable hline
         foreach line $hline {
            set y [lindex $line 0]
            set zero [expr $yplotmin-($scaley*$ymin)]
            set opt [lrange $line 1 end]
            set ypos [expr $yplotmin+($scaley*($y-$ymin))]
            if {${::nmwiz_MultiPlot::verbose}} {
              vmdcon -info "$c create line $xplotmin $ypos $xplotmax $ypos $opt"
            }
            eval $c create line $xplotmin $ypos $xplotmax $ypos $opt
         }

         # Draw user specified vertical lines
         variable vline
         foreach line $vline {
            set x [lindex $line 0]
            set opt [lrange $line 1 end]
            set xpos [expr $xplotmin+($scalex*($x-$xmin))]
            eval $c create line $xpos $yplotmin $xpos $yplotmax $opt
         }
      }

      proc print_datapoint {x y clr} {
         variable xplotmin
         variable yplotmin
         variable scalex
         variable scaley
         variable xmin
         variable ymin
         variable nmwizns
         variable fillcolor
         #set coords [format "%8g %8g" [expr ($x-$xplotmin)/$scalex+$xmin] [expr ($y-$yplotmin)/$scaley+$ymin]]
         ${nmwizns}::selectAtom [::tcl::mathfunc::int [expr ($x-$xplotmin)/$scalex+$xmin]] $clr  
         #puts $coords
      }

      proc draw_infobox {} {
         variable w
         variable infoFont
         labelframe $w.info -text "Info"
         label $w.info.headx -text [format "%10s" "x"] -font $infoFont
         label $w.info.heady -text [format "%10s" "y"] -font $infoFont
         grid $w.info.headx -row 1 -column 2
         grid $w.info.heady -row 1 -column 3
         label $w.info.xmint -text "X min: "       -font $infoFont
         label $w.info.xmin  -textvariable [namespace current]::xmin   -font $infoFont
         label $w.info.xminy -textvariable [namespace current]::xmin_y -font $infoFont
         grid $w.info.xmint  -row 2 -column 1
         grid $w.info.xmin   -row 2 -column 2
         grid $w.info.xminy  -row 2 -column 3
         label $w.info.xmaxt -text "X max: "       -font $infoFont
         label $w.info.xmax  -textvariable [namespace current]::xmax   -font $infoFont
         label $w.info.xmaxy -textvariable [namespace current]::xmax_y -font $infoFont
         grid $w.info.xmaxt  -row 3 -column 1
         grid $w.info.xmax   -row 3 -column 2
         grid $w.info.xmaxy  -row 3 -column 3
         label $w.info.ymint -text "Y min: "       -font $infoFont
         label $w.info.ymin  -textvariable [namespace current]::ymin   -font $infoFont
         label $w.info.yminx -textvariable [namespace current]::ymin_x -font $infoFont
         grid $w.info.ymint  -row 4 -column 1
         grid $w.info.ymin   -row 4 -column 2
         grid $w.info.yminx  -row 4 -column 3
         label $w.info.ymaxt -text "Y max: "       -font $infoFont
         label $w.info.ymax  -textvariable [namespace current]::ymax   -font $infoFont
         label $w.info.ymaxx -textvariable [namespace current]::ymax_x -font $infoFont
         grid $w.info.ymaxt  -row 5 -column 1
         grid $w.info.ymax   -row 5 -column 2
         grid $w.info.ymaxx  -row 5 -column 3
         pack  $w.info -side topm -ipadx 5m -ipady 2m
      }

      proc savedialog {} {
         variable w
         set types {
            {{Postscript files} {.ps}}
            {{All files}        *   }
         }
         variable postscript
         set newfile [tk_getSaveFile \
                             -title "Choose file name" -parent $w \
                             -initialdir [pwd] -filetypes $types -initialfile $postscript]
         variable c
         if {[llength $newfile]} {
            $c postscript -file $newfile
            vmdcon -info "Wrote plot postscript file to $newfile."
         }
      }

      proc xmgracedialog {} {
         variable w
         set types {
            {{xmgrace files} {.agr}}
            {{All files}        *  }
         }
         set newfile [tk_getSaveFile \
                      -title "Choose file name" -parent $w \
                     -initialdir [pwd] -filetypes $types -initialfile "multiplot.agr"]
         if {[llength $newfile]} {
           [namespace current]::plothandle export xmgrace $newfile
           vmdcon -info "Wrote plot to $newfile, running xmgrace..."
           ::ExecTool::exec xmgrace $newfile &
         }
      }

   } ; # END namespace $ns

   return "::nmwiz_MultiPlot::Plot${::nmwiz_MultiPlot::plotcount}::plothandle"
}

proc nmwiz_multiplot { args } {
   set keyword [lindex $args 0]
   if {![llength $keyword]} { return }
   if {$keyword=="list"} {
      set plist {}
      foreach plot [namespace children ::nmwiz_MultiPlot "Plot*"] { 
         lappend plist [subst $plot]::plothandle
      }
      return $plist
   } elseif {$keyword=="reset"} {
      # on reset we may delete only toplevel widgets
      foreach ploth [namespace children ::nmwiz_MultiPlot "Plot*"] {
         eval "set tl $[subst $ploth]::istoplevel"
         if { $tl } {
            if {$::nmwiz_MultiPlot::verbose} {
              vmdcon -info "deleting toplevel widget of $ploth"
            }
            eval "destroy $[subst $ploth]::w"
            namespace delete $ploth
         }
      }
      return
   }

   variable plothandle
   if {$keyword=="embed"} {
     set parent [lindex $args 1]
     if {$parent == ""} {
        return
     } else {
        set plothandle [::nmwiz_MultiPlot::init_plot $parent]
     }
   } else {
     set plothandle [::nmwiz_MultiPlot::init_plot ]
   } 
   #puts "$plothandle configure $args"
   eval $plothandle configure $args

   return $plothandle
}

#::nmwiz_MultiPlot::init_plot

# List of NMWiz functions
#   ::nmwiz::initGUI                 make main window
#   ::nmwiz::loadNMD                 load NMD file

#   ::nmguiX::deleteMolecules        delete molecules storing coordinates/graphics/animations/selection
#   ::nmguiX::prepareSelmol          prepare the molecule where selections are displayed
#   ::nmguiX::getPDBLines            return PDB files for the molecule
#   ::nmguiX::clearSelection         turn selection representation of and clear labels


namespace eval ::nmwiz:: {
  namespace export nmwizgui
  namespace export initialize

  variable guicount -1
  variable tmpdir
  variable titles [list]
  variable plothandles [list] 
  #variable namespaces [list]
  #variable nmwizguis [list]
  variable platform $tcl_platform(platform) 
  switch $platform {
    unix {
      set tmpdir "/tmp" ;  # or even $::env(TMPDIR), at times.
    } macintosh {
      set tmpdir $::env(TRASH_FOLDER)  ;# a better place?
    } default {
      set tmpdir [pwd]
      catch {set tmpdir $::env(TMP)}
      catch {set tmpdir $::env(TEMP)}
    }
  }
  
  proc showHelp {context} {
    set windowname nmwizhelp
    if {[winfo exists .$windowname] == 0} {
      set log [toplevel ".$windowname"]
      wm title $log "NMWiz Help"
      wm resizable $log 1 1
      incr logcount

      text $log.text -bg White -bd 2 -font Courier \
        -yscrollcommand ".$windowname.vscr set"
      scrollbar $log.vscr -command ".$windowname.text yview"
      pack $log.text -side left -fill both -expand 1  
      pack $log.vscr -side right -fill y
    } else {
      set log .$windowname
    }
    $log.text configure -state normal -wrap word
    $log.text delete 1.0 end
      
    if {$context == "wizard"} {
      $log.text insert end "ProDy Interface\n"
      $log.text insert end "===============\n\n"
      $log.text insert end "\nActive Mode\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "Select the active mode for which you want to draw arrows or make an animation. "
      $log.text insert end "Direction of arrows depicting the normal mode can be changed using +/- button. "
      $log.text insert end "Arrows can be drawn along both directions by changing the options Arrow Graphics Options panel. "
      $log.text insert end "The selected color effects both arrow graphics and square fluctuation plots."
      $log.text insert end "\n\n**RMSD**\n\n"
      $log.text insert end "The RMSD corresponding to the displacement described by the arrows is displayed. User can change the RMSD value to rescale the arrows. "
      $log.text insert end "The scaling factor that produces the specified RMSD is printed to the VMD console (along with the magnitude of the mode provided in NMD file). "
      $log.text insert end "\n\n**Selection**\n\n"
      $log.text insert end "Selection entry allows the user to display arrows for a subset of atoms.\n\n"
      $log.text insert end "*TIP*: If the arrow graphics are too crowded or the display is slow, draw arrows for an evenly spaced subset of residues, e.g try 'name CA and residue % 4 == 0', which will draw an arrow for every fourth residue."
      $log.text insert end "\n\n\n"
      $log.text insert end "Arrow Graphics\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "Id of the molecule that contains the arrow graphics of the active mode is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Draw: draw/redraw arrow graphics for the active mode\n"
      $log.text insert end " * Clean: remove most recently drawn arrow graphics\n"
      $log.text insert end " * Hide/Show: hide/show most recently drawn arrow graphics\n"
      $log.text insert end " * Options: show/hide arrow graphics option panel\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can change arrow graphics properties and how NMWiz behave upon such changes in this panel.\n"
      $log.text insert end "\nBy default:\n\n"
      $log.text insert end " * arrow graphics are set to automatically change when graphics properties are changed by the user\n"
      $log.text insert end " * current graphics are hidden the active mode is changed\n"
      $log.text insert end "\nOptionally:\n\n"
      $log.text insert end " * arrows can be drawn in both directions to look like a double headed arrow\n"
      $log.text insert end " * arrows shorter than a length (A) can be hidden\n"
      $log.text insert end "\nAdditionally, user can change:\n\n"
      $log.text insert end " * width of the arrow cylinder\n"
      $log.text insert end " * width/height of the arrow head code\n"
      $log.text insert end " * graphics material and resolution"
      $log.text insert end "\n\n\n"
      $log.text insert end "Animations\n"
      $log.text insert end "----------\n\n"
      $log.text insert end "Id of the molecule that contains the most recently generated animation is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Draw: animate fluctuations along the active mode\n"
      $log.text insert end " * Play : play/pause the animation\n"
      $log.text insert end " * Hide : hide/show the animation\n"
      $log.text insert end " * Options: show/hide animation option panel\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can elect automatic generation and continuous play of animations when the active mode changes. User can also select the number of frames in the animation."
      $log.text insert end "\n\n\n"
      $log.text insert end "Plotting\n"
      $log.text insert end "--------\n\n"
      $log.text insert end "Id of the molecule for displaying selected residues is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Plot: plot squared-fluctuations along the active mode\n"
      $log.text insert end " * Clear: clear all selections and selected atom labels\n"
      $log.text insert end " * Hide/Show: hide/show the selected residues\n"
      $log.text insert end " * Options: change plotting options\n"
      $log.text insert end "\n\n\n"
      $log.text insert end "Protein Graphics\n"
      $log.text insert end "----------------\n\n"
      $log.text insert end "Id of the molecule that contains the protein structure is shown in parentheses.\n\n"
      $log.text insert end "Buttons:\n\n"
      $log.text insert end " * Update: Uupdate protein representation\n"
      $log.text insert end " * Focus: reset view to focus on the structure\n"
      $log.text insert end " * Hide/Show: hide/show strudture\n"
      $log.text insert end " * Options: change molecular system representation\n"
      $log.text insert end "\nOptions:\n\n"
      $log.text insert end "User can select the representation and coloring scheme. User can change the protein representation settings manually, by setting 'Show protein as' to 'Custom'.\n\n"      
      $log.text insert end "Protein can be colored based on the `Mobility` of the residues in the active mode, based on 'Bfactors' that came in NMD file, or based on residue/atom 'Index'.\n\n"
      $log.text insert end "In addition to the standard representations (Tube/Trace/Licorice), protein can be represented as an elastic network."
      $log.text insert end "User can set the cutoff distance, width of dynamic bonds, and node spheres. Note that changing the cutoff distance distance only affects representation, not the precalculated normal mode data.\n\n"
      $log.text insert end "*TIP*: When visualizing a large system, display protein at lower resolutions and/or try displaying fewer atoms if all atoms are displayed."
    } elseif {$context == "prody"} {
      $log.text insert end "ProDy Interface\n"
      $log.text insert end "===============\n\n"
      $log.text insert end "ProDy interface allows users to perform the following calculations for molecules loaded in VMD:\n\n"
      $log.text insert end "* Anisotropic Network Model (ANM)\n"
      $log.text insert end "* Gaussian Network Model (GNM)\n"
      $log.text insert end "* Principal Component Analysis (PCA) a.k.a. Essential Dynamics Analysis (EDA)\n\n\n"
      $log.text insert end "Atom Selection\n"
      $log.text insert end "--------------\n\n"
      $log.text insert end "First thing you need to do is selecting the molecule and specifying the atoms that you want to include in the calculations. "
      $log.text insert end "If you do not see all molecules in the menu, click 'Update'."
      $log.text insert end "\n\n\n"
      $log.text insert end "ProDy Job Settings\n"
      $log.text insert end "------------------\n\n"
      $log.text insert end "Specify the calculation type and output options in this panel. "
      $log.text insert end "Coordinate data for selected atoms and the NMD data after calculations will be written into the 'Output directory'. "
      $log.text insert end "All output files will named after 'Output filename'.\n\n"
      $log.text insert end "**ProDy Scripts**\n\n"
      $log.text insert end "Note that you need to specify the path to the individual ProDy scripts in 'Settings' NMWiz will save and reload the path in the future sessions."
      $log.text insert end "\n\n\n"
      $log.text insert end "ANM/GNM Settings\n"
      $log.text insert end "----------------\n\n"
      $log.text insert end "Specify the following:\n\n"
      $log.text insert end " * number of modes to be calculated\n"
      $log.text insert end " * index of the frame (coordinate set) to be used in calculations\n" 
      $log.text insert end " * cutoff distance\n" 
      $log.text insert end " * force constant"
      $log.text insert end "\n\n\n"
      $log.text insert end "PCA/EDA Settings\n"
      $log.text insert end "----------------\n\n"
      $log.text insert end "Note that for PCA/EDA calculations molecule must have multiple frames. Specify the range of frames to be used in calculations. For large systems, prefer to write coordinates in DCD format to gain IO speed and save save disk space."
    } elseif {$context == "compare"} {
      $log.text insert end "Structure Comparison\n"
      $log.text insert end "====================\n\n"
      $log.text insert end "This interface allows the user to compare two molecules (or two frames of the same molecule) loaded in VMD. It can be used to align structures and draw deformation vector."
      $log.text insert end "Follow these steps to compare structures:\n\n"
      $log.text insert end "1) Load molecules into VMD\n\n"
      $log.text insert end "Select the molecules that you want to compare. If you don't see the molecule, click 'Update' button. "
      $log.text insert end "\n\n"
      $log.text insert end "2) Select atoms and specify frames\n\n"
      $log.text insert end "Make atom selections for each molecule (or frame). Note that the number of selected atoms must be the same. "
      $log.text insert end "\n\n"
      $log.text insert end "3) Align molecules (or frames)\n\n"
      $log.text insert end "Before deformation vector is calculated, you need to align the molecules (or frames) to have a meaningful depiction of structural change. "
      $log.text insert end "\n\n"
      $log.text insert end "Finally, click 'Calculate' button to generate depiction of structural change and NMWiz GUI."
    } elseif {$context == "frommolecule"} {
      $log.text insert end "From Molecule\n"
      $log.text insert end "=============\n\n"
      $log.text insert end "This interface allows the user to analyze normal mode data present in file formats that are recognized by VMD. "
      $log.text insert end "Follow these steps to analyze your data:\n\n"
      $log.text insert end "1) Load data into VMD\n\n"
      $log.text insert end "Normal mode data can be retrieved from a molecule with multiple frames in VMD. "
      $log.text insert end "The molecule must contain both coordinate and normal mode data. "
      $log.text insert end "First, you need to load the coordinate data as a new molecule and then load the normal mode data into the same molecule."
      $log.text insert end "\n\n"
      $log.text insert end "2) Select the molecule and atoms\n\n"
      $log.text insert end "Select the molecule with normal mode data. If you don't see the molecule, click 'Update' button. "
      $log.text insert end "You can also select a subset of atoms for which you want to display normal mode data graphics."
      $log.text insert end "\n\n"
      $log.text insert end "3) Specify data frames\n\n"
      $log.text insert end "Frames that contain coordinate and normal mode data must be specified. "
      $log.text insert end "Note that '0' is the index of the very first frame, and 'end' can be used to specify the last frame of the molecule. "
      $log.text insert end "\n\n"
      $log.text insert end "Finally, click load button to instantiate NMWiz window for selected data."
      $log.text insert end "\n\n\n"
      $log.text insert end "*TIP*: If normal mode data is calculated for all atom data for a large molecular system, select backbone or carbon alpha atoms "
      $log.text insert end "for more responsive visual analysis experience, e.g enter 'name CA' as the selection string."
    } elseif {$context == "main"} {
      $log.text insert end "NMWiz Main\n"
      $log.text insert end "==========\n\n\n"
      $log.text insert end "Load Normal Mode Data\n"
      $log.text insert end "---------------------\n\n"
      $log.text insert end "Main interface allows user to load data into NMWiz in two ways: "
      $log.text insert end "\n\n"
      $log.text insert end "**Load NMD File**\n\n"
      $log.text insert end "If you have an NMD file, click 'Load NMD file' button to select the file. "
      $log.text insert end "The contents will be loaded and a Wizard window associated with the data will appear."
      $log.text insert end "\n\n"
      $log.text insert end "**From Molecule**\n\n"
      $log.text insert end "Alternatively, when normal mode data is present in a file format that is recognized by VMD, "
      $log.text insert end "load the files into VMD as a molecule and click 'From Molecule' button. A window will appear to "
      $log.text insert end "facilitate selection of normal mode data from a molecule with multiple frames.\n\n"
      $log.text insert end "Note that, data obtained from a molecule can be saved in NMD format from the Main window and NMD files can be parsed with ProDy for further analysis."
      $log.text insert end "\n\n\n"
      $log.text insert end "Perform NMA Calculations\n"
      $log.text insert end "------------------------\n\n"
      $log.text insert end "You can use NMWiz to perform NMA calculations via ProDy for molecules loaded in VMD. "
      $log.text insert end "Click 'ProDy Interface' and follow the instructions therein for ANM, GNM, and PCA (EDA) calculations."
      $log.text insert end "\n\n\n"
      $log.text insert end "Settings and Options\n"
      $log.text insert end "--------------------\n\n"
      $log.text insert end "NMWiz saves some user settings in your home folder. These settings can be changed using 'Settings' window."
      $log.text insert end "\n\n"
      $log.text insert end "**Preserve View**\n\n"
      $log.text insert end "When NMWiz loads data, VMD will shift focus to the new molecule. Check this to preserve the current view when loading a new dataset."
    }
      
    $log.text yview moveto 0
    $log.text configure -state disabled
  }
  
  # Called by nmwiz_tk function
  # Makes the Main Window
  proc initGUI {} {
    variable w
    variable platform  
    if [winfo exists .nmwizgui] {
      wm deiconify .nmwizgui
      raise .nmwizgui
      return 
    }
    ::nmwiz::loadSettings
    set w [toplevel .nmwizgui]
    wm title $w "NMWiz 0.8 - Main"
    wm resizable $w 0 0

    set wmf [frame $w.mainframe -bd 2]
    
    grid [button $wmf.loadnmd -width 20 -text "Load NMD File" -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"NMD files" { .nmd .NMD }} {"Text files" { .txt .TXT }} {"All files" *}}]
        if {![string equal $tempfile ""]} {::nmwiz::loadNMD $tempfile}}] \
      -row 3 -column 0 -columnspan 3 -sticky we

    grid [button $wmf.fromol -width 20 -text "From Molecule" -command ::nmwiz::initFromMolecule] \
      -row 5 -column 0 -columnspan 3 -sticky we

    grid [button $wmf.prody -width 20 -text "ProDy Interface" -command ::nmwiz::initProdyGUI] \
      -row 6 -column 0 -columnspan 3 -sticky we

    grid [button $wmf.compare -width 20 -text "Structure Comparison" -command ::nmwiz::initStrComp] \
      -row 7 -column 0 -columnspan 3 -sticky we

   
    grid [button $wmf.showhelp -text "Help" \
        -command {::nmwiz::showHelp main}] \
      -row 8 -column 0 -sticky we
    grid [button $wmf.settings -text "Settings" \
        -command ::nmwiz::initSettingsGUI] \
      -row 8 -column 1 -sticky we
    grid [button $wmf.website -text "Website" \
        -command "vmd_open_url http://www.csb.pitt.edu/NMWiz/"] \
      -row 8 -column 2 -sticky we
    
    if {[molinfo num] > 0} {
      set ::nmwiz::preserview 1
    }
    
    grid [checkbutton $wmf.preserview -text "preserve current view" \
        -variable ::nmwiz::preserview] \
      -row 10 -column 0 -columnspan 3 -sticky w

    #pack $wmf.options -side top -fill x -expand 1
    pack $wmf -side top -fill x -expand 1

    if {$::nmwiz::guicount > -1} {    
      for {set i 0} {$i <= $::nmwiz::guicount} {incr i} {
        set ns "::nmgui$i"
        if {[namespace exists $ns]} {      
          set wgf [labelframe $w.{[string range $ns 2 end]}frame -text "[subst $${ns}::title]" -bd 2]
          grid [button $wgf.show -text "Show GUI" \
              -command "${ns}::nmwizgui" ] \
            -row 0 -column 0 -sticky we
          grid [button $wgf.remove -text "Remove" \
              -command "lset ::nmwiz::titles $::nmwiz::guicount NONE; pack forget $wgf; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
            -row 0 -column 1 -sticky we
          grid [button $wgf.save -text "Save" \
              -command "::nmwiz::writeNMD $ns"] \
            -row 0 -column 2 -sticky we
          pack $wgf -side top -fill x -expand 1
        }
      }
    }
  }
  
  
  variable nmwizColors "blue red gray orange yellow tan green white pink \
cyan purple black yellow2 yellow3 green2 green3 \
cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
orange3"
  switch [vmdinfo arch] {
    WIN64 -
    WIN32 {
      variable pybin "[::ExecTool::find python.exe]"
      variable pyANM ""
      variable pyPCA ""
      variable pyGNM ""
    }
    default {
      variable pybin "[::ExecTool::find python]"
      variable pyANM "[::ExecTool::find anm.py]"
      variable pyPCA "[::ExecTool::find pca.py]"
      variable pyGNM "[::ExecTool::find gnm.py]"
    }
  }
  variable outputdir [pwd]
  variable defaultColor "purple"
  variable settings [dict create anm $pyANM gnm $pyGNM pca $pyPCA color $defaultColor outputdir $outputdir pybin $pybin]
  proc saveSettings {} {
    vmdcon -info "Saving NMWiz settings"
    dict set ::nmwiz::settings anm $::nmwiz::pyANM
    dict set ::nmwiz::settings gnm $::nmwiz::pyGNM
    dict set ::nmwiz::settings pca $::nmwiz::pyPCA
    dict set ::nmwiz::settings color $::nmwiz::defaultColor
    dict set ::nmwiz::settings outputdir $::nmwiz::outputdir
    dict set ::nmwiz::settings pybin $::nmwiz::pybin
    global env
    set outf [open [file join $env(HOME) .nmwiz] w]
    puts $outf "$::nmwiz::settings"
    close $outf
  }
  
  proc loadSettings {} {
    vmdcon -info "Loading NMWiz settings"
    global env
    if {[file exists [file join $env(HOME) .nmwiz]]} {
      set inf [open [file join $env(HOME) .nmwiz]]
      gets $inf line
      close $inf
      foreach {key value} [split $line] {
        if {[dict exists $::nmwiz::settings $key]} {
          dict set ::nmwiz::settings $key $value
        }
      }
      variable ::nmwiz::pyANM "[dict get $::nmwiz::settings anm]"
      variable ::nmwiz::pyGNM "[dict get $::nmwiz::settings gnm]"
      variable ::nmwiz::pyPCA "[dict get $::nmwiz::settings pca]"
      variable ::nmwiz::defaultColor [dict get $::nmwiz::settings color]
      variable ::nmwiz::outputdir "[dict get $::nmwiz::settings outputdir]"
      variable ::nmwiz::pybin "[dict get $::nmwiz::settings pybin]"
    }
  }
  

  proc initSettingsGUI {} {
    variable settingsGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizsettings] {
      wm deiconify .nmwizsettings
      raise .nmwizsettings
      return 
    }    
    set settingsGUI [toplevel .nmwizsettings]
    wm title $settingsGUI "NMWiz - Settings"
    wm resizable $settingsGUI 0 0
    
    set wf [labelframe $settingsGUI.mainFrame -text "NMWiz Settings" -bd 2]

    grid [button $wf.clrHelp -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "The default color for arrow graphics."}] \
      -row 0 -column 0 -sticky w
    grid [label $wf.scriptLabel -text "Default color:"] \
      -row 0 -column 1 -sticky w
    grid [frame $wf.colorFrame] \
      -row 0 -column 2 -sticky ew
    tk_optionMenu $wf.colorFrame.list ::nmwiz::defaultColor "" 
    $wf.colorFrame.list.menu delete 0 last
    foreach acolor $::nmwiz::nmwizColors {
      $wf.colorFrame.list.menu add radiobutton -label $acolor \
          -variable ::nmwiz::defaultColor \
          -command ::nmwiz::saveSettings
    }
    pack $wf.colorFrame.list -side left -anchor w -fill x

    grid [button $wf.pyHelp -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Specify the path to the Python executable. If the folder\
                    that contains \"python\" or \"python.exe\"\
                    is included in your environment variable PATH, you may\
                    keep this as \"python\", otherwise\
                    specify the path to the executable,\
                    e.g. \"C:\\python27\\python.exe\""}] \
      -row 1 -column 0 -sticky w
    grid [label $wf.pyLabel -text "Python:"] \
      -row 1 -column 1 -sticky w
    grid [entry $wf.pyEntry -width 20 -textvariable ::nmwiz::pybin] \
      -row 1 -column 2 -sticky ew
    grid [button $wf.pyBrowse -text "Browse" \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"All files" *}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::python $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 1 -column 3 -sticky ew
      
    grid [button $wf.anmHelp -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Full path to ProDy ANM script (anm.py),\
                    e.g. C:\\python27\\Scripts\\anm.py"}] \
      -row 3 -column 0 -sticky w
    grid [label $wf.anmLabel -text "ANM script:"] \
      -row 3 -column 1 -sticky w
    grid [entry $wf.anmEntry -width 20 -textvariable ::nmwiz::pyANM] \
      -row 3 -column 2 -sticky ew
    grid [button $wf.anmBrowse -text "Browse" \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"ANM Script" { anm.py }}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::pyANM $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 3 -column 3 -sticky ew
      
    grid [button $wf.gnmHelp -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Full path to ProDy GNM script (gnm.py),\
                    e.g. C:\\python27\\Scripts\\gnm.py"}] \
      -row 4 -column 0 -sticky w
    grid [label $wf.gnmLabel -text "GNM script:"] \
      -row 4 -column 1 -sticky w
    grid [entry $wf.gnmEntry -width 20 -textvariable ::nmwiz::pyGNM] \
      -row 4 -column 2 -sticky ew
    grid [button $wf.gnmBrowse -text "Browse" \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"GNM Script" { gnm.py }}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::pyGNM $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 4 -column 3 -sticky ew

    grid [button $wf.pcaHelp -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Full path to the ProDy PCA script (pca.py),\
                    e.g. C:\\python27\\Scripts\\pca.py"}] \
      -row 5 -column 0 -sticky w
    grid [label $wf.pcaLabel -text "PCA script:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wf.pcaEntry -width 20 -textvariable ::nmwiz::pyPCA] \
      -row 5 -column 2 -sticky ew
    grid [button $wf.pcaBrowse -text "Browse" \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"PCA Script" { pca.py }}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::pyPCA $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 5 -column 3 -sticky ew

    grid [button $wf.prodySubmit -text "Save and Close" \
        -command "::nmwiz::saveSettings; destroy .nmwizsettings"] \
      -row 15 -column 0 -columnspan 4 -sticky we

    pack $wf -side top -fill x -expand 1
  }
  
  variable prodyMolecule
  variable prodyMolid -1
  variable prodyNFrames 0
  variable prodySelstr "protein name CA"
  variable prodySelAtoms 0
  variable prodyScript "ANM"
  variable prodyPrefix ""
  variable prodyRmCoords 0
  variable prodyPCAfile "DCD"
  variable prodyPCAfiletype "DCD file"
  variable prodyAllfig 0
  variable prodyAllnum 0
  variable prodyTask ""
  variable prodyFrame 0
  variable prodyCutoff 15
  variable prodyGNMCutoff 10
  variable prodyGamma 1
  variable prodyNModes 10
  variable prodyFirstFrame 0
  variable prodySkipFrame 0 
  variable prodyLastFrame end

  variable fromolMolecule
  variable fromolMolid -1
  variable fromolNFrames 0
  variable fromolSelstr "all"
  variable fromolSelAtoms 0
  variable fromolFrame 0
  variable fromolFirstFrame 1
  variable fromolLastFrame end
  
  variable strcompRefSelstr "name CA and protein"
  variable strcompTarSelstr "name CA and protein"
  variable strcompRefMol ""
  variable strcompTarMol ""
  variable strcompRefFrame 0
  variable strcompTarFrame 0
  variable strcompRefid -1
  variable strcompTarid -1
  variable strcompRefN 0
  variable strcompTarN 0
  
  proc initStrComp {} {
    variable strcompGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizstrcomp] {
      wm deiconify .nmwizstrcomp
      raise .nmwizstrcomp
      return 
    }    
    set strcompGUI [toplevel .nmwizstrcomp]
    wm title $strcompGUI "NMWiz - Structure Comparison"
    wm resizable $strcompGUI 0 0

    # Main frame (molecule and selection)
    set wmf [labelframe $strcompGUI.mainFrame -text "Molecule Selection" -bd 2]
    grid [label $wmf.molRLabel -text "Reference"] \
      -row 1 -column 2 -sticky w
    grid [label $wmf.molTLabel -text "Target"] \
      -row 1 -column 3 -sticky w

    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.refFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.refFrame.list ::nmwiz::strcompRefMol "" 
    grid [frame $wmf.tarFrame] \
      -row 2 -column 3 -sticky ew
    tk_optionMenu $wmf.tarFrame.list ::nmwiz::strcompTarMol "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::nmwiz::strcompUpdateMolList] \
      -row 2 -column 4 -sticky ew
      
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrRefEntry -width 12 -textvariable ::nmwiz::strcompRefSelstr] \
      -row 5 -column 2 -sticky ew
    grid [entry $wmf.selstrTarEntry -width 12 -textvariable ::nmwiz::strcompTarSelstr] \
      -row 5 -column 3 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::nmwiz::strcompUpdateMolinfo] \
      -row 5 -column 4 -sticky ew
      
    grid [label $wmf.frameLabel -text "Frame:"] \
      -row 7 -column 1 -sticky w
    grid [entry $wmf.frameRefEntry -width 4 -textvariable ::nmwiz::strcompRefFrame] \
      -row 7 -column 2 -sticky w
    grid [entry $wmf.frameTarEntry -width 4 -textvariable ::nmwiz::strcompTarFrame] \
      -row 7 -column 3 -sticky w

    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 8 -column 1 -sticky w
    grid [label $wmf.selinfoRefLabel -text ""] \
      -row 8 -column 2 -sticky w
    grid [label $wmf.selinfoTarLabel -text ""] \
      -row 8 -column 3 -sticky w
    grid [label $wmf.selinfoRMSDLabel -text ""] \
      -row 8 -column 4 -sticky w
      
    grid [button $wmf.showHelp -text "Help" \
        -command {::nmwiz::showHelp compare}] \
      -row 12 -column 1 -sticky we
    grid [button $wmf.rmsdUpdate -text "RMSD" \
        -command ::nmwiz::strcompRMSD] \
      -row 12 -column 2 -sticky ew
    grid [button $wmf.align -text "Align" \
        -command ::nmwiz::strcompAlign] \
      -row 12 -column 3 -sticky ew
    grid [button $wmf.prodySubmit -text "Calculate" \
        -command ::nmwiz::calcDeform] \
      -row 12 -column 4 -sticky we
      
    pack $wmf -side top -fill x -expand 1
    ::nmwiz::strcompUpdateMolList
    ::nmwiz::strcompUpdateMolinfo


  }

  proc initFromMolecule {} {
    variable fromolGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizfromol] {
      wm deiconify .nmwizfromol
      raise .nmwizfromol
      return 
    }    
    set fromolGUI [toplevel .nmwizfromol]
    wm title $fromolGUI "NMWiz - From Molecule"
    wm resizable $fromolGUI 0 0

    # Main frame (molecule and selection)
    set wmf [labelframe $fromolGUI.mainFrame -text "Molecule Selection" -bd 2]
    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.molFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.molFrame.list ::nmwiz::fromolMolecule "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::nmwiz::fromolUpdateMolList] \
      -row 2 -column 3 -sticky ew
      
    grid [label $wmf.molinfoLbl -text "Information:"] \
      -row 3 -column 1 -sticky w
    grid [label $wmf.molinfoLabel -text ""] \
      -row 3 -column 2 -columnspan 2 -sticky w
    
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrEntry -width 20 -textvariable ::nmwiz::fromolSelstr] \
      -row 5 -column 2 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::nmwiz::fromolUpdateSelection] \
      -row 5 -column 3 -sticky ew
      
    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 6 -column 1 -sticky w
    grid [label $wmf.selinfoLabel -text ""] \
      -row 6 -column 2 -columnspan 2 -sticky w

    grid [label $wmf.frameLabel -text "Coordinate frame:"] \
      -row 7 -column 1 -sticky w
    grid [entry $wmf.frameEntry -width 4 -textvariable ::nmwiz::fromolFrame] \
      -row 7 -column 2 -sticky w
    
    grid [label $wmf.firstLabel -text "First mode frame:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wmf.firstEntry -width 4 -textvariable ::nmwiz::fromolFirstFrame] \
      -row 8 -column 2 -sticky w
    
    grid [label $wmf.lastLabel -text "Last mode frame:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wmf.lastEntry -width 4 -textvariable ::nmwiz::fromolLastFrame] \
      -row 10 -column 2 -sticky w
      
    grid [button $wmf.showHelp -text "Help" \
        -command {::nmwiz::showHelp frommolecule}] \
      -row 12 -column 1 -sticky we
    grid [button $wmf.prodySubmit -text "Load data from molecule" \
        -command ::nmwiz::fromMolecule] \
      -row 12 -column 2 -columnspan 2 -sticky we
      
    pack $wmf -side top -fill x -expand 1
    ::nmwiz::fromolUpdateMolList
    ::nmwiz::fromolUpdateMolinfo

  }

  proc initProdyGUI {} {
    variable prodyGUI
    # If already initialized, just turn on
    if [winfo exists .nmwizprody] {
      wm deiconify .nmwizprody
      raise .nmwizprody
      return 
    }    
    set prodyGUI [toplevel .nmwizprody]
    wm title $prodyGUI "NMWiz - ProDy Interface"
    wm resizable $prodyGUI 0 0
    
    # Main frame (molecule and selection)
    set wmf [labelframe $prodyGUI.mainFrame -text "Atom Selection" -bd 2]
    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.molFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.molFrame.list ::nmwiz::prodyMolecule "" 
    grid [button $wmf.molUpdate -text "Update" \
        -command ::nmwiz::prodyUpdateMolList] \
      -row 2 -column 3 -sticky ew
    
    grid [label $wmf.molinfoLbl -text "Information:"] \
      -row 3 -column 1 -sticky w
    grid [label $wmf.molinfoLabel -text ""] \
      -row 3 -column 2 -columnspan 2 -sticky w
    
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrEntry -width 20 -textvariable ::nmwiz::prodySelstr] \
      -row 5 -column 2 -sticky ew
    grid [button $wmf.selUpdate -text "Select" \
        -command ::nmwiz::prodyUpdateSelection] \
      -row 5 -column 3 -sticky ew
      
    grid [label $wmf.selinfoLbl -text "Information:"] \
      -row 6 -column 1 -sticky w
    grid [label $wmf.selinfoLabel -text ""] \
      -row 6 -column 2 -columnspan 2 -sticky w
    
    pack $wmf -side top -fill x -expand 1
    ::nmwiz::prodyUpdateMolList
    ::nmwiz::prodyUpdateMolinfo
     
    # ProDy job frame
    set wf [labelframe $prodyGUI.jobFrame -text "ProDy Job Settings" -bd 2]
      
    grid [label $wf.scriptLabel -text "ProDy job:"] \
      -row 7 -column 1 -sticky w
    grid [frame $wf.scriptFrame] \
      -row 7 -column 2 -sticky ew
    tk_optionMenu $wf.scriptFrame.list ::nmwiz::prodyTask "ANM calculation" 
    $wf.scriptFrame.list.menu delete 0 last
    foreach script "ANM GNM PCA" {
      $wf.scriptFrame.list.menu add radiobutton -label "$script calculation" \
          -variable ::nmwiz::prodyTask \
          -command "set ::nmwiz::prodyScript $script; ::nmwiz::prodyChangeTask; ::nmwiz::prodyUpdatePrefix"
      incr counter  
    }
    pack $wf.scriptFrame.list -side left -anchor w -fill x
    variable prodyTask "ANM calculation"
    grid [button $wf.nmwizSettings -text "Settings" \
        -command "::nmwiz::initSettingsGUI"] \
      -row 7 -column 3 -sticky ew

    grid [label $wf.outdLabel -text "Output directory:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.outdEntry -width 16 -textvariable ::nmwiz::outputdir] \
      -row 8 -column 2 -sticky ew
    grid [button $wf.outdBrowse -text "Browse" \
        -command {
      set tempdir [tk_chooseDirectory -initialdir $::nmwiz::outputdir ]
        if {![string equal $tempdir ""]} {set ::nmwiz::outputdir $tempdir}
        ::nmwiz::saveSettings
        }] \
      -row 8 -column 3 -sticky ew

    grid [label $wf.filepLabel -text "Output filename:"] \
      -row 9 -column 1 -sticky w
    grid [entry $wf.filepEntry -width 20 -textvariable ::nmwiz::prodyPrefix] \
      -row 9 -column 2 -columnspan 2 -sticky we

    grid [checkbutton $wf.rmfileEntry -text "remove coordinate file upon job completion." \
        -variable ::nmwiz::prodyRmCoords] \
      -row 10 -column 1 -columnspan 3 -sticky w

    pack $wf -side top -fill x -expand 1
    

    # ANM frame
    set wf [labelframe $prodyGUI.anmFrame -text "ANM Settings" -bd 2]
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::nmwiz::prodyNModes] \
      -row 6 -column 2 -sticky w   
    pack $wf -side top -fill x -expand 1
    
    grid [label $wf.frameLabel -text "Frame number:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.frameEntry -width 4 -textvariable ::nmwiz::prodyFrame] \
      -row 6 -column 4 -sticky w

    grid [label $wf.cutoffLabel -text "Cutoff distance (A):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.cutoffEntry -width 4 -textvariable ::nmwiz::prodyCutoff] \
      -row 8 -column 2 -sticky w

    grid [label $wf.gammaLabel -text "Force constant:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.gammaEntry -width 4 -textvariable ::nmwiz::prodyGamma] \
      -row 8 -column 4 -sticky w    

    # GNM frame
    set wf [labelframe $prodyGUI.gnmFrame -text "GNM Settings" -bd 2]
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::nmwiz::prodyNModes] \
      -row 6 -column 2 -sticky w   
    
    grid [label $wf.frameLabel -text "Frame number:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.frameEntry -width 4 -textvariable ::nmwiz::prodyFrame] \
      -row 6 -column 4 -sticky w

    grid [label $wf.cutoffLabel -text "Cutoff distance (A):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.cutoffEntry -width 4 -textvariable ::nmwiz::prodyGNMCutoff] \
      -row 8 -column 2 -sticky w

    grid [label $wf.gammaLabel -text "Force constant:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.gammaEntry -width 4 -textvariable ::nmwiz::prodyGamma] \
      -row 8 -column 4 -sticky w    

    # PCA frame
    set wf [labelframe $prodyGUI.pcaFrame -text "PCA (EDA) Settings" -bd 2]
    
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::nmwiz::prodyNModes] \
      -row 6 -column 2 -sticky w
      
    grid [label $wf.skipLabel -text "Skip frame:"] \
      -row 6 -column 3 -sticky w
    grid [entry $wf.skipEntry -width 4 -textvariable ::nmwiz::prodySkipFrame] \
      -row 6 -column 4 -sticky w

    grid [label $wf.firstLabel -text "First frame:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.firstEntry -width 4 -textvariable ::nmwiz::prodyFirstFrame] \
      -row 8 -column 2 -sticky w    

    grid [label $wf.lastLabel -text "Last frame:"] \
      -row 8 -column 3 -sticky w
    grid [entry $wf.lastEntry -width 4 -textvariable ::nmwiz::prodyLastFrame] \
      -row 8 -column 4 -sticky w

    grid [label $wf.filetypeLabel -text "Trajectory file type:"] \
      -row 10 -column 1 -sticky w
    grid [frame $wf.filetypeFrame] \
      -row 10 -column 2 -columnspan 3 -sticky ew
    tk_optionMenu $wf.filetypeFrame.list ::nmwiz::prodyPCAfiletype "dcd" 
    $wf.filetypeFrame.list.menu delete 0 last
    foreach script "DCD PDB" {
      $wf.filetypeFrame.list.menu add radiobutton -label "$script file" \
          -variable ::nmwiz::prodyPCAfiletype \
          -command "set ::nmwiz::prodyPCAfile $script;"
      incr counter  
    }
    pack $wf.filetypeFrame.list -side left -anchor w -fill x
    variable prodyPCAfiletype "DCD file"



      
      
    # Submit button
    set wf [frame $prodyGUI.submitFrame -bd 2]
      
    grid [button $wf.showHelp -text "Help" \
        -command {::nmwiz::prodySubmitJob prody}] \
      -row 0 -column 0 -sticky we
    grid [button $wf.prodySubmit -text "Submit Job" \
        -command ::nmwiz::prodySubmitJob] \
      -row 0 -column 1 -sticky we
    grid [button $wf.prodyWebsite -text "ProDy Website" \
        -command "vmd_open_url http://www.csb.pitt.edu/ProDy"] \
      -row 0 -column 2 -sticky we
    pack $wf -side top -fill x -expand 1
         
  }
  
  proc prodyUpdateMolList {} {
    variable prodyGUI
    set wf $prodyGUI.mainFrame
    $wf.molFrame.list.menu delete 0 last
    set counter 0
    variable prodyMolid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 0} {
        if {$counter == 0} {
          variable prodyMolid $id
        }
        $wf.molFrame.list.menu add radiobutton -label "[::nmwiz::cleanMolName $id] ($id)" \
            -variable ::nmwiz::prodyMolecule \
            -command "set ::nmwiz::prodyMolid $id; ::nmwiz::prodyUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.molFrame.list -side left -anchor w -fill x
    variable prodyMolid
    if {$prodyMolid > -1} {
      variable prodyMolecule "[::nmwiz::cleanMolName $prodyMolid] ($prodyMolid)"
    }
    ::nmwiz::prodyUpdateMolinfo
  }

  proc strcompUpdateMolList {} {
    variable strcompGUI
    set wf $strcompGUI.mainFrame
    $wf.refFrame.list.menu delete 0 last
    $wf.tarFrame.list.menu delete 0 last
    set counter 0
    variable strcompRefid -1
    variable strcompTarid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 1} {
        if {$counter == 0} {
          variable strcompRefid $id
          variable strcompTarid $id
        }
        $wf.refFrame.list.menu add radiobutton -label "[::nmwiz::cleanMolName $id] ($id)" \
            -variable ::nmwiz::strcompRefMolecule \
            -command "set ::nmwiz::strcompRefid $id; ::nmwiz::strcompUpdateMolinfo"
        $wf.tarFrame.list.menu add radiobutton -label "[::nmwiz::cleanMolName $id] ($id)" \
            -variable ::nmwiz::strcompTarMolecule \
            -command "set ::nmwiz::strcompTarid $id; ::nmwiz::strcompUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.refFrame.list -side left -anchor w -fill x
    pack $wf.tarFrame.list -side left -anchor w -fill x
    if {$strcompTarid > -1} {
      variable strcompTarMol "[::nmwiz::cleanMolName $strcompTarid] ($strcompTarid)"
      variable strcompRefMol "[::nmwiz::cleanMolName $strcompRefid] ($strcompRefid)"
    }
    ::nmwiz::strcompUpdateMolinfo
  }
  
  proc strcompUpdateMolinfo {} {
    variable strcompRefid
    variable strcompTarid
    if {$strcompRefid > -1} {
      set ::nmwiz::strcompRefMol "[::nmwiz::cleanMolName $strcompRefid] ($strcompRefid)"
      set ref [atomselect $strcompRefid "$::nmwiz::strcompRefSelstr" frame $::nmwiz::strcompRefFrame]
      set ::nmwiz::strcompRefN [$ref num]
      .nmwizstrcomp.mainFrame.selinfoRefLabel configure \
        -text "$::nmwiz::strcompRefN selected"
      $ref delete
    } else {
      set ::nmwiz::strcompRefMol ""
      .nmwizstrcomp.mainFrame.selinfoRefLabel configure \
        -text "0 selected"
    }
    if {$strcompTarid > -1} {
      set ::nmwiz::strcompTarMol "[::nmwiz::cleanMolName $strcompTarid] ($strcompTarid)"
      set tar [atomselect $strcompTarid "$::nmwiz::strcompTarSelstr" frame $::nmwiz::strcompTarFrame]
      set ::nmwiz::strcompTarN [$tar num]
      .nmwizstrcomp.mainFrame.selinfoTarLabel configure \
        -text "$::nmwiz::strcompTarN selected"
      $tar delete
    } else {
      set ::nmwiz::strcompTarMol ""
      .nmwizstrcomp.mainFrame.selinfoTarLabel configure \
        -text "0 selected"
    }
    
  }
  proc strcompAlign {} {
    set ref [atomselect $::nmwiz::strcompRefid "$::nmwiz::strcompRefSelstr" frame $::nmwiz::strcompRefFrame]
    set tar [atomselect $::nmwiz::strcompTarid "$::nmwiz::strcompTarSelstr" frame $::nmwiz::strcompTarFrame]
    if {[$ref num] == [$tar num]} {
      set all [atomselect $::nmwiz::strcompTarid "all" frame $::nmwiz::strcompTarFrame]
      $all move [measure fit $tar $ref]
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "RMSD = [format %.2f [measure rmsd $ref $tar]] A"
      $all delete
    } else {      
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "Length mismatch"
    }
    $tar delete
    $ref delete
  }  
  proc strcompRMSD {} {
    set ref [atomselect $::nmwiz::strcompRefid "$::nmwiz::strcompRefSelstr" frame $::nmwiz::strcompRefFrame]
    set tar [atomselect $::nmwiz::strcompTarid "$::nmwiz::strcompTarSelstr" frame $::nmwiz::strcompTarFrame]
    if {[$ref num] == [$tar num]} {
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "RMSD = [format %.2f [measure rmsd $ref $tar]] A"
    } else {
      .nmwizstrcomp.mainFrame.selinfoRMSDLabel configure \
        -text "Length mismatch"
    }
    $tar delete
    $ref delete
  }


  proc fromolUpdateMolList {} {
    variable fromolGUI
    set wf $fromolGUI.mainFrame
    $wf.molFrame.list.menu delete 0 last
    set counter 0
    variable fromolMolid -1
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 1} {
        if {$counter == 0} {
          variable fromolMolid $id
        }
        $wf.molFrame.list.menu add radiobutton -label "[::nmwiz::cleanMolName $id] ($id)" \
            -variable ::nmwiz::fromolMolecule \
            -command "set ::nmwiz::fromolMolid $id; ::nmwiz::fromolUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.molFrame.list -side left -anchor w -fill x
    variable fromolMolid
    if {$fromolMolid > -1} {
      variable fromolMolecule "[::nmwiz::cleanMolName $fromolMolid] ($fromolMolid)"
    } 
    ::nmwiz::fromolUpdateMolinfo
  }

  
  proc prodyCheckMolecule {} {
    if {[lsearch [molinfo list] $::nmwiz::prodyMolid] > -1 && [molinfo $::nmwiz::prodyMolid get numframes] > 0} {
      return 1
    } 
    ::nmwiz::prodyUpdateMolList
    return 0
  }
  
  proc fromolCheckMolecule {} {
    if {[lsearch [molinfo list] $::nmwiz::fromolMolid] > -1 && [molinfo $::nmwiz::fromolMolid get numframes] > 1} {
      return 1
    }
    ::nmwiz::fromolUpdateMolList
    return 0
  }
  
  proc prodyUpdateMolinfo {} {
    variable prodyMolid
    if {$prodyMolid > -1} {
      variable prodyMolecule "[::nmwiz::cleanMolName $prodyMolid] ($prodyMolid)"
      set ::nmwiz::prodyNFrames [molinfo $::nmwiz::prodyMolid get numframes]
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "[molinfo $::nmwiz::prodyMolid get numatoms] atoms, $::nmwiz::prodyNFrames frames"
      ::nmwiz::prodyUpdateSelection
      ::nmwiz::prodyUpdatePrefix
    } else {
      set ::nmwiz::prodyNFrames 0
      variable fromolMolecule ""
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "Load a molecule and click Update."
      .nmwizprody.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
    }
  }
  
  proc fromolUpdateMolinfo {} {
    variable fromolMolid
    if {$fromolMolid > -1} {
      variable fromolMolecule "[::nmwiz::cleanMolName $fromolMolid] ($fromolMolid)"
      set ::nmwiz::fromolNFrames [molinfo $::nmwiz::fromolMolid get numframes]
      .nmwizfromol.mainFrame.molinfoLabel configure \
        -text "[molinfo $::nmwiz::fromolMolid get numatoms] atoms, $::nmwiz::fromolNFrames frames"
      ::nmwiz::fromolUpdateSelection
    } else {
      set ::nmwiz::fromolMolecule ""
      set ::nmwiz::fromolNFrames 0
      .nmwizfromol.mainFrame.molinfoLabel configure \
        -text "Load a molecule and click Update."
      .nmwizfromol.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
        
    }
  }
  
  proc prodyUpdatePrefix {} {
    if {[::nmwiz::prodyCheckMolecule]} {
      set prefix [::nmwiz::cleanMolName $::nmwiz::prodyMolid]
      
      if {[string range $prefix [expr [string length $prefix] - 4] end] == ".pdb"} {
        set prefix [string range $prefix 0 [expr [string length $prefix] - 5]]
      }
      if {$::nmwiz::prodyScript == "ANM"} {
        set ::nmwiz::prodyPrefix "$prefix\_anm"
      } elseif {$::nmwiz::prodyScript == "GNM"} {
        set ::nmwiz::prodyPrefix "$prefix\_gnm"
      } else {
        set ::nmwiz::prodyPrefix "$prefix\_pca"
      }
    }
  }
  
  proc cleanMolName {molid} {
    set name "[molinfo $molid get name]"
    foreach char {"\}" "\{"} { 
      set first [string first $char $name] 
      while {$first > -1} {
        set name [string replace $name $first $first ""] 
        set first [string first $char $name]
      }
    }
    foreach char {" " "."} { 
      set first [string first $char $name] 
      while {$first > -1} {
        set name [string replace $name $first $first "_"] 
        set first [string first $char $name]
      }
    }
    return $name
  }
  
  proc prodyUpdateSelection {} {
    ::nmwiz::prodyCheckMolecule
    variable prodyMolid
    variable prodySelstr
    set sel [atomselect $prodyMolid $prodySelstr]
    variable prodySelAtoms [$sel num]
    $sel delete
    variable prodyGUI
    $prodyGUI.mainFrame.selinfoLabel configure \
      -text "$prodySelAtoms atoms are selected"
  }
  
  proc fromolUpdateSelection {} {
    ::nmwiz::fromolCheckMolecule
    variable fromolMolid
    variable fromolGUI
    if {$fromolMolid > -1} {
      
      variable fromolSelstr
      set sel [atomselect $fromolMolid $fromolSelstr]
      variable fromolSelAtoms [$sel num]
      $sel delete
      $fromolGUI.mainFrame.selinfoLabel configure \
        -text "$fromolSelAtoms atoms are selected"
    } else {
      $fromolGUI.mainFrame.selinfoLabel configure \
        -text "Load a molecule and click Update."
    }
  }
  
  proc prodyChangeTask {} {
    variable prodyGUI
    variable prodyScript
    if {$prodyScript == "ANM"} {
      pack forget $prodyGUI.gnmFrame
      pack forget $prodyGUI.pcaFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.anmFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    } elseif {$prodyScript == "GNM"} {
      pack forget $prodyGUI.pcaFrame
      pack forget $prodyGUI.anmFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.gnmFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    } else {
      pack forget $prodyGUI.anmFrame
      pack forget $prodyGUI.gnmFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.pcaFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    }
  }
  
  proc prodySubmitJob {} {
    if {$::nmwiz::prodySelAtoms == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "You need to make an atom selection before you can submit a job."
      return 
    }
    if {!([string is digit $::nmwiz::prodyNModes] && $::nmwiz::prodyNModes > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Number of modes must be a number larger than 0."
      return 
    }
    if {![file isdirectory $::nmwiz::outputdir]} {
      tk_messageBox -type ok -title "ERROR" \
        -message "$::nmwiz::outputdir is not a valid directory."
      return 
    }
    if {$::nmwiz::pybin == "" || $::nmwiz::pybin == {}} {
      variable ::nmwiz::pybin "[::ExecTool::find -interactive python.exe]"
      ::nmwiz::saveSettings
    }
    if {$::nmwiz::prodyScript == "ANM"} {
      ::nmwiz::prodySubmitANMjob
    } elseif {$::nmwiz::prodyScript == "GNM"} {
      ::nmwiz::prodySubmitGNMjob
    } else {
      ::nmwiz::prodySubmitPCAjob
    }
  }
  proc prodySubmitANMjob {} {
    if {$::nmwiz::pyANM == "" || $::nmwiz::pyANM == {} || $::nmwiz::pyANM == "{}" ||
        ![file exists $::nmwiz::pyANM]} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Please specify the path to the ProDy ANM Script (anm.py) script."
      ::nmwiz::initSettingsGUI
      return
    }
    set n_frames [molinfo $::nmwiz::prodyMolid get numframes]
    if {!([string is digit $::nmwiz::prodyFrame] && $::nmwiz::prodyFrame >= 0 && 
        $::nmwiz::prodyFrame < $n_frames)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame number must be an integer from 0 to [expr $n_frames - 1]."
      return 
    }
    if {!([string is double $::nmwiz::prodyCutoff] && $::nmwiz::prodyCutoff > 4.5)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance (A) must be a number greater than 4.5."
      return 
    }
    if {!([string is double $::nmwiz::prodyGamma] && $::nmwiz::prodyGamma > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Force constant must be a positive number."
      return 
    }
    set pdbfn [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix.pdb]
    set sel [atomselect $::nmwiz::prodyMolid $::nmwiz::prodySelstr]
    $sel frame $::nmwiz::prodyFrame
    $sel writepdb $pdbfn
    $sel delete
    
    set allnum ""
    if {$::nmwiz::prodyAllnum} {
      set allnum "-a"
    } 
    set allfig ""
    if {$::nmwiz::prodyAllfig} {
      set allfig "-A"
    }    
    set prefix [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix]    
    vmdcon -info "Executing: $::nmwiz::pybin $::nmwiz::pyANM --quiet -s \"all\" -o \"$::nmwiz::outputdir\" -p \"$prefix\" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyCutoff -g $::nmwiz::prodyGamma \"$pdbfn\""
    set status [exec $::nmwiz::pybin $::nmwiz::pyANM --quiet -s all -o "$::nmwiz::outputdir" -p "$prefix" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyCutoff -g $::nmwiz::prodyGamma "$pdbfn"]

    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy ANM calculation is finished and results are being loaded."
      ::nmwiz::loadNMD "$prefix.nmd"
      if {$::nmwiz::prodyRmCoords} {
        file delete -force $pdbfn
      }
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
    }
  }  
  proc prodySubmitGNMjob {} {
    if {$::nmwiz::pyGNM == "" || $::nmwiz::pyGNM == {} || $::nmwiz::pyGNM == "{}" ||
        ![file exists $::nmwiz::pyGNM]} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Please specify the path to the ProDy GNM Script (gnm.py) script."
      ::nmwiz::initSettingsGUI
      return
    }
    set n_frames [molinfo $::nmwiz::prodyMolid get numframes]
    if {!([string is digit $::nmwiz::prodyFrame] && $::nmwiz::prodyFrame >= 0 && 
        $::nmwiz::prodyFrame < $n_frames)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame number must be an integer from 0 to [expr $n_frames - 1]."
      return 
    }
    if {!([string is double $::nmwiz::prodyGNMCutoff] && $::nmwiz::prodyGNMCutoff > 4.5)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance (A) must be a number greater than 4.5."
      return 
    }
    if {!([string is double $::nmwiz::prodyGamma] && $::nmwiz::prodyGamma > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Force constant must be a positive number."
      return 
    }
    set pdbfn [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix.pdb]
    set sel [atomselect $::nmwiz::prodyMolid $::nmwiz::prodySelstr]
    $sel frame $::nmwiz::prodyFrame
    $sel writepdb $pdbfn
    $sel delete
    
    set allnum ""
    if {$::nmwiz::prodyAllnum} {
      set allnum "-a"
    } 
    set allfig ""
    if {$::nmwiz::prodyAllfig} {
      set allfig "-A"
    }    
    set prefix [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix]    
    vmdcon -info "Executing: $::nmwiz::pybin $::nmwiz::pyGNM --quiet -s \"all\" -o \"$::nmwiz::outputdir\" -p \"$prefix\" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyGNMCutoff -g $::nmwiz::prodyGamma \"$pdbfn\""
    set status [exec $::nmwiz::pybin $::nmwiz::pyGNM --quiet -s all -o "$::nmwiz::outputdir" -p "$prefix" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyGNMCutoff -g $::nmwiz::prodyGamma "$pdbfn"]

    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy GNM calculation is finished and results are being loaded."
      ::nmwiz::loadNMD "$prefix.nmd"
      if {$::nmwiz::prodyRmCoords} {
        file delete -force $pdbfn
      } 
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
    }
  }  
  proc prodySubmitPCAjob {} {
    if {$::nmwiz::pyPCA == "" || $::nmwiz::pyPCA == {} || $::nmwiz::pyPCA == "{}" ||
        ![file exists $::nmwiz::pyPCA]} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Please specify the path to the ProDy PCA Script (pca.py) script."
      ::nmwiz::initSettingsGUI
      return
    }
    if {$::nmwiz::prodyNFrames < 2} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Selected molecule must have more than 1 frames for PCA calculations."
      return       
    }
    if {!([string is digit $::nmwiz::prodyFirstFrame] && $::nmwiz::prodyFirstFrame >= 0 && 
        $::nmwiz::prodyFirstFrame < [molinfo $::nmwiz::prodyMolid get numframes])} {
      tk_messageBox -type ok -title "ERROR" \
        -message "First frame must be a number and must be in the valid range."
      return 
    }
    if {!([string is digit $::nmwiz::prodySkipFrame] && $::nmwiz::prodySkipFrame >= 0 && 
        $::nmwiz::prodySkipFrame < [molinfo $::nmwiz::prodyMolid get numframes])} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Skip frame must be a number and must be in the valid range."
      return 
    }
    if {!($::nmwiz::prodyLastFrame == "end" || ([string is digit $::nmwiz::prodyLastFrame]
       && $::nmwiz::prodyLastFrame > 0 && $::nmwiz::prodyLastFrame < [molinfo $::nmwiz::prodyMolid get numframes]))} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Last frame may be \"end\" or a number in the valid range."
      return 
    }
    set sel [atomselect $::nmwiz::prodyMolid $::nmwiz::prodySelstr]
    set pdbfn [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix.[string tolower $::nmwiz::prodyPCAfile]]
    set end $::nmwiz::prodyLastFrame 
    if {$end == "end"} {
      set end [expr $::nmwiz::prodyNFrames -1]
    }
    #vmdcon -info "animate write pdb $pdbfn beg $::nmwiz::prodyFirstFrame end $end skip $::nmwiz::prodySkipFrame waitfor all sel $sel $::nmwiz::prodyMolid"
    if {$::nmwiz::prodyPCAfile == "DCD"} {
      set nwritten [animate write dcd $pdbfn beg $::nmwiz::prodyFirstFrame end $end skip $::nmwiz::prodySkipFrame waitfor all sel $sel $::nmwiz::prodyMolid]
    } else {
      set nwritten [animate write pdb $pdbfn beg $::nmwiz::prodyFirstFrame end $end skip $::nmwiz::prodySkipFrame waitfor all sel $sel $::nmwiz::prodyMolid]
    }
    vmdcon -info "$nwritten frames are written as $pdbfn"
    
    $sel delete
    set allnum ""
    if {$::nmwiz::prodyAllnum} {
      set allnum "-a"
    } 
    set allfig ""
    if {$::nmwiz::prodyAllfig} {
      set allfig "-A"
    } 
    set prefix [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix]
    vmdcon -info "Executing: $::nmwiz::pybin $::nmwiz::pyPCA --quiet -s all -o \"$::nmwiz::outputdir\" -p \"$prefix\" -n $::nmwiz::prodyNModes \"$pdbfn\""
    set status [exec $::nmwiz::pybin $::nmwiz::pyPCA --quiet -s all -o "$::nmwiz::outputdir" -p "$prefix" -n $::nmwiz::prodyNModes "$pdbfn"]
    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy PCA calculation is finished and results are being loaded."
      ::nmwiz::loadNMD "$prefix.nmd"
      if {$::nmwiz::prodyRmCoords} {
        file delete -force $pdbfn
      }  
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
      file delete -force $pdbfn
    }
  }
  proc init_anm_interface {} {
    variable a
    # If already initialized, just turn on
    if [winfo exists .nmwizanm] {
      wm deiconify .nmwizanm
      raise .nmwizanm
      return 
    }
    set a [toplevel .nmwizanm]
    wm title $a "NMWiz - ANM Server Interface"
    wm resizable $a 0 0


    variable anm_id
    variable anm_coormod "pdb"
    variable anm_chain "*"
    variable anm_model "all"
    variable anm_cutoff 15
    variable anm_pwr 0
    
    set wmf [frame $a.mainframe -bd 2]
    
    grid [button $wmf.id_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter 4 character PDB identifier."}] \
      -row 2 -column 0 -sticky w
    grid [label $wmf.id_label -text "PDB identifier:"] \
      -row 2 -column 1 -sticky w
    grid [entry $wmf.id_entry -width 4 -textvariable ::nmwiz::anm_id] \
      -row 2 -column 2 -sticky w
      
    grid [button $wmf.coordmod_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Choose whether to use PDB coordinates or biological unit coordinates."}] \
      -row 3 -column 0 -sticky w
    grid [radiobutton $wmf.coordmod_pdb -text "pdb coordinates" \
            -variable ::nmwiz::anm_coormod -value "pdb"] \
      -row 3 -column 1 -sticky w
    grid [radiobutton $wmf.coordmod_bio -text "biological unit" \
            -variable ::nmwiz::anm_coormod -value "bio"] \
      -row 3 -column 2 -sticky w
    
    grid [button $wmf.chain_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter 1 character chain identifier (default: all polypeptide chains)."}] \
      -row 4 -column 0 -sticky w
    grid [label $wmf.chain_label -text "Chain identifier:"] \
      -row 4 -column 1 -sticky w
    grid [entry $wmf.chain_entry -width 4 -textvariable ::nmwiz::anm_chain] \
      -row 4 -column 2 -sticky w

    grid [button $wmf.model_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter model (for multi-model files such as NMR structures)."}] \
      -row 6 -column 0 -sticky w
    grid [label $wmf.model_label -text "Model number:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wmf.model_entry -width 4 -textvariable ::nmwiz::anm_model] \
      -row 6 -column 2 -sticky w

    grid [button $wmf.cutoff_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter cutoff for interaction between C alpha atoms (Å)."}] \
      -row 8 -column 0 -sticky w
    grid [label $wmf.cutoff_label -text "Cutoff distance (Å):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wmf.cutoff_entry -width 4 -textvariable ::nmwiz::anm_cutoff] \
      -row 8 -column 2 -sticky w

    grid [button $wmf.pwr_help -text "?" \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter distance weight for interaction between C alpha atoms."}] \
      -row 10 -column 0 -sticky w
    grid [label $wmf.pwr_label -text "Distance weight:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wmf.pwr_entry -width 4 -textvariable ::nmwiz::anm_pwr] \
      -row 10 -column 2 -sticky w
      
    grid [button $wmf.retrieve -width 24 -text "Submit to ANM Server" \
      -command ::nmwiz::submit_anm_server] \
      -row 16 -column 0 -columnspan 3
      
    grid [button $wmf.website -width 24 -text "Go to ANM Server" \
        -command "vmd_open_url http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm1.cgi"] \
      -row 18 -column 0 -columnspan 3

    pack $wmf -side top -fill x -expand 1
  }
  
  proc submit_anm_server {} {
    variable anm_id
    if {[string length $anm_id] != 4} {
      tk_messageBox -type ok -title "ERROR" \
        -message "PDB identifier must be 4 characters."
      return
    }
    variable anm_coormod
    variable anm_chain
    if {! [string is alpha $anm_chain] && $anm_chain != "*"} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Chain identifier must be a letter or *."
      return
    }
    variable anm_model
    if {![string is integer $anm_model] && $anm_model != "all"} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Model identifier must be an integer or \"all\"."
      return
    }
    variable anm_cutoff
    if {![string is double $anm_cutoff] || $anm_cutoff < 6.0 || $anm_cutoff > 30.0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance must be between 6.0 and 30.0 A."
      return
    }
    variable anm_pwr
    if {![string is double $anm_pwr] || $anm_pwr < 0.0 || $anm_pwr > 3.0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Distance weight must be between 0.0 and 3.0 A."
      return
    }
    
    vmd_open_url http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm_blzpack.cgi?id=$anm_id&coormod=$anm_coormod&chain=$anm_chain&model=$anm_model&cutoff=$anm_cutoff&pwr=$anm_pwr
    
  
  }
  
  proc writeNMD {ns} {
    
    set tempfile [tk_getSaveFile -filetypes {{"NMD files" { .nmd .NMD }} {"All files" *}}]
    if {$tempfile == ""} {return}
    
    set fl [open $tempfile w]
    puts $fl "nmwiz_load $tempfile"
    puts $fl "name [subst $${ns}::title]"
    puts $fl "atomnames [subst $${ns}::atomnames]"
    puts $fl "resnames [subst $${ns}::resnames]"
    puts $fl "resids [subst $${ns}::resids]"
    puts $fl "chainids [subst $${ns}::chainids]"
    puts $fl "bfactors [subst $${ns}::bfactors]"
    puts $fl "coordinates [subst $${ns}::coordinates]"
    set indices [subst $${ns}::indices] 
    set lengths [subst $${ns}::lengths] 
    set modes [subst $${ns}::modes]
    for {set i 0} {$i < [llength $indices]} {incr i} {
      puts $fl "mode [lindex $indices $i] [lindex $lengths $i] [lindex $modes $i]"
    }
    close $fl    
  }
  
  proc loadNMD {fn} {
    if {![file isfile $fn]} {
      vmdcon -err "$fn is not a valid NMD file path."
      return
    }
    variable filename $fn
    vmdcon -info "NMWiz: Parsing file $filename"
    # Parse the file, and make sure coordinates are stored in the file
    #variable namespaces
    #variable nmwizguis
    
    set nmdfile [open $filename]
    set nmdlist [list]
    set coordinates 0 
    set n_dims 0
    while {[gets $nmdfile nmdline] != -1} { 
      if {[lindex $nmdline 0] == "coordinates"} {
        if {[expr [llength $nmdline] % 3] != 1} {
          tk_messageBox -type ok -title "ERROR" \
            -message "Length of the coordinate array in $filename must be a\
                      multiple of 3. An array of length\
                      [llength $coordinates] is provided."
          return
        }
        set coordinates [lrange $nmdline 1 end]
        set n_atoms [expr [llength $coordinates] / 3]
      }
      lappend nmdlist $nmdline
    }
    if {$coordinates == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Coordinate data was not found in the input file.\
                  NMD files must contain system coordinate data."
      return
    }
    
    # Evaluate each line
    set title ""
    set atomnames ""
    set resnames ""
    set resids ""
    set bfactors ""
    set chainids ""
    set modes [list]
    set modelength 0
    set modecounter 0
    set dof 0
    foreach nmdline $nmdlist {
      switch -exact [lindex $nmdline 0] {
        name {
          set title [lrange $nmdline 1 end]
        }
        coordinates {
        }
        atomnames {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of atomnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set atomnames [lrange $nmdline 1 end]
          }
        }
        resnames {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of resnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set resnames [lrange $nmdline 1 end]
          }
        }
        chainids {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of chainids array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set chainids [lrange $nmdline 1 end]
          }
        }
        resids {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of resids array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set resids [lrange $nmdline 1 end]
          }
        }
        bfactors {
          if {[llength $nmdline] != [expr $n_atoms + 1]} {
            vmdcon -info "NMWiz WARNING: Length of bfactors array must be $n_atoms, not [expr [llength $nmdline] -1]."
          } else {
            set bfactors [lrange $nmdline 1 end]
          }      
        }
        mode {
          incr modecounter
          set l [expr [llength $nmdline] - 1]
          if {$n_dims == 0} {
            set modelength $l
            if {$l >= [llength $coordinates]} {
              set diff [expr $l - [llength $coordinates]] 
              if {$diff > 2} {
                tk_messageBox -type ok -title "ERROR" \
                  -message "First mode line is not formatted correctly."
                return
              }
              set dof [llength $coordinates]
              set n_dims 3
              vmdcon -info "NMWiz INFO: File contains a 3D model."
            } else {
              set diff [expr $l - $n_atoms]
              if {$diff > 2} {
                tk_messageBox -type ok -title "ERROR" \
                  -message "First mode line is not formatted correctly."
                return
              }
              set dof $n_atoms
              set n_dims 1
              vmdcon -info "NMWiz INFO: File contains a 1D model."
            }
          }
          if {$modelength > 0 && $l != $modelength} {
            vmdcon -info "NMWiz WARNING: Mode line $modecounter does not have the same length as the first mode line."
          } else {
            switch -exact [expr $l - $dof] {
              0 {
                lappend modes [lrange $nmdline 1 end]
                lappend indices [llength $modes]
                lappend lengths 1
              }
              1 {
                lappend modes [lrange $nmdline 2 end]
                if {[string is integer [lindex $nmdline 1]]} {
                  lappend indices [lindex $nmdline 1]
                  lappend lengths 1
                } else {
                  lappend lengths [lindex $nmdline 1]
                  lappend indices [llength $modes]
                }
              }
              2 {
                lappend modes [lrange $nmdline 3 end]
                if {[string is integer [lindex $nmdline 1]]} {
                  lappend indices [lindex $nmdline 1]
                  lappend lengths [lindex $nmdline 2]
                } else {
                  lappend indices [lindex $nmdline 2]
                  lappend lengths [lindex $nmdline 1]
                }
              } 
              default {
                vmdcon -info "NMWiz WARNING: Mode data was not understood. Line starts with [lrange $nmdline 0 4]."
              }
            } 
          }
        }
        default {
          vmdcon -info "NMWiz WARNING: Unrecognized line starting with \"[lindex $nmdline 0]\""
        }
      }
      
    } 
    
    if {[llength $modes] == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Mode data was not found in the input file."
      return
    }

    if {$title == ""} {
      set title "Untitled ($::nmwiz::guicount)"
      vmdcon -info "NMWiz INFO: Dataset is named as \"$title\"."
    }

    set ns [::nmwiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes $n_dims $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::nmwiz::appendGUIcontrols $ns
    
  }

  proc calcDeform {} {
    variable strcompRefN
    variable strcompTarN
    variable strcompRefid
    variable strcompTarid
    variable strcompRefFrame 
    variable strcompTarFrame
    variable strcompRefMol 
    variable strcompTarMol
    variable strcompRefSelstr 
    variable strcompTarSelstr
    if {$strcompRefid == $strcompTarid && $strcompRefFrame == $strcompTarFrame} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Reference and target cannot be the same frame of the same molecule."
      return
    }
    if {$strcompRefN == 0 || $strcompRefN != $strcompTarN} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Reference and target selections must have same number of atoms."
      return
    }
    set title "Deformation $strcompRefMol -> $strcompTarMol"
    set ref [atomselect $strcompRefid $strcompRefSelstr frame $strcompRefFrame]
    set coordinates [concat {*}[$ref get {x y z}]]
    set modes [list]
    set lengths [list]
    set indices [list]
    set tar [atomselect $strcompTarid $strcompTarSelstr frame $strcompTarFrame]
    set target [concat {*}[$tar get {x y z}]]
    $tar delete
    set m [vecsub $target $coordinates]
    set len [veclength $m]
    lappend indices 1
    lappend modes [vecscale $m [expr 1 / $len]]
    lappend lengths $len
    set atomnames [$ref get name]
    set resnames [$ref get resname]
    set resids [$ref get resid]
    set chainids [$ref get chain]
    set bfactors [$ref get beta]
    $ref delete
    
    set ns [::nmwiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes 3 $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::nmwiz::appendGUIcontrols $ns
    
  }

  proc fromMolecule {} {  
    variable fromolMolecule
    variable fromolMolid
    variable fromolSelstr
    variable fromolFrame
    variable fromolFirstFrame
    variable fromolLastFrame
    set lastframe $fromolLastFrame
    if {$lastframe == "end"} {
      set lastframe [expr [molinfo $fromolMolid get numframes] - 1]
    }
    if {$fromolFrame >= $fromolFirstFrame && $fromolFrame <= $lastframe} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Coordinate frame cannot be one of the mode frames."
      return
    }
    set title $fromolMolecule
    set n_dims 3
    set sel [atomselect $fromolMolid $fromolSelstr frame $fromolFrame]
    set coordinates [concat {*}[$sel get {x y z}]]
    set modes [list]
    set lengths [list]
    set indices [list]
    for {set i $fromolFirstFrame} {$i <= $lastframe} {incr i} {
      $sel frame $i
      set m [concat {*}[$sel get {x y z}]]
      set length [veclength $m]
      lappend lengths $length  
      lappend modes [vecscale $m [expr 1 / $length]]
      lappend indices $i
    }
    $sel frame $fromolFrame
    set atomnames [$sel get name]
    set resnames [$sel get resname]
    set resids [$sel get resid]
    set chainids [$sel get chain]
    set bfactors [$sel get beta]
    
    $sel delete
    
    set ns [::nmwiz::makeNMWizGUI]
    ${ns}::initialize $coordinates $modes $n_dims $title $lengths $indices $atomnames $resnames $resids $chainids $bfactors
    ${ns}::nmwizgui
    ::nmwiz::appendGUIcontrols $ns
    
  }
  
  proc appendGUIcontrols {ns} {
    
    set w .nmwizgui
    set wgf [labelframe $w.{[string range $ns 2 end]}frame -text "[subst $${ns}::title]" -bd 2]
    
    grid [button $wgf.show -text "Show GUI" \
        -command "${ns}::nmwizgui" ] \
      -row 0 -column 0 -sticky we
    grid [button $wgf.remove -text "Remove" \
        -command "lset ::nmwiz::titles $::nmwiz::guicount NONE; pack forget $wgf; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
      -row 0 -column 1 -sticky we
    grid [button $wgf.save -text "Save" \
        -command "::nmwiz::writeNMD $ns"] \
      -row 0 -column 2 -sticky we

    pack $wgf -side top -fill x -expand 1
    
  }
  
  proc array2xyz {arr} {
    set coords [list]
    foreach {x y z} $arr {
      lappend coords "$x $y $z"
    }
    return $coords
  }
  
  proc makeNMWizGUI {} {

    variable guicount
    incr guicount
    set ns "::nmgui$guicount"
    
    namespace eval $ns {
      variable tempfn ".[string range [namespace current] 2 end].pdb"
      variable w
      variable name
      variable molid -1
      variable arrid -1
      variable arridlist [list]
      variable animid -1
      variable animidlist [list]
      variable stopped 1 
      variable arrids
      variable scalearrows 1
      variable sense +1
      variable color $::nmwiz::defaultColor
      variable colorlist [list]
      variable materials on
      variable material "HardPlastic"
      variable material_protein "HardPlastic"
      variable resolution 10
      variable resolution_protein 10
      variable ndim 3
      variable bothdirections 0
      
      #GNM option
      variable msformode "Mobility"
      
      variable selstr "all"
      variable selrep 0
      
      variable autoupdate 1
      variable autoanimate 0
      
      variable scalearrows_list [list]
      variable hide_shorter_list [list] 
      variable cylinder_radius_list [list]
      variable cone_radius_list [list]
      variable cone_height_list [list]
      variable resolution_list [list]
      variable material_list [list]

      variable hide_shorter 0.0 
      variable cylinder_radius 0.3
      variable cone_radius 0.5
      variable cone_height 1.0
      
      variable showproteinas "Tube"
      variable tuberadius 0.4
      variable bondradius 0.3   
      variable spherescale 0.6
      variable cutoffdistance 8.0
      variable proteincolor "Mobility" 
      
      variable nframes 50
      
      variable betalist {}
      variable betamin
      variable betamax
      
      variable overwrite 1
      variable autoplay 1

      variable arropt 0
      variable anmopt 0
      variable prtopt 0
      variable pltopt 0

      variable selid -1
      variable overplot 0
      variable selectscale 0.8
      variable plotwidth 800
      variable plotheight 600
      variable linewidth 1
      variable mradius 2
      variable dash "-"
      variable lornol "lines"
      variable marker "circle"
      variable plothandles [list]
      
      proc initialize {xyz m d t l i an rn ri ci bf} {
        variable arridlist
        variable animidlist
        variable colorlist
        variable coordinates $xyz
        variable modes $m
        variable ndim $d
        variable title $t
        variable lengths $l
        variable indices $i
        variable atomnames $an
        variable resnames $rn
        variable chainids $ci
        variable resids $ri
        variable bfactors $bf
        variable plotrids [list]
        variable bfactormin
        variable bfactormax
        variable n_dims 0
        variable n_atoms [expr [llength $coordinates] / 3]
        
        if {[lsearch $::nmwiz::titles $title] > -1} {
          set title "$title ($::nmwiz::guicount)"
        }
        lappend ::nmwiz::titles $title

        
        if {$atomnames == ""} {
          set atomnames [string repeat "CA " $n_atoms]
          vmdcon -info "NMWiz INFO: All atom names are set as \"CA\"."
        } else {
          variable showproteinas
          foreach an $atomnames {
            if {$an != "CA"} {
              set showproteinas "Licorice"
              break
            }
          }
        }
        if {$resnames == ""} {
          set resnames [string repeat "GLY " $n_atoms]
          vmdcon -info "NMWiz INFO: All residue names are named set as \"GLY\"."
        }
        if {$chainids == ""} {
          set chainids [string repeat "A " $n_atoms]
          vmdcon -info "NMWiz INFO: All chain identifiers are set as \"A\"."
        }
        if {[llength $resids] == 0} {
          for {set i 1} {$i <= $n_atoms} {incr i} {lappend resids $i}
          vmdcon -info "NMWiz INFO: Residues are numbered starting from 1."
        }

        foreach i [lrange $resids 0 end-1] j [lrange $resids 1 end] {
          if {$i >= $j} {
            for {set i 1} {$i <= $n_atoms} {incr i} {lappend plotrids $i}
            break
          }
        }
        if {[llength $plotrids] == 0} {
          set plotrids $ri
          vmdcon -info "NMWiz INFO: Residue numbers will be used for plotting."
        }

        if {[llength $bfactors] == 0} {
          vmdcon -info "NMWiz INFO: Experimental bfactors were not found in the data file."
          for {set i 1} {$i <= $n_atoms} {incr i} {lappend bfactors 0.0}
          set bfactormin 0.0
          set bfactormax 0.0
        } else {
          set bfactormin [lindex $bfactors 0]
          set bfactormax [lindex $bfactors 0]
          foreach i $bfactors {
            if {$i < $bfactormin} {set bfactormin $i}
            if {$i > $bfactormax} {set bfactormax $i}
          }
        }
        
        variable prefix $title
        while {[string first " " $prefix] > -1} {
          set i [string first " " $prefix]
          set prefix [string replace $prefix $i $i "_"]
        }
        while {[string first "." $prefix] > -1} {
          set i [string first "." $prefix]
          set prefix [string replace $prefix $i $i "_"]
        }

        variable numofmodes [llength $modes] 
        variable activeindex 0
        variable drawlength [lindex $lengths $activeindex]
        variable drawlengthstr [format "%.1f" [lindex $lengths $activeindex]]
        variable activemode [lindex $indices $activeindex]
        variable betalist
        foreach atnm $atomnames {
          lappend betalist 0
        } 
        variable scalearrows_list
        variable rmsd 2
        variable rmsdprev 0
        variable activeindexprev 0
        variable rmsd_list
        set one_over_root_of_n_atoms [expr 1 / $n_atoms ** 0.5]
        foreach len $lengths mode $modes {
            set l [::tcl::mathfunc::abs $len] 
            lappend scalearrows_list [expr $rmsd / $one_over_root_of_n_atoms / [veclength $mode] / $l]
            lappend rmsd_list $rmsd  
        }
        variable scalearrows
        set scalearrows [lindex $scalearrows_list 0]
        
        variable arridlist
        variable animidlist
        variable colorlist
        
        variable hide_shorter_list 
        variable cylinder_radius_list
        variable cone_radius_list
        variable cone_height_list
        variable resolution_list
        variable material_list
        
        variable hide_shorter 
        variable cylinder_radius
        variable cone_radius
        variable cone_height
        variable resolution
        variable material
        set shift 0
        for {set i 0} {$i < [llength $modes]} {incr i} {
          lappend arridlist -1
          lappend animidlist -1
          set curcolor [lindex $::nmwiz::nmwizColors [expr ($i + $shift + $::nmwiz::guicount) % [llength $::nmwiz::nmwizColors]]]
          if {$curcolor == $::nmwiz::defaultColor} {
            incr shift 
            set curcolor [lindex $::nmwiz::nmwizColors [expr ($i + $shift + $::nmwiz::guicount) % [llength $::nmwiz::nmwizColors]]]
          }
          lappend colorlist $curcolor
          lappend hide_shorter_list $hide_shorter 
          lappend cylinder_radius_list $cylinder_radius
          lappend cone_radius_list $cone_radius
          lappend cone_height_list $cone_height
          lappend resolution_list $resolution
          lappend material_list $material
        }
        if {$::nmwiz::guicount == 0} {
          lset colorlist 0 $::nmwiz::defaultColor  
        }
        variable color [lindex $colorlist 0]
        
      }
      
      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      proc deleteMolecules {} {
        variable molid
        variable animidlist
        variable selid
        variable arridlist
        if {[lsearch [molinfo list] $molid] > -1} {
          mol delete $molid
        }
        foreach arrid $arridlist {
          if {[lsearch [molinfo list] $arrid] > -1} {
            mol delete $arrid
          }
        }
        foreach animid $animidlist {
          if {[lsearch [molinfo list] $animid] > -1} {
            mol delete $animid
          }
        }
        if {[lsearch [molinfo list] $selid] > -1} {
          mol delete $selid
        }
      }
      
      proc Plot {ns} {
        
        #puts [subst $${ns}::overplot]
        set plothandle 0
        if {[subst $${ns}::overplot]} {
          for {set i [expr [llength [subst $${ns}::plothandles]] - 1]} {$i >= 0} {incr i -1} {
            set temphandle [lindex [subst $${ns}::plothandles] $i]
            if {[namespace exists [string range $temphandle 0 end-12]]} {
              set plothandle $temphandle
              break
            } else {
              set ${ns}::plothandles [lrange [subst $${ns}::plothandles] 0 $i]
            }
          }
          
        }
        if {$plothandle != 0} {
          $plothandle add \
            [subst $${ns}::plotrids] [subst $${ns}::betalist] \
            -title "[subst $${ns}::title] square fluctuations" \
            -linewidth [subst $${ns}::linewidth] \
            -legend "Mode [subst $${ns}::activemode]" \
            -[subst $${ns}::lornol] -linecolor [subst $${ns}::color] \
            -xsize [subst $${ns}::plotwidth] -ysize [subst $${ns}::plotheight] \
            -radius [subst $${ns}::mradius] \
            -fillcolor [subst $${ns}::color] -marker [subst $${ns}::marker] \
            -xlabel "Residue #" \
            -nmwizns $ns \
            -plot
             
        } else {
          lappend ${ns}::plothandles [nmwiz_multiplot \
            -x [subst $${ns}::plotrids] -y [subst $${ns}::betalist] \
            -title "[subst $${ns}::title] square fluctuations" \
            -linewidth [subst $${ns}::linewidth] \
            -legend "Mode [subst $${ns}::activemode]" \
            -[subst $${ns}::lornol] -linecolor [subst $${ns}::color] \
            -xsize [subst $${ns}::plotwidth] -ysize [subst $${ns}::plotheight] \
            -radius [subst $${ns}::mradius] \
            -fillcolor [subst $${ns}::color] -marker [subst $${ns}::marker] \
            -xlabel "Residue #" \
            -nmwizns $ns \
            -plot]
        }
        vmdcon -info "Plot handle: [lindex [subst $${ns}::plothandles] end]"
        #-dash [subst $${ns}::dash] \
      }
      
      proc prepareSelmol {} {
        variable selid
        variable molid
        # make sure selmol exists
        if {$selid == -1 || [lsearch [molinfo list] $selid] == -1} {
          set selid [mol new]
          variable title
          variable w
          mol rename $selid "$title selections"
          $w.draw_arrows.plot_label configure -text "Plotting ($selid):"
        }
        # make sure coordinates are loaded
        if {[molinfo $selid get numframes] == 0} {
          variable molid
          set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]

          variable coordinates
          variable tempfn
          set outfile [open [file join $::nmwiz::tmpdir $tempfn] w]
          foreach line [[namespace current]::getPDBLines $coordinates] {
            puts $outfile $line
          } 
          close $outfile
          mol addfile [file join $::nmwiz::tmpdir $tempfn] molid $selid
          
          foreach id [molinfo list] {
            molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
          }
        }
        
        variable resids
        if {[molinfo $selid get numreps] != [llength $resids]} {
          for {set i [molinfo $selid get numreps]} {$i >= 0} {incr i -1} {
            mol delrep $i $selid
          }
          for {set i 0} {$i < [llength $resids]} {incr i} {
            mol addrep $selid
            mol modstyle $i $selid VDW
            mol showrep $selid $i off
          }
        }
        mol top $molid
      }

      proc clearSelection {} {
        variable selid
        if {$selid > -1 && [lsearch [molinfo list] $selid] > -1} {
          for {set i [expr [molinfo $selid get numreps] - 1]} {$i >= 0} {incr i -1} {
            mol showrep $selid $i off
          }
        }
        set labels [label list Atoms]
        for {set i [expr [llength $labels] - 1]} {$i >= 0} {incr i -1} {
          if {[lindex [lindex [lindex $labels $i] 0] 0] == $selid} {
            label delete Atoms $i
          }
        }
      }

      proc selectAtom {resid color} {

        variable plotrids
        set which [lsearch $plotrids $resid]
        if {$which == -1} {return 0}
        [namespace current]::prepareSelmol
        
        variable selid
        variable selection
        variable selectscale
        variable resolution
        variable material
        variable chainids
        variable resnames
        variable resids

        label add Atoms $selid/$which

        #set i [molinfo $selid get numreps]
        #mol addrep $selid
        if {[mol showrep $selid $which]} {
          mol showrep $selid $which off
          vmdcon -info "Deselected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
        } else {
          vmdcon -info "Selected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
          mol showrep $selid $which on
          mol modstyle $which $selid VDW $selectscale $resolution
          mol modmaterial $which $selid $material
          mol modselect $which $selid "index $which"
          mol modcolor $which $selid ColorID [lsearch "blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3" $color]
        }
      }


      proc updateProtRep {targetid} {
        variable molid

        if {[lsearch [molinfo list] $targetid] == -1} {
          if {$targetid == $molid} {
            [namespace current]::loadCoordinates
          } else {
            return 0
          }
        }
        variable showproteinas
        if {$showproteinas == "Custom" && $targetid == $molid} {
          return
        }
        variable tuberadius
        variable bondradius
        variable cutoffdistance
        variable resolution_protein
        variable material_protein
        variable betamin
        variable betamax
        variable spherescale
        variable proteincolor
        variable bfactors
        variable bfactormin
        variable bfactormax
        variable betalist
        variable selstr
        variable selrep
        
        for {set i [molinfo $targetid get numreps]} {$i >= 0} {incr i -1} {
          mol delrep $i $targetid
        }
        
        #if {$proteincolor != "Mobility" && $proteincolor != "Bfactors"} {
        #  set vmdcolorid [lsearch "blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3" $proteincolor]
        #} else {
        #  set vmdcolorid 0
        #}
        
        if {$proteincolor == "Bfactors"} {
          [atomselect $targetid "all"] set beta $bfactors
        } elseif {$proteincolor == "Mobility"} {
          [atomselect $targetid "all"] set beta $betalist
        } 
        
        if {$showproteinas == "Network"} {
          mol addrep $targetid
          mol modstyle 0 $targetid VDW $spherescale $resolution_protein
          mol modmaterial 0 $targetid $material_protein
          mol addrep $targetid
          mol modstyle 1 $targetid DynamicBonds $cutoffdistance $bondradius $resolution_protein
          mol modmaterial 1 $targetid $material_protein
          if {$proteincolor == "Mobility"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $betamin $betamax 
            mol modcolor 1 $targetid Beta
            mol scaleminmax $targetid 1 $betamin $betamax
            color scale midpoint 0.1 
          } elseif {$proteincolor == "Bfactors"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $bfactormin $bfactormax 
            mol modcolor 1 $targetid Beta
            mol scaleminmax $targetid 1 $bfactormin $bfactormax
            color scale midpoint 0.1 
          } else {
            mol modcolor 0 $targetid $proteincolor
            mol modcolor 1 $targetid $proteincolor
            #mol modcolor 0 $targetid ColorID $vmdcolorid
            #mol modcolor 1 $targetid ColorID $vmdcolorid
            color scale midpoint 0.5
          }
          if {$selrep} {
            mol modselect 0 $targetid $selstr
            mol modselect 1 $targetid $selstr
          }
        } else {
          mol addrep $targetid
          mol modstyle 0 $targetid $showproteinas $tuberadius $resolution_protein
          mol modmaterial 0 $targetid $material_protein
          if {$selrep} {
            mol modselect 0 $targetid $selstr
          }
          if {$proteincolor == "Mobility"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $betamin $betamax
            color scale midpoint 0.1 
          } elseif {$proteincolor == "Bfactors"} {
            mol modcolor 0 $targetid Beta
            mol scaleminmax $targetid 0 $bfactormin $bfactormax
            color scale midpoint 0.1 
          } else {
            mol modcolor 0 $targetid $proteincolor
            #mol modcolor 0 $targetid ColorID $vmdcolorid
            color scale midpoint 0.5
          }
        }
      }

      proc calcMSF {} {
        variable molid
        variable activemode
        variable indices
        variable lengths
        variable scalearrows
        variable modes
        variable animid
        variable material
        variable selstr
        variable ndim

        set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [lindex $modes [lsearch $indices $activemode]]
        variable msformode
        if {$msformode == "Mobility"} {
          set mode [vecscale [expr $length * $length] [vecmul $mode $mode]]  
        } else {
          set mode [vecscale $length $mode]
        }
        

        set index 0
        variable betalist {}
        variable betamin 10000
        variable betamax -10000
        if {$ndim == 3} {
          foreach {mx my mz} $mode {
            set beta [expr $mx + $my + $mz]
            lappend betalist $beta
            if {$beta < $betamin} {set betamin $beta}
            if {$beta > $betamax} {set betamax $beta}
            incr index
          }
        } else {
          foreach mx $mode {
            set beta [expr $mx]
            lappend betalist $beta
            if {$beta < $betamin} {set betamin $beta}
            if {$beta > $betamax} {set betamax $beta}
            incr index
          }
        }
        [atomselect $molid "all"] set beta $betalist 
        #mol scaleminmax $molid 0 [::tcl::mathfunc::min $betalist] [::tcl::mathfunc::max $betalist]
        color scale midpoint 0.1
        color scale method BWR
        
      }

      proc drawAction {} {
        variable overwrite
        variable arrid
        variable arrids
        if {!$overwrite} {
          set arrid [mol new]
          lappend arrids $arrid
        }
        [namespace current]::drawArrows
      }
      
      proc autoUpdate {} {
        variable autoupdate
        if {$autoupdate} {
          variable overwrite
          variable arrid
          variable arrids
          [namespace current]::drawArrows
        }
      }
      
      proc evalRMSD {} {
        variable rmsd
        if {![string is double $rmsd]} {
          set rmsd 2
        } elseif {$rmsd <= 0} {
          set rmsd 0.1
        }
        variable scalearrows
        variable lengths
        variable modes
        variable activemode
        variable indices
        variable n_atoms
        set whichmode [lsearch $indices $activemode]
        set length [lindex $lengths $whichmode]
        set scalearrows [expr [sign $scalearrows] * $rmsd / $length / [veclength [lindex $modes $whichmode]] * $n_atoms ** 0.5 ]
        vmdcon -info "Mode $activemode is scaled by [format %.2f $scalearrows](x[format %.2f $length]) for [format %.2f $rmsd] A RMSD."
        variable rmsdprev $rmsd
      }

      proc drawArrows {} {
        variable rmsd
        variable rmsdprev
        variable activeindex
        variable activeindexprev
        if {$rmsdprev != $rmsd || $activeindex!= $activeindexprev} {[namespace current]::evalRMSD}
        variable color
        variable material
        variable resolution
        variable coordinates
        variable scalearrows
        variable sense
        variable arrid
        variable arrids
        variable activemode
        variable indices
        variable modes
        variable title
        variable prefix
        variable lengths
        variable molid
        variable w
        variable selstr
        
        variable hide_shorter 
        variable cylinder_radius
        variable cone_radius
        variable cone_height
        
        set whichmode [lsearch $indices $activemode] 

        if {[lsearch [molinfo list] $molid] == -1} {
          [namespace current]::loadCoordinates
        }

        if {[lsearch [molinfo list] $arrid] == -1} {
          set arrid [mol new]
          lappend arrids $arrid
        } else {
          graphics $arrid delete all
        }
        graphics $arrid color $color
        graphics $arrid materials on
        graphics $arrid material $material
        set length [lindex $lengths $whichmode]
        set mode [vecscale [expr $length * $scalearrows] [lindex $modes $whichmode]]
        
        variable bothdirections
  
        foreach index [[atomselect $molid $selstr] get index] {
          set from [expr $index * 3]
          set to  [expr $from + 2]
          set xyz [lrange $coordinates $from $to ] 
          set v [lrange $mode $from $to ] 
          if {$hide_shorter < [veclength $v]} {
            set temp [vecadd $xyz $v]
            graphics $arrid cylinder $xyz $temp radius $cylinder_radius resolution $resolution
            set temp2 [vecadd $temp [vecscale $v [expr $cone_height / [veclength $v]]]]
            graphics $arrid cone $temp $temp2 radius $cone_radius resolution $resolution
            
            if {$bothdirections} {
              set v [vecscale -1 $v]
              set temp [vecadd $xyz $v]
              graphics $arrid cylinder $xyz $temp radius $cylinder_radius resolution $resolution
              set temp2 [vecadd $temp [vecscale $v [expr $cone_height / [veclength $v]]]]
              graphics $arrid cone $temp $temp2 radius $cone_radius resolution $resolution
            }
          }
        }
        mol rename $arrid "$title mode $activemode arrows"
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]
        display resetview
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
        [namespace current]::calcMSF
        $w.draw_arrows.arrowbuttons_label configure -text "Arrows ($arrid):"
        
        variable arridlist
        lset arridlist $whichmode $arrid
      }

      proc getPDBLines {coords} {
        variable atomnames
        variable resnames
        variable chainids
        variable resids
        variable betalist
        set pdblines ""
        set i 0
        foreach an $atomnames rn $resnames ci $chainids ri $resids {x y z} $coords b $betalist {
          incr i
          if {[string length $an] < 4} {
            #lappend pdblines [format "ATOM  %5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
            #                     $i  $an $rn  $ci $ri $x $y $z 1.0 $b]
            lappend pdblines [format "ATOM  %5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f" \
                                 $i  $an [string range $rn 0 4] $ci $ri $x $y $z]
          } else {
            set an [string range $an 0 3]
            lappend pdblines [format "ATOM  %5d %4s %-4s%1s%4d    %8.3f%8.3f%8.3f" \
                                 $i  $an [string range $rn 0 4] $ci $ri $x $y $z]
          }
        }
        return $pdblines
      }

      proc locateCoordinates {} {
        variable pdbfile
        variable molid
        if {$molid > -1} {
          tk_messageBox -type ok -title "WARNING" \
            -message "Coordinate data from $pdbfile is already loaded."
          return
        }
        
        set tempfile [tk_getOpenFile \
          -filetypes {{"PDB files" { .pdb .PDB }} {"All files" *}}]
        if {![string equal $tempfile ""]} { 
          set pdbfile $tempfile 
        }
        [namespace current]::loadCoordinates
      }


      proc Animate {} {
        variable nframes
        variable activemode
        variable coordinates
        variable prefix
        variable lengths
        variable indices
        variable scalearrows
        variable modes
        variable animid
        variable title
        variable molid
        variable autoplay
        variable w
        variable tempfn
        variable selstr
        variable selrep
        variable betalist
        variable betamin
        variable betamax
        puts [namespace current]
        set whichmode [lsearch $indices $activemode] 
        animate pause
        set animfn [file join $::nmwiz::tmpdir $tempfn]
        [namespace current]::loadCoordinates
        variable molid
        set sel [atomselect $molid "all"]
        $sel writepdb $animfn
        $sel delete
        
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]
        if {[lsearch [molinfo list] $animid] == -1} {
          set animid [mol new $animfn]
        } else {
          animate delete beg 0 end -1 skip 0 $animid
          mol addfile $animfn waitfor all $animid
        }
        mol delrep 0 $animid
        mol off $animid

        mol top $animid
        $w.draw_arrows.animbuttons_label configure -text "Animation ($animid):"
        mol rename $animid "$title mode $activemode animation"
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
        for {set i 1} {$i <= $nframes} {incr i} {
          
        }

        set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [vecscale [expr $length * [::tcl::mathfunc::abs $scalearrows]] [lindex $modes [lsearch $indices $activemode]]]
        set coords [vecadd $coordinates $mode]
        set mode [::nmwiz::array2xyz [vecscale $mode [expr  -2.0 / $nframes]]]
        set sel [atomselect $animid "all"]
        $sel set {x y z} [::nmwiz::array2xyz $coords]
        for {set i 1} {$i <= $nframes} {incr i} {
          animate dup frame [expr $i - 1] $animid
          $sel frame $i
          $sel lmoveby $mode

        }
        $sel delete
        [namespace current]::updateProtRep $animid 
        mol off $animid
        if {$autoplay} {
          animate speed 0.96
          animate style rock
          animate forward
        }
        mol on $animid
        if {$selrep} {
          mol modselect 0 $animid $selstr
        }
        eval "\$[namespace current]::w.draw_arrows.animbuttons_showhide configure -text Hide"
        eval "\$[namespace current]::w.draw_arrows.animbuttons_stop configure -text Pause"
        set [namespace current]::stopped 0
        
        variable animidlist
        lset animidlist $whichmode $animid

      }


      proc loadCoordinates {} {
        variable molid
        variable coordinates
        variable title
        variable w
        variable tempfn
        if {[lsearch [molinfo list] $molid] != -1} {
          return 0
        }
        
        set outfile [open [file join $::nmwiz::tmpdir $tempfn] w]
        foreach line [[namespace current]::getPDBLines $coordinates] {
          puts $outfile $line
        } 
        close $outfile
        if {[molinfo num] > 0 && $::nmwiz::preserview} {
          set currentview [molinfo [lindex [molinfo list] 0] get {rotate_matrix center_matrix scale_matrix global_matrix}]
        }

        set molid [mol new [file join $::nmwiz::tmpdir $tempfn]]

        $w.draw_arrows.protbuttons_label configure -text "Protein ($molid):"
        mol rename $molid "$title coordinates"
        [namespace current]::calcMSF
        [namespace current]::updateProtRep $molid

        if {[molinfo num] > 0 && $::nmwiz::preserview} {
          foreach id [molinfo list] {
            molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
          }
        }
      }

      proc checkCoordinates {} {
        variable molid
        variable modes
        variable coordinates
        if {[[atomselect $molid "all"] num] != [llength [lindex $modes 0]]} {
          mol delete $molid
          set molid -1
          if {"ok" == [tk_messageBox -type okcancel -title "ERROR" \
              -message "[[atomselect $molid all] num] atoms are loaded. Coordinate data file must contain [llength [lindex $modes 0]] atoms. Please locate the correct file."]} {
            [namespace current]::locateCoordinates
          } 
        } else {
          set coordinates [[atomselect $molid all] get {x y z}]
        }
      }
      
      proc changeColor {} {
        variable color
        variable colorlist
        variable indices
        variable activemode
        lset colorlist [lsearch $indices $activemode] $color
      }
      
      proc prevMode {} {
        variable activeindex
        if {$activeindex > 0} {
          variable activemode
          variable indices
          set activemode [lindex $indices [expr $activeindex - 1]];
          [namespace current]::changeMode;
        }
      }
      
      proc nextMode {} {
        variable activeindex
        variable numofmodes
        if {$activeindex < [expr $numofmodes - 1]} {
          variable activemode
          variable indices
          set activemode [lindex $indices [expr $activeindex + 1]];
          [namespace current]::changeMode;
        }
      }

      
      proc changeMode {} {
        variable w
        variable activemode
        variable indices
        variable arrid
        variable arridlist
        variable animid
        variable animidlist
        variable color
        variable colorlist
        variable drawlengthstr
        variable lengths
        variable activeindex 
        variable ndim
        variable activeindexprev $activeindex
        set inactiveindex $activeindex 
        set activeindex [lsearch $indices $activemode]
        if {$ndim == 3} {
          variable scalearrows_list
          variable rmsd_list
          variable hide_shorter_list 
          variable cylinder_radius_list
          variable cone_radius_list
          variable cone_height_list
          variable resolution_list
          variable material_list
          
          variable scalearrows
          variable hide_shorter 
          variable cylinder_radius
          variable cone_radius
          variable cone_height
          variable resolution
          variable material
          variable rmsd
          
          lset scalearrows_list $inactiveindex $scalearrows
          lset hide_shorter_list $inactiveindex $hide_shorter  
          lset cylinder_radius_list $inactiveindex $cylinder_radius
          lset cone_radius_list $inactiveindex $cone_radius
          lset cone_height_list $inactiveindex $cone_height
          lset resolution_list $inactiveindex $resolution
          lset material_list $inactiveindex $material
          lset rmsd_list $inactiveindex $rmsd

          set drawlengthstr [format "%.1f" [lindex $lengths $activeindex]];
          
          set scalearrows [lindex $scalearrows_list $activeindex]
          set hide_shorter [lindex $hide_shorter_list $activeindex]
          set cylinder_radius [lindex $cylinder_radius_list $activeindex]
          set cone_radius [lindex $cone_radius_list $activeindex]
          set cone_height [lindex $cone_height_list $activeindex]
          set resolution [lindex $resolution_list $activeindex]
          set material [lindex $material_list $activeindex]
          set rmsd [lindex $rmsd_list $activeindex]

      
          variable overwrite
          if {$overwrite} {
            if {$arrid > -1 && [lsearch [molinfo list] $arrid] > -1} {
              mol off $arrid
            }
          }
          if {$animid > -1 && [lsearch [molinfo list] $animid] > -1} {
            mol off $animid
          }
          
          set which [lsearch $indices $activemode]
          set arrid [lindex $arridlist $which]      
          set animid [lindex $animidlist $which]

          
          set color [lindex $colorlist $which]

          if {$arrid > -1 && [lsearch [molinfo list] $arrid] > -1} {
            mol on $arrid
            [namespace current]::calcMSF
          } else {
            [namespace current]::drawAction
          }
          
          $w.draw_arrows.arrowbuttons_showhide configure -text Hide

          if {$animid > -1 && [lsearch [molinfo list] $animid] > -1} {
            mol on $animid
            mol top $animid
          } else {
            variable autoanimate  
            if {$autoanimate} {
              [namespace current]::Animate
            }
          }
          $w.draw_arrows.animbuttons_showhide configure -text Hide
          $w.draw_arrows.animbuttons_stop configure -text Play
        } else {
          [namespace current]::calcMSF          
        }
        
      }
      
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      
      proc nmwizgui {} {
        variable w
        set ns [namespace current]
        if [winfo exists .[string range $ns 2 end]] {
          wm deiconify .[string range $ns 2 end]
          raise .[string range $ns 2 end]
          return
        }
        set w [toplevel .[string range $ns 2 end]]
        variable title

        wm title $w "NMWiz - $title"
        wm resizable $w 0 0
        set wam [labelframe $w.active_mode -text "$title" -bd 2] 

        variable ndim
        grid [label $wam.active_label -text "Active mode:"] \
            -row 0 -column 1 -sticky w

        grid [frame $wam.active] \
            -row 0 -column 2 -columnspan 3 -sticky w

        tk_optionMenu $wam.active.list ${ns}::activemode 0
        $wam.active.list.menu delete 0
        variable indices
        variable lengths
        foreach index $indices length $lengths {
          $wam.active.list.menu add radiobutton -label $index \
              -variable ${ns}::activemode \
              -command "${ns}::changeMode;"
        }
        
        button $wam.active.prev -text "<=" -command "${ns}::prevMode"
        button $wam.active.negate -text "+/-" -command "set ${ns}::scalearrows \[expr - \$${ns}::scalearrows]; ${ns}::autoUpdate"  
        button $wam.active.next -text "=>" -command "${ns}::nextMode"
        
        variable ndim
        if {$ndim == 3} {
          tk_optionMenu $wam.active.color ${ns}::color "blue"
          $wam.active.color.menu delete 0
          foreach acolor "blue red gray orange yellow tan green white pink \
        cyan purple black yellow2 yellow3 green2 green3 \
        cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
        orange3" {
            $wam.active.color.menu add radiobutton -label $acolor \
                -variable ${ns}::color \
                -command "${ns}::changeColor; ${ns}::autoUpdate; "
          }
          pack $wam.active.list $wam.active.prev $wam.active.negate $wam.active.next $wam.active.color -side left -anchor w -fill x
        } else {
          tk_optionMenu $wam.active.color ${ns}::msformode "Mobility"
          $wam.active.color.menu delete 0
          foreach acolor "Mobility Eigenvector" {
            $wam.active.color.menu add radiobutton -label $acolor \
                -variable ${ns}::msformode \
                -command "${ns}::calcMSF;"
          }
          pack $wam.active.list $wam.active.prev $wam.active.next $wam.active.color -side left -anchor w -fill x
        }
        
        
        
        #blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3

        if {$ndim == 3} {
          grid [label $wam.scale_label -text "RMSD (A):"] \
            -row 2 -column 1 -sticky w
          grid [frame $wam.adjust_frame] \
            -row 2 -column 2 -columnspan 3 -sticky w
          entry $wam.adjust_frame.entry -width 4 -textvariable ${ns}::rmsd
          button $wam.adjust_frame.incr1 -text "+0.1" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd + 0.1]; ${ns}::autoUpdate"
          button $wam.adjust_frame.incr5 -text "+0.5" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd + 0.5]; ${ns}::autoUpdate"
          button $wam.adjust_frame.decr5 -text "-0.5" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd - 0.5]; ${ns}::autoUpdate"
          button $wam.adjust_frame.decr1 -text "-0.1" -command \
            "set ${ns}::rmsd \[expr \$${ns}::rmsd - 0.1]; ${ns}::autoUpdate"
          pack $wam.adjust_frame.entry \
            $wam.adjust_frame.incr1 $wam.adjust_frame.incr5 \
            $wam.adjust_frame.decr5 $wam.adjust_frame.decr1 \
            -side left -anchor w -fill x

          grid [label $wam.selstr_label -text "Selection:"] \
            -row 5 -column 1 -sticky w
          grid [entry $wam.selstr_entry \
            -textvariable ${ns}::selstr] \
            -row 5 -column 2 -columnspan 2 -sticky we
          grid [button $wam.selstr_draw -width 4 -text "Redraw" \
              -command ${ns}::drawArrows] \
            -row 5 -column 4 -sticky we
        }
        grid [button $wam.showmain -text "Main" \
            -command nmwiz_tk] \
          -row 8 -column 1 -sticky we
        grid [button $wam.save -text "Save" \
            -command "::nmwiz::writeNMD $ns"] \
          -row 8 -column 2 -sticky we
        grid [button $wam.remove -text "Remove" \
            -command "lset ::nmwiz::titles $::nmwiz::guicount NONE; pack forget .nmwizgui.{[string range $ns 2 end]}frame; ${ns}::deleteMolecules; namespace delete $ns; destroy .[string range $ns 2 end]"] \
          -row 8 -column 3 -sticky we
        grid [button $wam.showhelp -text "Help" \
            -command {::nmwiz::showHelp wizard}] \
          -row 8 -column 4 -sticky we

        pack $wam -side top -ipadx 10 -ipady 5 -fill x -expand 1

        set wda [labelframe $w.draw_arrows -text "Actions" -bd 2]
        
        if {$ndim == 3} {
          grid [label $wda.arrowbuttons_label -text "Arrows:"] \
            -row 5 -column 0 -sticky w
          grid [button $wda.arrowbuttons_draw -text "Draw" \
              -command ${ns}::drawArrows] \
            -row 5 -column 2 -sticky ew
          grid [button $wda.arrowbuttons_clean -text "Clean" \
              -command "foreach anarrid \$${ns}::arrids {if {\$anarrid != \$${ns}::arrid && \[lsearch \[molinfo list] \$${ns}::arrid] != -1} {mol delete \$anarrid}; if {\[lsearch \[molinfo list] \$${ns}::arrid] != -1} {graphics \$${ns}::arrid delete all}}"] \
            -row 5 -column 3
          grid [button $wda.arrowbuttons_showhide -text "Hide" \
              -command "if {\[molinfo \$${ns}::arrid get displayed]} {mol off \$${ns}::arrid;\
                        \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Show} else {mol on \$${ns}::arrid;\
                        \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Hide}"] \
            -row 5 -column 4 -sticky ew
          grid [button $wda.arrowbuttons_options -text "Options" \
              -command "if {\$${ns}::arropt} {pack forget \$${ns}::w.graphics_options;\
                        set ${ns}::arropt 0; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief raised} else {pack \$${ns}::w.graphics_options -side top -ipadx 10 -ipady 5 -fill x -expand 1;\
                        set ${ns}::arropt 1; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief sunken}"] \
            -row 5 -column 5 -sticky ew

          grid [label $wda.animbuttons_label -text "Animation:"] \
            -row 6 -column 0 -sticky w
          grid [button $wda.animbuttons_animate -text "Make" \
              -command ${ns}::Animate] \
            -row 6 -column 2 -sticky ew
          grid [button $wda.animbuttons_stop -text "Play" \
              -command "if {\$${ns}::animid == -1} {${ns}::Animate} else {if {\$${ns}::stopped} {mol top \$${ns}::animid; animate forward; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Pause; set ${ns}::stopped 0} else {animate pause; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Play; set ${ns}::stopped 1}}"] \
            -row 6 -column 3 -sticky ew
          grid [button $wda.animbuttons_showhide -text "Hide" \
              -command "if {\$${ns}::animid > -1 && \[lsearch \[molinfo list] \$${ns}::animid] > -1} {if {\[molinfo \$${ns}::animid get displayed]} {animate pause; mol off \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Show} else {mol on \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Hide; animate forward}}"] \
            -row 6 -column 4 -sticky ew
          grid [button $wda.animbuttons_options -text "Options" \
              -command "if {\$${ns}::anmopt} {pack forget \$${ns}::w.animation_options; set ${ns}::anmopt 0; \$${ns}::w.draw_arrows.animbuttons_options configure -relief raised} else {pack \$${ns}::w.animation_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::anmopt 1; \$${ns}::w.draw_arrows.animbuttons_options configure -relief sunken}"] \
            -row 6 -column 5 -sticky ew
        }

        grid [label $wda.plot_label -text "Plotting:"] \
          -row 8 -column 0 -sticky w
        grid [button $wda.plot_plot -text "Plot" \
            -command "${ns}::Plot ${ns}"] \
          -row 8 -column 2 -sticky ew
        grid [button $wda.plot_clear -text "Clear" \
            -command "${ns}::clearSelection"] \
          -row 8 -column 3 -sticky ew
        grid [button $wda.plot_showhide -text "Hide" \
            -command "if {\$${ns}::selid > -1 && \[lsearch \[molinfo list] \$${ns}::selid] > -1} {if {\[molinfo \$${ns}::selid get displayed]} {mol off \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Show} else {mol on \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Hide}}"] \
          -row 8 -column 4 -sticky ew
        grid [button $wda.plot_options -text "Options" \
            -command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.plotting_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plotting_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \
          -row 8 -column 5 -sticky ew
         
        ##-command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.animation_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plot_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \

        grid [label $wda.protbuttons_label -text "Protein:"] \
          -row 9 -column 0 -sticky w
        grid [button $wda.prt_update -text "Update" \
            -command "${ns}::updateProtRep \$${ns}::molid"] \
          -row 9 -column 2 -sticky ew
        grid [button $wda.protbuttons_focus -text "Focus" \
            -command "mol top \$${ns}::molid; display resetview"] \
          -row 9 -column 3  -sticky ew
        grid [button $wda.protbuttons_showhide -text "Hide" \
            -command "if {\[molinfo \$${ns}::molid get displayed]} {mol off \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Show;} else {mol on \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Hide;}"] \
          -row 9 -column 4 -sticky ew
        grid [button $wda.protbuttons_repoptions -text "Options" \
            -command "if {\$${ns}::prtopt} {pack forget \$${ns}::w.prograph_options; set ${ns}::prtopt 0; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief raised} else {pack \$${ns}::w.prograph_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::prtopt 1; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief sunken}"] \
          -row 9 -column 5 -sticky ew

        pack $wda -side top -fill x -expand 1

        set wgo [labelframe $w.graphics_options -text "Arrow Graphics Options" -bd 2]
        
        grid [checkbutton $wgo.auto_check -text "auto update arrows" \
            -variable ${ns}::autoupdate] \
          -row 0 -column 1 -sticky w

        grid [checkbutton $wgo.overwrite_check -text "auto hide inactive mode" \
            -variable ${ns}::overwrite] \
          -row 0 -column 2 -sticky w

        grid [checkbutton $wgo.both_check -text "draw in both directions" \
            -variable ${ns}::bothdirections] \
          -row 7 -column 1 -columnspan 2 -sticky w

        grid [label $wgo.hide_label -text "Draw if longer than:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wgo.hide_frame] \
          -row 9 -column 2 -sticky w
        entry $wgo.hide_frame.entry -width 4 -textvariable ${ns}::hide_shorter
        button $wgo.hide_frame.decr -text "-0.5" \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter - 0.5]]; ${ns}::autoUpdate"
        button $wgo.hide_frame.incr -text "+0.5" \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter + 0.5]]; ${ns}::autoUpdate"
        label $wgo.hide_frame.angstrom -text "A"
        pack $wgo.hide_frame.entry $wgo.hide_frame.decr $wgo.hide_frame.incr \
          $wgo.hide_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.cylinder_label -text "Arrow cylinder radius:"] \
          -row 10 -column 1 -sticky w
        grid [frame $wgo.cylinder_frame] \
          -row 10 -column 2 -sticky w
        entry $wgo.cylinder_frame.entry -width 4 -textvariable ${ns}::cylinder_radius
        button $wgo.cylinder_frame.decr -text "-0.1" \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius - 0.1]]; ${ns}::autoUpdate"
        button $wgo.cylinder_frame.incr -text "+0.1" \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius + 0.1]]; ${ns}::autoUpdate"
        label $wgo.cylinder_frame.angstrom -text "A"
        pack $wgo.cylinder_frame.entry $wgo.cylinder_frame.decr \
          $wgo.cylinder_frame.incr $wgo.cylinder_frame.angstrom \
          -side left -anchor w -fill x

        grid [label $wgo.coner_label -text "Arrow cone radius:"] \
          -row 11 -column 1 -sticky w
        grid [frame $wgo.coner_frame] \
          -row 11 -column 2 -sticky w
        entry $wgo.coner_frame.entry -width 4 -textvariable ${ns}::cone_radius
        button $wgo.coner_frame.decr -text "-0.1" \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius - 0.1]]; ${ns}::autoUpdate"
        button $wgo.coner_frame.incr -text "+0.1" \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius + 0.1]]; ${ns}::autoUpdate"
        label $wgo.coner_frame.angstrom -text "A"
        pack $wgo.coner_frame.entry $wgo.coner_frame.decr $wgo.coner_frame.incr \
          $wgo.coner_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.coneh_label -text "Arrow cone height:"] \
          -row 12 -column 1 -sticky w
        grid [frame $wgo.coneh_frame] \
          -row 12 -column 2 -sticky w
        entry $wgo.coneh_frame.entry -width 4 -textvariable ${ns}::cone_height
        button $wgo.coneh_frame.decr -text "-0.2" \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height - 0.2]]; ${ns}::autoUpdate"
        button $wgo.coneh_frame.incr -text "+0.2" \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height + 0.2]]; ${ns}::autoUpdate"
        label $wgo.coneh_frame.angstrom -text "A"
        pack $wgo.coneh_frame.entry $wgo.coneh_frame.decr $wgo.coneh_frame.incr \
          $wgo.coneh_frame.angstrom -side left -anchor w -fill x

        grid [label $wgo.material_label -text "Graphics material:"] \
          -row 20 -column 1 -sticky w
        grid [frame $wgo.material_frame] \
          -row 20 -column 2 -sticky w
        tk_optionMenu $wgo.material_frame.list ${ns}::material "Opaque"
        $wgo.material_frame.list.menu delete 0
        foreach mtrl "Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy" {
          $wgo.material_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::material \
              -command "${ns}::autoUpdate"    
        }
        pack $wgo.material_frame.list -side left -anchor w -fill x

        grid [label $wgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wgo.resolution_frame.list ${ns}::resolution 6 
        $wgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution \
              -command "${ns}::autoUpdate"  
        } 
        pack $wgo.resolution_frame.list -side left -anchor w -fill x

        set wpgo [labelframe $w.prograph_options -text "Protein Graphics Options" -bd 2]
        
        grid [checkbutton $wpgo.selstr_check -text "show selected atoms" \
            -variable ${ns}::selrep -command "${ns}::autoUpdate"] \
          -row 0 -column 1 -columnspan 2 -sticky w
        
        grid [label $wpgo.protas_label -text "Show protein as:"] \
          -row 13 -column 1 -sticky w
        grid [frame $wpgo.protas_frame] \
          -row 13 -column 2 -sticky w
        tk_optionMenu $wpgo.protas_frame.list ${ns}::showproteinas "Tube"
        $wpgo.protas_frame.list.menu delete 0
        $wpgo.protas_frame.list.menu add radiobutton -label "Custom" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Licorice" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Network" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Trace" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Tube" -variable ${ns}::showproteinas -command "${ns}::updateProtRep \$${ns}::molid"
        pack $wpgo.protas_frame.list -side left -anchor w -fill x

        grid [label $wpgo.procolor_label -text "Color protein:"] \
          -row 14 -column 1 -sticky w
        grid [frame $wpgo.procolor_frame] \
          -row 14 -column 2 -sticky w
        tk_optionMenu $wpgo.procolor_frame.list ${ns}::proteincolor "Mobility"
        $wpgo.procolor_frame.list.menu delete 0
        $wpgo.procolor_frame.list.menu add radiobutton -label "Mobility" -variable ${ns}::proteincolor -command "${ns}::updateProtRep \$${ns}::molid"
        $wpgo.procolor_frame.list.menu add radiobutton -label "Bfactors" -variable ${ns}::proteincolor -command "${ns}::updateProtRep \$${ns}::molid"
        foreach acolor "Index Chain ResName ResType" {
          $wpgo.procolor_frame.list.menu add radiobutton -label $acolor \
              -variable ${ns}::proteincolor \
              -command "${ns}::updateProtRep \$${ns}::molid"
        }
        pack $wpgo.procolor_frame.list -side left -anchor w -fill x


        grid [label $wpgo.cutoffdistance_label -text "Network cutoff distance:"] \
          -row 16 -column 1 -sticky w
        grid [frame $wpgo.cutoffdistance_frame] \
          -row 16 -column 2 -sticky w
        entry $wpgo.cutoffdistance_frame.entry -width 4 -textvariable ${ns}::cutoffdistance
        button $wpgo.cutoffdistance_frame.decr -text "-1.0" \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance - 1.0]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.cutoffdistance_frame.incr -text "+1.0" \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance + 1.0]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.cutoffdistance_frame.angstrom -text "A"
        pack $wpgo.cutoffdistance_frame.entry $wpgo.cutoffdistance_frame.decr $wpgo.cutoffdistance_frame.incr \
          $wpgo.cutoffdistance_frame.angstrom -side left -anchor w -fill x

        grid [label $wpgo.nodescale_label -text "Scale node spheres:"] \
          -row 17 -column 1 -sticky w
        grid [frame $wpgo.nodescale_frame] \
          -row 17 -column 2 -sticky w
        entry $wpgo.nodescale_frame.entry -width 4 -textvariable ${ns}::spherescale
        button $wpgo.nodescale_frame.decr -text "-0.1" \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.nodescale_frame.incr -text "+0.1" \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        pack $wpgo.nodescale_frame.entry $wpgo.nodescale_frame.decr $wpgo.nodescale_frame.incr \
          -side left -anchor w -fill x
          
        grid [label $wpgo.bondradius_label -text "Network bond radius:"] \
          -row 18 -column 1 -sticky w
        grid [frame $wpgo.bondradius_frame] \
          -row 18 -column 2 -sticky w
        entry $wpgo.bondradius_frame.entry -width 4 -textvariable ${ns}::bondradius
        button $wpgo.bondradius_frame.decr -text "-0.1" \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.bondradius_frame.incr -text "+0.1" \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.bondradius_frame.angstrom -text "A"
        pack $wpgo.bondradius_frame.entry $wpgo.bondradius_frame.decr $wpgo.bondradius_frame.incr \
          $wpgo.bondradius_frame.angstrom -side left -anchor w -fill x


        grid [label $wpgo.tuberadius_label -text "Tube/licorice radius:"] \
          -row 19 -column 1 -sticky w
        grid [frame $wpgo.tuberadius_frame] \
          -row 19 -column 2 -sticky w
        entry $wpgo.tuberadius_frame.entry -width 4 -textvariable ${ns}::tuberadius
        button $wpgo.tuberadius_frame.decr -text "-0.1" \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius - 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        button $wpgo.tuberadius_frame.incr -text "+0.1" \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius + 0.1]]; ${ns}::updateProtRep \$${ns}::molid"
        label $wpgo.tuberadius_frame.angstrom -text "A"
        pack $wpgo.tuberadius_frame.entry $wpgo.tuberadius_frame.decr $wpgo.tuberadius_frame.incr \
          $wpgo.tuberadius_frame.angstrom -side left -anchor w -fill x


        grid [label $wpgo.material_label -text "Graphics material:"] \
          -row 20 -column 1 -sticky w
        grid [frame $wpgo.material_frame] \
          -row 20 -column 2 -sticky w
        tk_optionMenu $wpgo.material_frame.list ${ns}::material_protein "Opaque"
        $wpgo.material_frame.list.menu delete 0
        foreach mtrl "Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy" {
          $wpgo.material_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::material_protein \
              -command "${ns}::updateProtRep \$${ns}::molid"    
        }
        pack $wpgo.material_frame.list -side left -anchor w -fill x

        grid [label $wpgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wpgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wpgo.resolution_frame.list ${ns}::resolution_protein 6 
        $wpgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wpgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution_protein \
              -command "${ns}::updateProtRep \$${ns}::molid"  
        } 
        pack $wpgo.resolution_frame.list -side left -anchor w -fill x

        set wao [labelframe $w.animation_options -text "Animation Options" -bd 2]
        
        grid [checkbutton $wao.auto_check -text "auto animate" \
            -variable ${ns}::autoanimate] \
          -row 0 -column 1 -columnspan 2 -sticky w
        
        grid [checkbutton $wao.autoplay_check -text "continuous autoplay" -variable ${ns}::autoplay] \
          -row 0 -column 3 -columnspan 2 -sticky w
          
        grid [label $wao.nframes_label -text "Number of frames:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wao.nframes_frame] \
          -row 9 -column 2 -columnspan 3 -sticky w
        entry $wao.nframes_frame.entry -width 4 -textvariable ${ns}::nframes
        button $wao.nframes_frame.decr5 -text "-5" \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes - 5]"
        button $wao.nframes_frame.incr5 -text "+5" \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes + 5]"
        pack $wao.nframes_frame.entry $wao.nframes_frame.decr5 \
          $wao.nframes_frame.incr5 -side left -anchor w -fill x

        set wpo [labelframe $w.plotting_options -text "Plotting Options" -bd 2]

        grid [checkbutton $wpo.overplot_check -text "reuse figure for new plots" \
            -variable ${ns}::overplot] \
          -row 0 -column 1 -columnspan 2 -sticky w
  
        if {$ndim == 1} {        
          grid [label $wpo.plotcolor_label -text "Line color:"] \
            -row 0 -column 5 -sticky w
            
          grid [frame $wpo.pcf] \
            -row 0 -column 6 -sticky w
          tk_optionMenu $wpo.pcf.color ${ns}::color "blue"
            $wpo.pcf.color.menu delete 0
            foreach acolor "blue red gray orange yellow tan green white pink \
          cyan purple black yellow2 yellow3 green2 green3 \
          cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
          orange3" {
              $wpo.pcf.color.menu add radiobutton -label $acolor \
                  -variable ${ns}::color
            }
          pack $wpo.pcf.color -side left -anchor w -fill x
        }
        grid [label $wpo.plotwidth_label -text "Plot width:"] \
          -row 1 -column 1 -sticky w
        grid [entry $wpo.plotwidth_entry -width 4 -textvariable ${ns}::plotwidth] \
          -row 1 -column 2 -sticky w
        grid [label $wpo.spacing_label -text "  "] \
          -row 1 -column 3 -sticky w

        grid [label $wpo.plotheight_label -text "Plot height:"] \
          -row 1 -column 5 -sticky w
        grid [entry $wpo.plotheight_entry -width 4 -textvariable ${ns}::plotheight] \
          -row 1 -column 6 -sticky w

        grid [label $wpo.line_label -text "Lines:"] \
          -row 3 -column 1 -sticky w
        grid [frame $wpo.line_frame] \
          -row 3 -column 2 -sticky w
        tk_optionMenu $wpo.line_frame.list ${ns}::lornol "lines"
        $wpo.line_frame.list.menu delete 0
        foreach lnl "lines nolines" {
          $wpo.line_frame.list.menu add radiobutton -label $lnl \
              -variable ${ns}::lornol 
        }
        pack $wpo.line_frame.list -side left -anchor w -fill x  

        grid [label $wpo.linewidth_label -text "Line width:"] \
          -row 3 -column 5 -sticky w
        grid [entry $wpo.linewidth_entry -width 4 -textvariable ${ns}::linewidth] \
          -row 3 -column 6 -sticky w
          
        grid [label $wpo.marker_label -text "Marker:"] \
          -row 5 -column 1 -sticky w
        grid [frame $wpo.marker_frame] \
          -row 5 -column 2 -sticky w
        tk_optionMenu $wpo.marker_frame.list ${ns}::marker "circle"
        $wpo.marker_frame.list.menu delete 0
        foreach mrkr "none point circle square" {
          $wpo.marker_frame.list.menu add radiobutton -label $mrkr \
              -variable ${ns}::marker     
        }
        pack $wpo.marker_frame.list -side left -anchor w -fill x
        
        grid [label $wpo.radius_label -text "Marker size:"] \
          -row 5 -column 5 -sticky w
        grid [entry $wpo.radius_entry -width 4 -textvariable ${ns}::mradius] \
          -row 5 -column 6 -sticky w

        #grid [button $wpo.dash_help -text "?" \
        #    -command {tk_messageBox -type ok -title "HELP" \
        #      -message "Draw dashed lines."}] \
        #  -row 7 -column 0 -sticky w
        #grid [label $wpo.dash_label -text "Dashed lines:"] \
        #  -row 7 -column 1 -sticky w
        #grid [frame $wpo.dash_frame] \
        #  -row 7 -column 2 -sticky w
        #tk_optionMenu $wpo.dash_frame.list ${ns}::dash "-"
        #$wpo.dash_frame.list.menu delete 0
        #foreach dsh "no - , . _" {
        #  $wpo.dash_frame.list.menu add radiobutton -label $dsh \
        #      -variable ${ns}::dash
        #}
        #pack $wpo.dash_frame.list -side left -anchor w -fill x  
        
        ${ns}::loadCoordinates
        if {$ndim == 3} {
          ${ns}::drawArrows
        }

        return $w
      }

    }

    #lappend namespaces $ns
    #lappend nmwizguis [string range $ns 2 end]
    return $ns
  }
}

proc nmwiz_tk {} {
  ::nmwiz::initGUI
  return $::nmwiz::w
}

proc nmwiz_load {filename} {
  nmwiz_tk
  ::nmwiz::loadNMD $filename
} 

#nmwiz_tk
