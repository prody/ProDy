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
package provide nmwiz 0.6

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

proc ::nmwiz_MultiPlot::init_plot {args} {
   variable parent
   variable verbose

   set  parent [lindex $args 0]
   incr ::nmwiz_MultiPlot::plotcount
   set ns "::nmwiz_MultiPlot::Plot${::nmwiz_MultiPlot::plotcount}"

   if {$verbose} {
     if {[namespace exists $ns]} {
       puts "Reinitializing namespace $ns."
     } else {
       puts "Creating namespace $ns"
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
         ${nmwizns}::Select_residue [::tcl::mathfunc::int [expr ($x-$xplotmin)/$scalex+$xmin]] $clr  
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
         pack  $w.info -side top -pady 2m -ipadx 5m -ipady 2m
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

namespace eval ::nmwiz:: {
  namespace export nmwizgui
  namespace export initialize
  
  variable license "NMWiz: Normal Mode Visualization, Animation, and Plotting 
University of Illinois Open Source License
Copyright 2010-2011 Ahmet Bakan,
All rights reserved.

Developed by:
Ahmet Bakan
http://www.csb.pitt.edu/People/abakan/

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the Software), to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution.
Neither the names of Theoretical and Computational Biophysics Group, University of Illinois at Urbana-Champaign, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.
THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR\ 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR\ 
OTHER DEALINGS WITH THE SOFTWARE.\
\n\n\
NMWiz makes use of a modified version of VMD plugin MultiPlot, which is also distributed under UIUC Open Source License.\
"
  
  variable guicount -1
  variable tmpdir
  variable titles [list]
  variable plothandles [list] 
  variable openfiles [list]
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
  
  proc init_gui {} {
    variable w
    variable platform  
    # If already initialized, just turn on
    if [winfo exists .nmwizgui] {
      wm deiconify .nmwizgui
      raise .nmwizgui
      return 
    }
    set w [toplevel .nmwizgui]
    wm title $w "NMWiz 0.6 - Main"
    wm resizable $w 0 0

    set wmf [frame $w.mainframe -bd 2]
    
    grid [button $wmf.loadnmd -width 20 -pady 2 -text "Load NMD file" -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"NMD files" { .nmd .NMD }} {"Text files" { .txt .TXT }} {"All files" *}}]
        if {![string equal $tempfile ""]} {::nmwiz::load_nmd $tempfile}}] \
      -row 5 -column 0 -columnspan 2

    grid [button $wmf.prody -width 20 -pady 2 -text "ProDy Interface" -command ::nmwiz::initProdyGUI] \
      -row 6 -column 0 -columnspan 2
   
    if {$platform != "windows"} {
      grid [button $wmf.retrieve -width 20 -pady 2 -text "ANM Server Interface" -command ::nmwiz::init_anm_interface] \
        -row 7 -column 0 -columnspan 2
    }

    grid [button $wmf.settings -width 8 -pady 2 -text "Settings" \
        -command ::nmwiz::initSettingsGUI] \
      -row 8 -column 0 -sticky we
    grid [button $wmf.website -width 8 -pady 2 -text "Website" \
        -command "vmd_open_url http://www.csb.pitt.edu/NMWiz/"] \
      -row 8 -column 1 -sticky we
    #grid [button $wmf.license -pady 2 -text "License" \
    #    -command {tk_messageBox -type ok -title "NMWiz License" \
    #      -message $::nmwiz::license} ] \
    #  -row 8 -column 2 -sticky we
    
    grid [frame $wmf.options] -row 10 -column 0 -columnspan 2
    
    
    grid [button $wmf.options.pcv_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "If checked, loading a new dataset will not change the current view."}] \
      -row 0 -column 0 -sticky w
    grid [label $wmf.options.pcv_label -text "Preserve current view:"] \
      -row 0 -column 1
    grid [checkbutton $wmf.options.preserview -text "" \
        -variable ::nmwiz::preserview] \
      -row 0 -column 2

    #pack $wmf.options -side top -fill x -expand 1
    pack $wmf -side top -fill x -expand 1
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
    }
    default {
      variable pybin "[::ExecTool::find python]"
      variable pyANM "[::ExecTool::find anm.py]"
      variable pyPCA "[::ExecTool::find pca.py]"
    }
  }
  variable outputdir [pwd]
  variable defaultColor "purple"
  variable settings [dict create anm $pyANM pca $pyPCA color $defaultColor outputdir $outputdir pybin $pybin]
  proc saveSettings {} {
    puts "Saving NMWiz settings"
    dict set ::nmwiz::settings anm $::nmwiz::pyANM
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
    puts "Loading NMWiz settings"
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
      variable ::nmwiz::pyPCA "[dict get $::nmwiz::settings pca]"
      variable ::nmwiz::defaultColor [dict get $::nmwiz::settings color]
      variable ::nmwiz::outputdir "[dict get $::nmwiz::settings outputdir]"
      variable ::nmwiz::pybin "[dict get $::nmwiz::settings pybin]"
    }
  }
  loadSettings

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

    grid [button $wf.pyHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Path to the Python executable. If the folder that contains \"python\" or \"python.exe\"\
is included in your environment variable PATH, you may keep this as \"python\", otherwise\
specify the path to the executable, e.g. \"C:\\python27\\python.exe\""}] \
      -row 1 -column 0 -sticky w
    grid [label $wf.pyLabel -text "Python:"] \
      -row 1 -column 1 -sticky w
    grid [entry $wf.pyEntry -width 20 -textvariable ::nmwiz::pybin] \
      -row 1 -column 2 -sticky ew
    grid [button $wf.pyBrowse -text "Browse" -pady 2 \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"All files" *}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::python $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 1 -column 3 -sticky ew
      
    grid [button $wf.anmHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Full path to ProDy ANM script (anm.py), e.g. C:\\python27\\Scripts\\anm.py"}] \
      -row 3 -column 0 -sticky w
    grid [label $wf.anmLabel -text "ANM script:"] \
      -row 3 -column 1 -sticky w
    grid [entry $wf.anmEntry -width 20 -textvariable ::nmwiz::pyANM] \
      -row 3 -column 2 -sticky ew
    grid [button $wf.anmBrowse -text "Browse" -pady 2 \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"ANM Script" { anm.py }}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::pyANM $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 3 -column 3 -sticky ew

    grid [button $wf.pcaHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Full path to the ProDy PCA script (pca.py) e.g. C:\\python27\\Scripts\\pca.py"}] \
      -row 5 -column 0 -sticky w
    grid [label $wf.pcaLabel -text "PCA script:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wf.pcaEntry -width 20 -textvariable ::nmwiz::pyPCA] \
      -row 5 -column 2 -sticky ew
    grid [button $wf.pcaBrowse -text "Browse" -pady 2 \
        -command {
      set tempfile [tk_getOpenFile \
        -filetypes {{"PCA Script" { pca.py }}}]
        if {![string equal $tempfile ""]} {set ::nmwiz::pyPCA $tempfile}
        ::nmwiz::saveSettings
        }] \
      -row 5 -column 3 -sticky ew

    grid [button $wf.clrHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "The default color for arrow graphics."}] \
      -row 7 -column 0 -sticky w
    grid [label $wf.scriptLabel -text "Default color:"] \
      -row 7 -column 1 -sticky w
    grid [frame $wf.colorFrame] \
      -row 7 -column 2 -sticky ew
    tk_optionMenu $wf.colorFrame.list ::nmwiz::defaultColor "" 
    $wf.colorFrame.list.menu delete 0 last
    foreach acolor $::nmwiz::nmwizColors {
      $wf.colorFrame.list.menu add radiobutton -label $acolor \
          -variable ::nmwiz::defaultColor \
          -command ::nmwiz::saveSettings
    }
    pack $wf.colorFrame.list -side left -anchor w -fill x

    grid [button $wf.prodySubmit -text "Save and Close" -pady 2 \
        -command "::nmwiz::saveSettings; destroy .nmwizsettings"] \
      -row 15 -column 0 -columnspan 4 -sticky we

    pack $wf -side top -fill x -expand 1
  }
  
  variable prodyMolecule
  variable prodyMolid -1
  variable prodyNFrames 0
  variable prodySelstr "name CA"
  variable prodySelAtoms 0
  variable prodyScript "ANM"
  variable prodyPrefix ""
  variable prodyAllfig 0
  variable prodyAllnum 0
  variable prodyTask ""
  variable prodyFrame 0
  variable prodyCutoff 15
  variable prodyGamma 1
  variable prodyNModes 10
  variable prodyFirstFrame 0
  variable prodySkipFrame 0 
  variable prodyLastFrame end
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
    grid [button $wmf.molHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Select a molecule to be used in calculations."}] \
      -row 2 -column 0 -sticky w
    grid [label $wmf.molLabel -text "Molecule:"] \
      -row 2 -column 1 -sticky w
    grid [frame $wmf.molFrame] \
      -row 2 -column 2 -sticky ew
    tk_optionMenu $wmf.molFrame.list ::nmwiz::prodyMolecule "" 
    grid [button $wmf.molUpdate -text "Update" -pady 2 \
        -command ::nmwiz::prodyUpdateMolList] \
      -row 2 -column 3 -sticky ew
    
    
    grid [button $wmf.molinfoHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Shows number of atoms and frames for the selected molecule."}] \
      -row 3 -column 0 -sticky w
    grid [label $wmf.molinfoLbl -text "Information:"] \
      -row 3 -column 1 -sticky w
    grid [label $wmf.molinfoLabel -text ""] \
      -row 3 -column 2 -columnspan 2 -sticky w
    
    grid [button $wmf.selstr -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Atom selection string. Click Select to update your \
selection."}] \
      -row 5 -column 0 -sticky w
    grid [label $wmf.selstrLabel -text "Selection:"] \
      -row 5 -column 1 -sticky w
    grid [entry $wmf.selstrEntry -width 20 -textvariable ::nmwiz::prodySelstr] \
      -row 5 -column 2 -sticky ew
    grid [button $wmf.selUpdate -text "Select" -pady 2 \
        -command ::nmwiz::prodyUpdateSelection] \
      -row 5 -column 3 -sticky ew
      
    grid [button $wmf.selinfoHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Shows number of selected atoms."}] \
      -row 6 -column 0 -sticky w
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
    foreach script "ANM PCA" {
      $wf.scriptFrame.list.menu add radiobutton -label "$script calculation" \
          -variable ::nmwiz::prodyTask \
          -command "set ::nmwiz::prodyScript $script; ::nmwiz::prodyChangeTask; ::nmwiz::prodyUpdatePrefix"
      incr counter  
    }
    pack $wf.scriptFrame.list -side left -anchor w -fill x
    variable prodyTask "ANM calculation"
    grid [button $wf.nmwizSettings -text "Settings" -pady 2 \
        -command "::nmwiz::initSettingsGUI"] \
      -row 7 -column 3 -sticky ew

    grid [button $wf.outdHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Select a directory for writing output files."}] \
      -row 8 -column 0 -sticky w
    grid [label $wf.outdLabel -text "Output dir:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.outdEntry -width 20 -textvariable ::nmwiz::outputdir] \
      -row 8 -column 2 -sticky w
    grid [button $wf.outdBrowse -text "Browse" -pady 2 \
        -command {
      set tempdir [tk_chooseDirectory -initialdir $::nmwiz::outputdir ]
        if {![string equal $tempdir ""]} {set ::nmwiz::outputdir $tempdir}
        ::nmwiz::saveSettings
        }] \
      -row 8 -column 3 -sticky ew

    grid [button $wf.filepHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Select the prefix for output files."}] \
      -row 9 -column 0 -sticky w
    grid [label $wf.filepLabel -text "File prefix:"] \
      -row 9 -column 1 -sticky w
    grid [entry $wf.filepEntry -width 20 -textvariable ::nmwiz::prodyPrefix] \
      -row 9 -column 2 -sticky w


    #grid [button $wf.numoutHelp -text "?" -padx 0 -pady 0 \
    #    -command {tk_messageBox -type ok -title "HELP" \
    #      -message "Check this if you want all default numerical output files to be written."}] \
    #  -row 10 -column 0 -sticky w
    #grid [label $wf.numoutLabel -text "Save all numerical output files:"] \
    #  -row 10 -column 1 -columnspan 2 -sticky w
    #grid [checkbutton $wf.numoutCheck -text "" \
    #    -variable ::nmwiz::prodyAllnum] \
    #  -row 10 -column 3 -sticky w


    #grid [button $wf.figoutHelp -text "?" -padx 0 -pady 0 \
    #    -command {tk_messageBox -type ok -title "HELP" \
    #      -message "Check this if you want all default graphical output files to be written."}] \
    #  -row 12 -column 0 -sticky w
    #grid [label $wf.figoutLabel -text "Save all graphical output files:"] \
    #  -row 12 -column 1 -columnspan 2 -sticky w
    #grid [checkbutton $wf.figoutCheck -text "" \
    #    -variable ::nmwiz::prodyAllfig] \
    #  -row 12 -column 3 -sticky w

    pack $wf -side top -fill x -expand 1
    

    # ANM frame
    set wf [labelframe $prodyGUI.anmFrame -text "ANM Settings" -bd 2]
    grid [button $wf.frameHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter index of the frame for the selected molecule to be used\
in the calculations. Index of the very first frame is 0,"}] \
      -row 6 -column 0 -sticky w
    grid [label $wf.frameLabel -text "Frame number:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.frameEntry -width 4 -textvariable ::nmwiz::prodyFrame] \
      -row 6 -column 2 -sticky w

    grid [button $wf.cutoffHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter cutoff distance for interactions between selected atoms."}] \
      -row 8 -column 0 -sticky w
    grid [label $wf.cutoffLabel -text "Cutoff distance:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.cutoffEntry -width 4 -textvariable ::nmwiz::prodyCutoff] \
      -row 8 -column 2 -sticky w

    grid [button $wf.gammaHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the force constant value."}] \
      -row 10 -column 0 -sticky w
    grid [label $wf.gammaLabel -text "Force constant:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wf.gammaEntry -width 4 -textvariable ::nmwiz::prodyGamma] \
      -row 10 -column 2 -sticky w    
    pack $wf -side top -fill x -expand 1
    
    grid [button $wf.modesHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the number of non-zero eigenvalues/vectors to calculate."}] \
      -row 12 -column 0 -sticky w
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 12 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::nmwiz::prodyNModes] \
      -row 12 -column 2 -sticky w   
    pack $wf -side top -fill x -expand 1
    
    # PCA frame
    set wf [labelframe $prodyGUI.pcaFrame -text "PCA Settings" -bd 2]
    
    grid [button $wf.firstHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the index of the first frame to be used in calculations.\
Index of the very first frame is 0."}] \
      -row 6 -column 0 -sticky w
    grid [label $wf.firstLabel -text "First frame:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wf.firstEntry -width 4 -textvariable ::nmwiz::prodyFirstFrame] \
      -row 6 -column 2 -sticky w
    
    grid [button $wf.skipHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the number of frames to skip after each frame."}] \
      -row 8 -column 0 -sticky w
    grid [label $wf.skipLabel -text "Skip frame:"] \
      -row 8 -column 1 -sticky w
    grid [entry $wf.skipEntry -width 4 -textvariable ::nmwiz::prodySkipFrame] \
      -row 8 -column 2 -sticky w

    grid [button $wf.lastHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the index of the last frame number to be used in calculations.\
\"end\" can be used as the index of the very last frame."}] \
      -row 10 -column 0 -sticky w
    grid [label $wf.lastLabel -text "Last frame:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wf.lastEntry -width 4 -textvariable ::nmwiz::prodyLastFrame] \
      -row 10 -column 2 -sticky w


    grid [button $wf.modesHelp -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter the number of non-zero eigenvalues/vectors to calculate."}] \
      -row 12 -column 0 -sticky w
    grid [label $wf.modesLabel -text "Number of modes:"] \
      -row 12 -column 1 -sticky w
    grid [entry $wf.modesEntry -width 4 -textvariable ::nmwiz::prodyNModes] \
      -row 12 -column 2 -sticky w
      
      
    # Submit button
    set wf [frame $prodyGUI.submitFrame -bd 2]
      
    grid [button $wf.prodySubmit -text "Submit ProDy Job" -pady 2 \
        -command ::nmwiz::prodySubmitJob] \
      -row 0 -column 0 -sticky we
    grid [button $wf.prodyWebsite -text "Go to ProDy Website" -pady 2 \
        -command "vmd_open_url http://www.csb.pitt.edu/ProDy"] \
      -row 0 -column 1 -sticky we
    pack $wf -side top -fill x -expand 1
         
  }
  
  proc prodyUpdateMolList {} {
    variable prodyGUI
    set wf $prodyGUI.mainFrame
    $wf.molFrame.list.menu delete 0 last
    set counter 0
    foreach id [molinfo list] {
      if {[molinfo $id get numatoms] > 0 && [molinfo $id get numframes] > 0} {
        if {$counter == 0} {
          variable prodyMolid $id
        }
        $wf.molFrame.list.menu add radiobutton -label "[molinfo $id get name]" \
            -variable ::nmwiz::prodyMolecule \
            -command "set ::nmwiz::prodyMolid $id; ::nmwiz::prodyUpdateMolinfo"
        incr counter  
      }
    }
    pack $wf.molFrame.list -side left -anchor w -fill x
    variable prodyMolid
    if {$prodyMolid > -1} {
      variable prodyMolecule "[molinfo $prodyMolid get name]"
      ::nmwiz::prodyUpdateMolinfo
    }
  }
  
  proc prodyCheckMolecule {} {
    if {[lsearch [molinfo list] $::nmwiz::prodyMolid] == -1} {
      ::nmwiz::prodyUpdateMolList
      return 0
    } 
    return 1
  }
  
  proc prodyUpdateMolinfo {} {
    if {$::nmwiz::prodyMolid > -1} {
      ::nmwiz::prodyCheckMolecule
      set ::nmwiz::prodyNFrames [molinfo $::nmwiz::prodyMolid get numframes]
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "[molinfo $::nmwiz::prodyMolid get numatoms] atoms, $::nmwiz::prodyNFrames frames"
      ::nmwiz::prodyUpdateSelection
      ::nmwiz::prodyUpdatePrefix
    } else {
      set ::nmwiz::prodyNFrames 0
      .nmwizprody.mainFrame.molinfoLabel configure \
        -text "Load a molecule and click Update."
    }
  }
  
  proc prodyUpdatePrefix {} {
    if {[::nmwiz::prodyCheckMolecule]} {
      set prefix [molinfo $::nmwiz::prodyMolid get name]
      if {[string range $prefix [expr [string length $prefix] - 4] end] == ".pdb"} {
        set prefix [string range $prefix 0 [expr [string length $prefix] - 5]]
      }
      if {$::nmwiz::prodyScript == "ANM"} {
        set ::nmwiz::prodyPrefix "$prefix\_anm"
      } else {
        set ::nmwiz::prodyPrefix "$prefix\_pca"
      }
    }
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
  
  proc prodyChangeTask {} {
    variable prodyGUI
    variable prodyScript
    if {$prodyScript == "ANM"} {
      pack forget $prodyGUI.pcaFrame
      pack forget $prodyGUI.submitFrame
      pack $prodyGUI.anmFrame -side top -fill x -expand 1
      pack $prodyGUI.submitFrame -side top -fill x -expand 1
    } else {
      pack forget $prodyGUI.anmFrame
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
    if {$::nmwiz::prodyScript == "ANM"} {
      ::nmwiz::prodySubmitANMjob
    } else {
      ::nmwiz::prodySubmitPCAjob
    }
    if {$::nmwiz::pybin == "" || $::nmwiz::pybin == {}} {
      variable ::nmwiz::pyPCA "[::ExecTool::find -interactive python.exe]"
      ::nmwiz::saveSettings
    }
  }
  proc prodySubmitANMjob {} {
    if {$::nmwiz::pyANM == "" || $::nmwiz::pyANM == {} || $::nmwiz::pyANM == "{}"} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Please specify the path to the ANM Script (anm.py) script."
      ::nmwiz::initSettingsGUI
      return
    }
    if {!([string is digit $::nmwiz::prodyFrame] && $::nmwiz::prodyFrame >= 0 && 
        $::nmwiz::prodyFrame < [molinfo $::nmwiz::prodyMolid get numframes])} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Frame number must be a number and must be in the valid range."
      return 
    }
    if {!([string is double $::nmwiz::prodyCutoff] && $::nmwiz::prodyCutoff > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Cutoff distance must be a number and must be greater than 0."
      return 
    }
    if {!([string is double $::nmwiz::prodyGamma] && $::nmwiz::prodyGamma > 0)} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Force constant must be a number and must be greater than 0."
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
    puts "Executing: $::nmwiz::pybin $::nmwiz::pyANM --quiet -s \"all\" -o \"$::nmwiz::outputdir\" -p \"$prefix\" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyCutoff -g $::nmwiz::prodyGamma \"$pdbfn\""
    set status [exec $::nmwiz::pybin $::nmwiz::pyANM --quiet -s all -o "$::nmwiz::outputdir" -p "$prefix" -n $::nmwiz::prodyNModes -c $::nmwiz::prodyCutoff -g $::nmwiz::prodyGamma "$pdbfn"]

    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy ANM calculation is finished and results are being loaded."
      ::nmwiz::load_nmd "$prefix.nmd" 
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
    }
  }  
  proc prodySubmitPCAjob {} {
    if {$::nmwiz::pyPCA == "" || $::nmwiz::pyPCA == {} || $::nmwiz::pyPCA == "{}"} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Please specify the path to the PCA Script (pca.py) script."
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
    set pdbfn [file join $::nmwiz::outputdir $::nmwiz::prodyPrefix.pdb]
    set end $::nmwiz::prodyLastFrame 
    if {$end == "end"} {
      set end [expr $::nmwiz::prodyNFrames -1]
    }
    #puts "animate write pdb $pdbfn beg $::nmwiz::prodyFirstFrame end $end skip $::nmwiz::prodySkipFrame waitfor all sel $sel $::nmwiz::prodyMolid"
    set nwritten [animate write pdb $pdbfn beg $::nmwiz::prodyFirstFrame end $end skip $::nmwiz::prodySkipFrame waitfor all sel $sel $::nmwiz::prodyMolid]
    puts "$nwritten frames are written as $pdbfn"
    
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
    puts "Executing: $::nmwiz::pybin $::nmwiz::pyPCA --quiet -s all -o \"$::nmwiz::outputdir\" -p \"$prefix\" -n $::nmwiz::prodyNModes \"$pdbfn\""
    set status [exec $::nmwiz::pybin $::nmwiz::pyPCA --quiet -s all -o "$::nmwiz::outputdir" -p "$prefix" -n $::nmwiz::prodyNModes "$pdbfn"]
    if {$status != -1} {
      tk_messageBox -type ok -title "INFO" \
        -message "ProDy PCA calculation is finished and results are being loaded."
      ::nmwiz::load_nmd "$prefix.nmd" 
    }  else {
      tk_messageBox -type ok -title "ERROR" \
        -message "An error occured."
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
    
    grid [button $wmf.id_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter 4 character PDB identifier."}] \
      -row 2 -column 0 -sticky w
    grid [label $wmf.id_label -text "PDB identifier:"] \
      -row 2 -column 1 -sticky w
    grid [entry $wmf.id_entry -width 4 -textvariable ::nmwiz::anm_id] \
      -row 2 -column 2 -sticky w
      
    grid [button $wmf.coordmod_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Choose whether to use PDB coordinates or biological unit coordinates."}] \
      -row 3 -column 0 -sticky w
    grid [radiobutton $wmf.coordmod_pdb -text "pdb coordinates" \
            -variable ::nmwiz::anm_coormod -value "pdb"] \
      -row 3 -column 1 -sticky w
    grid [radiobutton $wmf.coordmod_bio -text "biological unit" \
            -variable ::nmwiz::anm_coormod -value "bio"] \
      -row 3 -column 2 -sticky w
    
    grid [button $wmf.chain_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter 1 character chain identifier (default: all polypeptide chains)."}] \
      -row 4 -column 0 -sticky w
    grid [label $wmf.chain_label -text "Chain identifier:"] \
      -row 4 -column 1 -sticky w
    grid [entry $wmf.chain_entry -width 4 -textvariable ::nmwiz::anm_chain] \
      -row 4 -column 2 -sticky w

    grid [button $wmf.model_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter model (for multi-model files such as NMR structures)."}] \
      -row 6 -column 0 -sticky w
    grid [label $wmf.model_label -text "Model number:"] \
      -row 6 -column 1 -sticky w
    grid [entry $wmf.model_entry -width 4 -textvariable ::nmwiz::anm_model] \
      -row 6 -column 2 -sticky w

    grid [button $wmf.cutoff_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter cutoff for interaction between C alpha atoms (Å)."}] \
      -row 8 -column 0 -sticky w
    grid [label $wmf.cutoff_label -text "Cutoff distance (Å):"] \
      -row 8 -column 1 -sticky w
    grid [entry $wmf.cutoff_entry -width 4 -textvariable ::nmwiz::anm_cutoff] \
      -row 8 -column 2 -sticky w

    grid [button $wmf.pwr_help -text "?" -padx 0 -pady 0 \
        -command {tk_messageBox -type ok -title "HELP" \
          -message "Enter distance weight for interaction between C alpha atoms."}] \
      -row 10 -column 0 -sticky w
    grid [label $wmf.pwr_label -text "Distance weight:"] \
      -row 10 -column 1 -sticky w
    grid [entry $wmf.pwr_entry -width 4 -textvariable ::nmwiz::anm_pwr] \
      -row 10 -column 2 -sticky w
      
    grid [button $wmf.retrieve -width 24 -pady 2 -text "Submit to ANM Server" \
      -command ::nmwiz::submit_anm_server] \
      -row 16 -column 0 -columnspan 3
      
    grid [button $wmf.website -width 24 -pady 2 -text "Go to ANM Server" \
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
  
  proc load_nmd {fn} {
    variable filename $fn
    puts "NMWiz: Parsing file $filename"
    # Parse the file, and make sure coordinates are stored in the file
    variable openfiles
    #variable namespaces
    #variable nmwizguis
    
    if {[lsearch $openfiles $filename] > -1} {
      tk_messageBox -type ok -title "WARNING" \
        -message "Content of this file is already loaded in NMWiz."
      return
    }
    
    set nmdfile [open $filename]
    variable nmdlist [list]
    set coordinates 0 
    while {[gets $nmdfile nmdline] != -1} { 
      if {[lindex $nmdline 0] == "coordinates"} {
        if {[expr [llength $nmdline] % 3] != 1} {
          tk_messageBox -type ok -title "ERROR" \
            -message "Length of the coordinate array in $filename must be a multiple of 3. An array of length [llength $coordinates] is provided."
          return
        }
        set coordinates 1
      }
      lappend nmdlist $nmdline
    }
    if {$coordinates == 0} {
      tk_messageBox -type ok -title "ERROR" \
        -message "Coordinates were not found in the input file."
      return
    }
    variable guicount
    incr guicount
    set ns "::nmgui$guicount"
    
    namespace eval $ns {
      puts "DEBUG: evaluating namespace"
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
      variable resolution 10
      
      variable selstr "all"
      variable selrep 0
      
      variable autoupdate 1
      variable autoanimate 0
      
      variable hide_shorter 0.0 
      variable cylinder_radius 0.4
      variable cone_radius 0.6
      variable cone_height 1.0
      
      variable showproteinas "Backbone"
      variable tuberadius 0.4
      variable bondradius 0.3   
      variable spherescale 0.6
      variable cutoffdistance 8.0
      variable proteincolor "Mobility" 
      
      variable nframes 50
      
      variable betalist {}
      variable betamin
      variable betamax
      
      variable keepfiles 0
      variable overanim 1
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
      
      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      proc Close_mols {} {
        variable molid
        variable arrid
        variable animid
        variable selid
        if {[lsearch [molinfo list] $molid] > -1} {
          mol delete $molid
        }
        if {[lsearch [molinfo list] $arrid] > -1} {
          mol delete $arrid
        }
        if {[lsearch [molinfo list] $animid] > -1} {
          mol delete $animid
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
        puts "Plot handle: [lindex [subst $${ns}::plothandles] end]"
        #-dash [subst $${ns}::dash] \
      }
      
      proc Prepare_selmol {} {
        variable selid
        variable tempfn
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
          set outfile [open [file join $::nmwiz::tmpdir $tempfn] w]
          foreach line [[namespace current]::Get_pdb_lines $coordinates] {
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

      proc Clear_selection {} {
        variable selid
        if {$selid > -1 && [lsearch [molinfo list] $selid] > -1} {
          for {set i [molinfo $selid get numreps]} {$i >= 0} {incr i -1} {
            mol molrep $i $selid off
          }
        }
        label delete Atoms all
      }

      proc Select_residue {resid color} {

        variable plotrids
        set which [lsearch $plotrids $resid]
        if {$which == -1} {return 0}
        [namespace current]::Prepare_selmol
        
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
          puts "Deselected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
        } else {
          puts "Selected [lindex $chainids $which]:[lindex $resnames $which][lindex $resids $which]"
          mol showrep $selid $which on
          mol modstyle $which $selid VDW $selectscale $resolution
          mol modmaterial $which $selid $material
          mol modselect $which $selid "residue $which"
          mol modcolor $which $selid ColorID [lsearch "blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3" $color]
        }
      }


      proc Protein_representation {targetid} {
        variable molid

        if {[lsearch [molinfo list] $targetid] == -1} {
          if {$targetid == $molid} {
            [namespace current]::Load_coordinates
          } else {
            return 0
          }
        }
        
        
        variable showproteinas
        variable tuberadius
        variable bondradius
        variable cutoffdistance
        variable resolution
        variable material
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
          mol modstyle 0 $targetid VDW $spherescale $resolution
          mol modmaterial 0 $targetid $material
          mol addrep $targetid
          mol modstyle 1 $targetid DynamicBonds $cutoffdistance $bondradius $resolution
          mol modmaterial 1 $targetid $material
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
          mol modstyle 0 $targetid Tube $tuberadius $resolution
          mol modmaterial 0 $targetid $material
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

      proc Beta_msf {} {
        variable molid
        variable activemode
        variable indices
        variable lengths
        variable scalearrows
        variable modes
        variable animid
        variable material
        variable selstr

        set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [lindex $modes [lsearch $indices $activemode]]
        set mode [vecscale [expr $length * $length] [vecmul $mode $mode]]

        set index 0
        variable betalist {}
        variable betamin 10000
        variable betamax -10000
        foreach {mx my mz} $mode {
          set beta [expr $mx + $my + $mz]
          lappend betalist $beta
          if {$beta < $betamin} {set betamin $beta}
          if {$beta > $betamax} {set betamax $beta}
          incr index
        }
        [atomselect $molid "all"] set beta $betalist 
        #mol scaleminmax $molid 0 [::tcl::mathfunc::min $betalist] [::tcl::mathfunc::max $betalist]
        color scale midpoint 0.1
        color scale method BWR
        
      }

      proc Draw_action {} {
        variable overwrite
        variable arrid
        variable arrids
        if {!$overwrite} {
          set arrid [mol new]
          lappend arrids $arrid
        }
        [namespace current]::Draw_arrows
      }
      
      proc Auto_update {} {
        variable autoupdate
        if {$autoupdate} {
          variable overwrite
          variable arrid
          variable arrids
          [namespace current]::Draw_arrows
        }
      }

      proc Draw_arrows {} {
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
          [namespace current]::Load_coordinates
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
        set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [vecscale [expr $length * $scalearrows] [lindex $modes [lsearch $indices $activemode]]]
        
        foreach index [[atomselect $molid $selstr] get index] {
          set from [expr $index * 3]
          set to  [expr $from + 2]
          set xyz [lrange $coordinates $from $to ] 
          set v [lrange $mode $from $to ] 
          #puts "$xyz $v"
          #set xyz "$x $y $z"
          #set v "$mx $my $mz"
          if {$hide_shorter < [veclength $v]} {
            set temp [vecadd $xyz $v]
            graphics $arrid cylinder $xyz $temp radius $cylinder_radius resolution $resolution
            set temp2 [vecadd $temp [vecscale $v [expr $cone_height / [veclength $v]]]]
            graphics $arrid cone $temp $temp2 radius $cone_radius resolution $resolution
          }
        }
        mol rename $arrid "$title mode $activemode arrows"
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]
        display resetview
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
        [namespace current]::Beta_msf
        $w.draw_arrows.arrowbuttons_label configure -text "Arrows ($arrid):"
        
        variable arridlist
        lset arridlist $whichmode $arrid
      }

      proc Get_pdb_lines {coords} {
        variable atomnames
        variable resnames
        variable chainids
        variable resids
        variable betalist
        set pdblines ""
        set i 0
        foreach an $atomnames rn $resnames ci $chainids ri $resids {x y z} $coords b $betalist {
          incr i
          lappend pdblines [format "ATOM  %5d  %-3s %-4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f" \
                               $i  $an $rn  $ci $ri $x $y $z 1.0 $b]
        
        }
        return $pdblines
        
      }

      proc Locate_coordinates {} {
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
        [namespace current]::Load_coordinates
      }


      proc Animate {} {
        variable nframes
        variable keepfiles
        variable activemode
        variable coordinates
        variable prefix
        variable lengths
        variable indices
        variable scalearrows
        variable modes
        variable overanim
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
        if {$keepfiles} {
          set animfn "$prefix\_mode_$activemode\_animation.pdb"
        } else {
          set animfn [file join $::nmwiz::tmpdir $tempfn]
        }
        
        set outfile [open $animfn w]
        
        set length [lindex $lengths [lsearch $indices $activemode]]
        set mode [vecscale [expr $length * [::tcl::mathfunc::abs $scalearrows]] [lindex $modes [lsearch $indices $activemode]]]
        
        # Write the initial coordinates as model 0
        # so that VMD sets bonds correctly
        
        puts $outfile "MODEL      0"
        foreach line [[namespace current]::Get_pdb_lines $coordinates] {
          puts $outfile $line
        } 
        puts $outfile "ENDMDL"
        
        set coords [vecadd $coordinates $mode]
        set mode [vecscale $mode [expr  -2.0 / $nframes]]
        for {set i 0} {$i <= $nframes} {incr i} {
          puts $outfile [format "MODEL %6d" [expr $i + 1]]
          foreach line [[namespace current]::Get_pdb_lines $coords] {
            puts $outfile $line
          } 
          puts $outfile "ENDMDL"
          set coords [vecadd $coords $mode]
        }
        
        close $outfile
        set currentview [molinfo $molid get {rotate_matrix center_matrix scale_matrix global_matrix}]

        if {[lsearch [molinfo list] $animid] == -1 || !$overanim} {
          set animid [mol new $animfn]
        } else {
          animate delete beg 0 end -1 skip 0 $animid
          mol addfile $animfn waitfor all $animid
        }
        
        animate delete beg 0 end 0 skip 0 $animid

        mol top $animid
        #mol modstyle 0 $animid Tube 0.3
        #mol modcolor 0 $animid Beta
        #mol scaleminmax $animid 0 $betamin $betamax
        [namespace current]::Protein_representation $animid 
        $w.draw_arrows.animbuttons_label configure -text "Animation ($animid):"
        mol rename $animid "$title mode $activemode animation"
        foreach id [molinfo list] {
          molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
        }
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


      proc Load_coordinates {} {
        variable molid
        variable coordinates
        variable title
        variable w
        variable tempfn
        if {[lsearch [molinfo list] $molid] != -1} {
          return 0
        }
        
        set outfile [open [file join $::nmwiz::tmpdir $tempfn] w]
        foreach line [[namespace current]::Get_pdb_lines $coordinates] {
          puts $outfile $line
        } 
        close $outfile
        if {[molinfo num] > 0 && $::nmwiz::preserview} {
          set currentview [molinfo [lindex [molinfo list] 0] get {rotate_matrix center_matrix scale_matrix global_matrix}]
        }

        set molid [mol new [file join $::nmwiz::tmpdir $tempfn]]

        $w.draw_arrows.protbuttons_label configure -text "Protein ($molid):"
        mol rename $molid "$title coordinates"
        [namespace current]::Beta_msf
        [namespace current]::Protein_representation $molid

        if {[molinfo num] > 0 && $::nmwiz::preserview} {
          foreach id [molinfo list] {
            molinfo $id set {rotate_matrix center_matrix scale_matrix global_matrix} $currentview
          }
        }
      }

      proc Check_coordinates {} {
        variable molid
        variable modes
        variable coordinates
        if {[[atomselect $molid "all"] num] != [llength [lindex $modes 0]]} {
          mol delete $molid
          set molid -1
          if {"ok" == [tk_messageBox -type okcancel -title "ERROR" \
              -message "[[atomselect $molid all] num] atoms are loaded. Coordinate data file must contain [llength [lindex $modes 0]] atoms. Please locate the correct file."]} {
            [namespace current]::Locate_coordinates
          } 
        } else {
          set coordinates [[atomselect $molid all] get {x y z}]
        }
      }
      
      proc Change_color {} {
        variable color
        variable colorlist
        variable indices
        variable activemode
        lset colorlist [lsearch $indices $activemode] $color
      }
      
      proc Change_mode {} {
        variable w
        variable activemode
        variable indices
        variable arrid
        variable arridlist
        variable animid
        variable animidlist
        variable color
        variable colorlist
        
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
          [namespace current]::Beta_msf
        } else {
          [namespace current]::Draw_action
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
      }
      
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
      proc initialize {} {
        set nmdlist $::nmwiz::nmdlist
        set filename $::nmwiz::filename
        variable arridlist
        variable animidlist
        variable colorlist
        variable coordinates ""
        variable atomnames ""
        variable resnames ""
        variable chainids ""
        variable resids [list]
        variable plotrids [list]
        variable bfactors [list]
        variable bfactormin
        variable bfactormax
        variable modes [list]
        variable lengths [list]
        variable indices [list]
        variable title ""
        variable n_dims 0
        variable n_atoms
        
        foreach nmdline $nmdlist {
          switch -exact [lindex $nmdline 0] {
            coordinates {
              set coordinates [lrange $nmdline 1 end]
              set n_atoms [expr [llength $coordinates] / 3]
            }
          }
        }     

        
        # Evaluate each line
        foreach nmdline $nmdlist {
          switch -exact [lindex $nmdline 0] {
            name {
              variable title [lrange $nmdline 1 end]
              if {[lsearch $::nmwiz::titles $title] > -1 && [lsearch $::nmwiz::openfiles $filename] > -1} {
                set title "$title ($::nmwiz::guicount)"
              }
              lappend ::nmwiz::titles $title
              puts $::nmwiz::titles
            }
            coordinates {
            }
            atomnames {
              if {[llength $nmdline] != [expr $n_atoms + 1]} {
                puts "NMWiz WARNING: Length of atomnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
              } else {
                set atomnames [lrange $nmdline 1 end]
              }
            }
            resnames {
              if {[llength $nmdline] != [expr $n_atoms + 1]} {
                puts "NMWiz WARNING: Length of resnames array must be $n_atoms, not [expr [llength $nmdline] -1]."
              } else {
                set resnames [lrange $nmdline 1 end]
              }
            }
            chainids {
              if {[llength $nmdline] != [expr $n_atoms + 1]} {
                puts "NMWiz WARNING: Length of chainids array must be $n_atoms, not [expr [llength $nmdline] -1]."
              } else {
                set chainids [lrange $nmdline 1 end]
              }
            }
            resids {
              if {[llength $nmdline] != [expr $n_atoms + 1]} {
                puts "NMWiz WARNING: Length of resids array must be $n_atoms, not [expr [llength $nmdline] -1]."
              } else {
                set resids [lrange $nmdline 1 end]
              }
            }
            bfactors {
              if {[llength $nmdline] != [expr $n_atoms + 1]} {
                puts "NMWiz WARNING: Length of bfactors array must be $n_atoms, not [expr [llength $nmdline] -1]."
              } else {
                set bfactors [lrange $nmdline 1 end]
              }      
            }
            mode {
              set l [expr [llength $nmdline] - 1]
              if {$l >= [llength $coordinates]} {
                if {$n_dims == 0} {
                  set n_dims 3
                  puts "NMWiz INFO: File contains a 3D model."
                } elseif {$n_dims != 3} {
                  tk_messageBox -type ok -title "ERROR" \
                    -message "All modes must have the same dimensions."
                  return
                }
                switch -exact [expr $l - [llength $coordinates]] {
                  0 {
                    lappend modes [lrange $nmdline 1 end]
                    lappend indices [llength $modes]
                    lappend lengths 1
                    lappend arridlist -1
                    lappend animidlist -1
                    lappend colorlist $::nmwiz::defaultColor
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
                    lappend arridlist -1
                    lappend animidlist -1
                    lappend colorlist $::nmwiz::defaultColor
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
                    lappend arridlist -1
                    lappend animidlist -1
                    lappend colorlist $::nmwiz::defaultColor
                  } 
                  default {
                    puts "NMWiz WARNING: Mode data was not understood. Line starts with [lrange $nmdline 0 4]."
                  }
                }
              } else {
                if {$n_dims == 0} {
                  set n_dims 1
                  puts "NMWiz INFO: File contains a 1D model."
                } elseif {$n_dims != 1} {
                  tk_messageBox -type ok -title "ERROR" \
                    -message "All modes must have the same dimensions."
                  return
                }
                  tk_messageBox -type ok -title "ERROR" \
                    -message "1D models are not handled yet."
                  return
              } 
            }
            default {
              puts "NMWiz WARNING: Unrecognized line \"[lindex $nmdline 0]\""
            }
          }
          
        } 
        
        if {[llength $modes] == 0} {
          tk_messageBox -type ok -title "ERROR" \
            -message "Mode data was NOT found in the input file."
          return
        }
        
        
        if {$title == ""} {
          set title "Untitled ($::nmwiz::guicount)"
          puts "NMWiz INFO: Dataset is named as \"$title\"."
        }
        if {$atomnames == ""} {
          set atomnames [string repeat "CA " $n_atoms]
          puts "NMWiz INFO: Atoms are named as \"CA\"."
        }
        if {$resnames == ""} {
          set resnames [string repeat "GLY " $n_atoms]
          puts "NMWiz INFO: Residues are named as \"GLY\"."
        }
        if {$chainids == ""} {
          set chainids [string repeat "A " $n_atoms]
          puts "NMWiz INFO: Chains are named as \"A\"."
        }
        if {[llength $resids] == 0} {
          for {set i 1} {$i <= $n_atoms} {incr i} {lappend resids $i}
          puts "NMWiz INFO: Residues are numbered starting from 1."
        }
        foreach i [lrange $resids 0 end-1] j [lrange $resids 1 end] {
          if {$i >= $j} {
            for {set i 1} {$i <= $n_atoms} {incr i} {lappend plotrids $i}
            puts "NMWiz INFO: Residue numbers will NOT be used for plotting."
            break
          }
        }
        if {[llength $plotrids] == 0} {
          set plotrids $resids
          puts "NMWiz INFO: Residue numbers will be used for plotting."
        }
        if {[llength $bfactors] == 0} {
          puts "NMWiz INFO: Experimental bfactors were not found in the data file."
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
        

        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        variable drawlength [lindex $lengths 0]
        variable drawlengthstr [format "%.1f" [lindex $lengths 0]]
        variable activemode [lindex $indices 0]
        
        variable betalist
        foreach atnm $atomnames {
          lappend betalist 0
        } 
        variable scalearrows 
        if {[lindex $lengths 0] < 60.0} {
          set scalearrows [::tcl::mathfunc::int [expr 60.0 / [lindex $lengths 0]] ]
        }
      }
      initialize
      
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

        grid [frame $wam.active] \
          -row 0 -column 0 -columnspan 3
        label $wam.active.label -text "Active mode"

        button $wam.active.help -text "?" -padx 0 -pady 0 -command {
          tk_messageBox -type ok -title "HELP" \
            -message "Select the active mode for which you want to draw arrows\
or make an animation. The selected color effects both arrow graphics and plots."}
        tk_optionMenu $wam.active.list ${ns}::activemode 0
        $wam.active.list.menu delete 0
        variable indices
        variable lengths
        foreach index $indices length $lengths {
          $wam.active.list.menu add radiobutton -label $index \
              -variable ${ns}::activemode \
              -command "set ${ns}::drawlengthstr \[format \"%.1f\" \[lindex \$${ns}::lengths \[lsearch \$${ns}::indices \$${ns}::activemode]]]; ${ns}::Change_mode;"
        }
        
        tk_optionMenu $wam.active.color ${ns}::color "blue"
        $wam.active.color.menu delete 0
        foreach acolor "blue red gray orange yellow tan green white pink \
      cyan purple black yellow2 yellow3 green2 green3 \
      cyan2 cyan3 blue2 blue3 violet magenta magenta2 red2 red3 orange2 \
      orange3" {
          $wam.active.color.menu add radiobutton -label $acolor \
              -variable ${ns}::color \
              -command "${ns}::Change_color; ${ns}::Auto_update; "
        }
        pack $wam.active.help $wam.active.label $wam.active.list $wam.active.color -side left -anchor w -fill x
        
        #blue red gray orange yellow tan silver green white pink cyan purple lime mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet violet2 magenta magenta2 red2 red3 orange2 orange3
        
        grid [button $wam.scale_help -text "?" -padx 0 -pady 0 -command {
          tk_messageBox -type ok -title "HELP"\
            -message "Mode will be multiplied with the scalar value and the length\
of the mode (shown in the first box). If these modes are from a PCA, the\
length of the mode is the standard deviation along the principal mode. If the\
modes are from a NMA, the length of the mode is the square-root of the inverse\
of the eigenvalue corresponding to this mode."}] \
          -row 2 -column 0 -sticky w
        grid [label $wam.scale_label -text "Scale by:"] \
          -row 2 -column 1 -sticky w
        grid [frame $wam.scale_frame] \
          -row 2 -column 2 -sticky w
        entry $wam.scale_frame.length -width 5 \
          -textvariable ${ns}::drawlengthstr \
          -state disabled -disabledbackground white -disabledforeground black
        label $wam.scale_frame.angstrom -text "A"
        label $wam.scale_frame.product -text "x"
        button $wam.scale_frame.negate -text "+/-" -padx 0 -pady 0 -command "set ${ns}::scalearrows \[expr - \$${ns}::scalearrows]; ${ns}::Draw_action"  
        entry $wam.scale_frame.entry -width 4 -textvariable ${ns}::scalearrows
        button $wam.scale_frame.decr5 -text "-5" -padx 0 -pady 0 -command \
          "set ${ns}::scalearrows \[expr \$${ns}::scalearrows - 5]; ${ns}::Auto_update"
        button $wam.scale_frame.decr1 -text "-1" -padx 0 -pady 0 -command \
          "set ${ns}::scalearrows \[expr \$${ns}::scalearrows - 1]; ${ns}::Auto_update"
        button $wam.scale_frame.one -text "1" -padx 0 -pady 0 -command \
          "set ${ns}::scalearrows 1; ${ns}::Auto_update"
        button $wam.scale_frame.incr1 -text "+1" -padx 0 -pady 0 -command \
          "set ${ns}::scalearrows \[expr \$${ns}::scalearrows + 1]; ${ns}::Auto_update"
        button $wam.scale_frame.incr5 -text "+5" -padx 0 -pady 0 -command \
          "set ${ns}::scalearrows \[expr \$${ns}::scalearrows + 5]; ${ns}::Auto_update"
        pack $wam.scale_frame.length $wam.scale_frame.angstrom \
          $wam.scale_frame.product \
          $wam.scale_frame.negate \
          $wam.scale_frame.entry \
          $wam.scale_frame.decr5 $wam.scale_frame.decr1 \
          $wam.scale_frame.one \
          $wam.scale_frame.incr1 $wam.scale_frame.incr5 \
          -side left -anchor w -fill x
        pack $wam -side top -ipadx 10 -ipady 5 -fill x -expand 1

        grid [button $wam.selstr_help -text "?" -padx 0 -pady 0 -command {
          tk_messageBox -type ok -title "HELP" \
            -message "Atoms and arrows for them in the selection string will be displayed."}] \
          -row 3 -column 0 -sticky w
        grid [label $wam.selstr_label -text "Selection:"] \
          -row 3 -column 1 -sticky w
        grid [entry $wam.selstr_entry \
          -textvariable ${ns}::selstr] \
          -row 3 -column 2 -sticky we


        set wda [labelframe $w.draw_arrows -text "Actions" -bd 2]
        
        grid [label $wda.arrowbuttons_label -text "Arrows:"] \
          -row 5 -column 0 -sticky w
        grid [button $wda.ab_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Molecule id for the current arrow graphics is shown in parentheses.\n\nDraw : draw/redraw arrows for the active mode\nClean : remove most recently drawn arrows\nHide : hide/show most recently drawn arrows\nOptions : change arrow properties and drawing options"}] \
          -row 5 -column 1 -sticky w
        grid [button $wda.arrowbuttons_draw -width 4 -pady 1 -text "Draw" \
            -command ${ns}::Draw_arrows] \
          -row 5 -column 2
        grid [button $wda.arrowbuttons_clean -width 4 -pady 1 -text "Clean" \
            -command "foreach anarrid \$${ns}::arrids {if {\$anarrid != \$${ns}::arrid && \[lsearch \[molinfo list] \$${ns}::arrid] != -1} {mol delete \$anarrid}; if {\[lsearch \[molinfo list] \$${ns}::arrid] != -1} {graphics \$${ns}::arrid delete all}}"] \
          -row 5 -column 3
        grid [button $wda.arrowbuttons_showhide -width 4 -pady 1 -text "Hide" \
            -command "if {\[molinfo \$${ns}::arrid get displayed]} {mol off \$${ns}::arrid; \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Show} else {mol on \$${ns}::arrid; \$${ns}::w.draw_arrows.arrowbuttons_showhide configure -text Hide}"] \
          -row 5 -column 4
        grid [button $wda.arrowbuttons_options -width 5 -pady 1 -text "Options" \
            -command "if {\$${ns}::arropt} {pack forget \$${ns}::w.graphics_options; set ${ns}::arropt 0; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief raised} else {pack \$${ns}::w.graphics_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::arropt 1; \$${ns}::w.draw_arrows.arrowbuttons_options configure -relief sunken}"] \
          -row 5 -column 5

        grid [label $wda.animbuttons_label -text "Animation:"] \
          -row 6 -column 0 -sticky w
        grid [button $wda.anb_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Molecule id for the most recent animation is shown in parentheses.\n\nMake : animate fluctuations along the active mode\nPlay : play/pause the animation\nHide : hide/show the animation\nOptions : change animation options"}] \
          -row 6 -column 1 -sticky w
        grid [button $wda.animbuttons_animate -width 4 -pady 1 -text "Make" \
            -command ${ns}::Animate] \
          -row 6 -column 2
        grid [button $wda.animbuttons_stop -width 4 -pady 1 -text "Play" \
            -command "if {\$${ns}::animid == -1} {${ns}::Animate} else {if {\$${ns}::stopped} {mol top \$${ns}::animid; animate forward; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Pause; set ${ns}::stopped 0} else {animate pause; \$${ns}::w.draw_arrows.animbuttons_stop configure -text Play; set ${ns}::stopped 1}}"] \
          -row 6 -column 3
        grid [button $wda.animbuttons_showhide -width 4 -pady 1 -text "Hide" \
            -command "if {\$${ns}::animid > -1 && \[lsearch \[molinfo list] \$${ns}::animid] > -1} {if {\[molinfo \$${ns}::animid get displayed]} {animate pause; mol off \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Show} else {mol on \$${ns}::animid; \$${ns}::w.draw_arrows.animbuttons_showhide configure -text Hide; animate forward}}"] \
          -row 6 -column 4
        grid [button $wda.animbuttons_options -width 5 -pady 1 -text "Options" \
            -command "if {\$${ns}::anmopt} {pack forget \$${ns}::w.animation_options; set ${ns}::anmopt 0; \$${ns}::w.draw_arrows.animbuttons_options configure -relief raised} else {pack \$${ns}::w.animation_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::anmopt 1; \$${ns}::w.draw_arrows.animbuttons_options configure -relief sunken}"] \
          -row 6 -column 5

        grid [label $wda.plot_label -text "Plotting:"] \
          -row 8 -column 0 -sticky w
        grid [button $wda.plt_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Molecule id for displaying selected residues is shown in parentheses.\n\nPlot : plot squared-fluctuations along the active mode\nClear : clear all selections\nHide : hide/show the selected residues\nOptions : change plotting options"}] \
          -row 8 -column 1 -sticky w
        grid [button $wda.plot_plot -width 4 -pady 1 -text "Plot" \
            -command "${ns}::Plot ${ns}"] \
          -row 8 -column 2
        grid [button $wda.plot_clear -width 4 -pady 1 -text "Clear" \
            -command "${ns}::Clear_selection"] \
          -row 8 -column 3
        grid [button $wda.plot_showhide -width 4 -pady 1 -text "Hide" \
            -command "if {\$${ns}::selid > -1 && \[lsearch \[molinfo list] \$${ns}::selid] > -1} {if {\[molinfo \$${ns}::selid get displayed]} {mol off \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Show} else {mol on \$${ns}::selid; \$${ns}::w.draw_arrows.plot_showhide configure -text Hide}}"] \
          -row 8 -column 4
        grid [button $wda.plot_options -width 5 -pady 1 -text "Options" \
            -command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.plotting_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plotting_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \
          -row 8 -column 5
         
        ##-command "if {\$${ns}::pltopt} {pack forget \$${ns}::w.animation_options; set ${ns}::pltopt 0; \$${ns}::w.draw_arrows.plot_options configure -relief raised} else {pack \$${ns}::w.plot_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::pltopt 1; \$${ns}::w.draw_arrows.plot_options configure -relief sunken}"] \

        grid [label $wda.protbuttons_label -text "Protein:"] \
          -row 9 -column 0 -sticky w
        grid [button $wda.prt_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Molecule id of the molecular system is shown in parentheses.\n\nUpdate : Update protein representation\nFocus : reset view to focus on the molecular system\nHide : hide/show the molecular system\nOptions : change molecular system representation\n"}] \
          -row 9 -column 1 -sticky w
        grid [button $wda.prt_update -width 4 -pady 1 -text "Update" \
            -command "${ns}::Protein_representation \$${ns}::molid"] \
          -row 9 -column 2
        grid [button $wda.protbuttons_focus -width 4 -pady 1 -text "Focus" \
            -command "mol top \$${ns}::molid; display resetview"] \
          -row 9 -column 3  
        grid [button $wda.protbuttons_showhide -width 4 -pady 1 -text "Hide" \
            -command "if {\[molinfo \$${ns}::molid get displayed]} {mol off \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Show;} else {mol on \$${ns}::molid; \$${ns}::w.draw_arrows.protbuttons_showhide configure -text Hide;}"] \
          -row 9 -column 4
        grid [button $wda.protbuttons_repoptions -width 5 -pady 1 -text "Options" \
            -command "if {\$${ns}::prtopt} {pack forget \$${ns}::w.prograph_options; set ${ns}::prtopt 0; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief raised} else {pack \$${ns}::w.prograph_options -side top -ipadx 10 -ipady 5 -fill x -expand 1; set ${ns}::prtopt 1; \$${ns}::w.draw_arrows.protbuttons_repoptions configure -relief sunken}"] \
          -row 9 -column 5

        pack $wda -side top -fill x -expand 1

        set wgo [labelframe $w.graphics_options -text "Arrow Graphics Options" -bd 2]
        
        grid [button $wgo.auto_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, arrow graphics will be updated automatically when color selection, arrow scaling factor, arrow cone height, etc. variables change."}] \
          -row 0 -column 0 -sticky w
        grid [label $wgo.auto_label -text "Auto update:"] \
          -row 0 -column 1 -sticky w
        grid [checkbutton $wgo.auto_check -text "" \
            -variable ${ns}::autoupdate] \
          -row 0 -column 2 -sticky w

        grid [button $wgo.overwrite_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, when the active mode is changed, previously drawn mode will be hidden."}] \
          -row 1 -column 0 -sticky w
        grid [label $wgo.overwrite_label -text "Auto hide inactive mode:"] \
          -row 1 -column 1 -sticky w
        grid [checkbutton $wgo.overwrite_check -text "" \
            -variable ${ns}::overwrite] \
          -row 1 -column 2 -sticky w

        grid [button $wgo.hide_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Arrows shorter than the specified value will not be drawn."}] \
          -row 9 -column 0 -sticky w
        grid [label $wgo.hide_label -text "Draw if longer than:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wgo.hide_frame] \
          -row 9 -column 2 -sticky w
        entry $wgo.hide_frame.entry -width 4 -textvariable ${ns}::hide_shorter
        button $wgo.hide_frame.decr -text "-0.5" -padx 0 -pady 0 \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter - 0.5]]; ${ns}::Auto_update"
        button $wgo.hide_frame.incr -text "+0.5" -padx 0 -pady 0 \
          -command "set ${ns}::hide_shorter \[::tcl::mathfunc::abs \[expr \$${ns}::hide_shorter + 0.5]]; ${ns}::Auto_update"
        label $wgo.hide_frame.angstrom -text "A"
        pack $wgo.hide_frame.entry $wgo.hide_frame.decr $wgo.hide_frame.incr \
          $wgo.hide_frame.angstrom -side left -anchor w -fill x

        grid [button $wgo.cylinder_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Radius of the arrow cylinders."}] \
          -row 10 -column 0 -sticky w
        grid [label $wgo.cylinder_label -text "Arrow cylinder radius:"] \
          -row 10 -column 1 -sticky w
        grid [frame $wgo.cylinder_frame] \
          -row 10 -column 2 -sticky w
        entry $wgo.cylinder_frame.entry -width 4 -textvariable ${ns}::cylinder_radius
        button $wgo.cylinder_frame.decr -text "-0.1" -padx 0 -pady 0 \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius - 0.1]]; ${ns}::Auto_update"
        button $wgo.cylinder_frame.incr -text "+0.1" -padx 0 -pady 0 \
          -command "set ${ns}::cylinder_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cylinder_radius + 0.1]]; ${ns}::Auto_update"
        label $wgo.cylinder_frame.angstrom -text "A"
        pack $wgo.cylinder_frame.entry $wgo.cylinder_frame.decr \
          $wgo.cylinder_frame.incr $wgo.cylinder_frame.angstrom \
          -side left -anchor w -fill x

        grid [button $wgo.coner_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Radius of the arrow cones. For a better representation,\
this value should be larger than arrow cylinder radius."}] \
          -row 11 -column 0 -sticky w
        grid [label $wgo.coner_label -text "Arrow cone radius:"] \
          -row 11 -column 1 -sticky w
        grid [frame $wgo.coner_frame] \
          -row 11 -column 2 -sticky w
        entry $wgo.coner_frame.entry -width 4 -textvariable ${ns}::cone_radius
        button $wgo.coner_frame.decr -text "-0.1" -padx 0 -pady 0 \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius - 0.1]]; ${ns}::Auto_update"
        button $wgo.coner_frame.incr -text "+0.1" -padx 0 -pady 0 \
          -command "set ${ns}::cone_radius \[::tcl::mathfunc::abs \[expr \$${ns}::cone_radius + 0.1]]; ${ns}::Auto_update"
        label $wgo.coner_frame.angstrom -text "A"
        pack $wgo.coner_frame.entry $wgo.coner_frame.decr $wgo.coner_frame.incr \
          $wgo.coner_frame.angstrom -side left -anchor w -fill x

        grid [button $wgo.coneh_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Height of the arrow cones. This height is not counted towards arrow length."}] \
          -row 12 -column 0 -sticky w
        grid [label $wgo.coneh_label -text "Arrow cone height:"] \
          -row 12 -column 1 -sticky w
        grid [frame $wgo.coneh_frame] \
          -row 12 -column 2 -sticky w
        entry $wgo.coneh_frame.entry -width 4 -textvariable ${ns}::cone_height
        button $wgo.coneh_frame.decr -text "-0.2" -padx 0 -pady 0 \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height - 0.2]]; ${ns}::Auto_update"
        button $wgo.coneh_frame.incr -text "+0.2" -padx 0 -pady 0 \
          -command "set ${ns}::cone_height \[::tcl::mathfunc::abs \[expr \$${ns}::cone_height + 0.2]]; ${ns}::Auto_update"
        label $wgo.coneh_frame.angstrom -text "A"
        pack $wgo.coneh_frame.entry $wgo.coneh_frame.decr $wgo.coneh_frame.incr \
          $wgo.coneh_frame.angstrom -side left -anchor w -fill x

        grid [button $wgo.resolution_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "The quality of arrow and protein graphics."}] \
          -row 21 -column 0 -sticky w
        grid [label $wgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wgo.resolution_frame.list ${ns}::resolution 6 
        $wgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution \
              -command "${ns}::Protein_representation \$${ns}::molid; ${ns}::Auto_update"  
        } 
        pack $wgo.resolution_frame.list -side left -anchor w -fill x

        set wpgo [labelframe $w.prograph_options -text "Protein Graphics Options" -bd 2]
        
        grid [button $wpgo.selstr_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, selection string will be effective for\
protein and animation representations."}] \
          -row 0 -column 0 -sticky w
        grid [label $wpgo.selstr_label -text "Show selected atoms:"] \
          -row 0 -column 1 -sticky w
        grid [checkbutton $wpgo.selstr_check -text "" \
            -variable ${ns}::selrep -command "${ns}::Protein_representation \$${ns}::molid"] \
          -row 0 -column 2 -sticky w
        
        grid [button $wpgo.protas_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Protein representation."}] \
          -row 13 -column 0 -sticky w
        grid [label $wpgo.protas_label -text "Show protein as:"] \
          -row 13 -column 1 -sticky w
        grid [frame $wpgo.protas_frame] \
          -row 13 -column 2 -sticky w
        tk_optionMenu $wpgo.protas_frame.list ${ns}::showproteinas "Backbone"
        $wpgo.protas_frame.list.menu delete 0
        $wpgo.protas_frame.list.menu add radiobutton -label "Backbone" -variable ${ns}::showproteinas -command "${ns}::Protein_representation \$${ns}::molid"
        $wpgo.protas_frame.list.menu add radiobutton -label "Network" -variable ${ns}::showproteinas -command "${ns}::Protein_representation \$${ns}::molid"
        pack $wpgo.protas_frame.list -side left -anchor w -fill x

        grid [button $wpgo.procolor_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Color scheme for the protein.\nMobility : protein is colored based on the mobility of residues in the active mode\nBfactors : protein is colored based on the PDB bfactors\nIndex : protein is colored based on residue/atom index"}] \
          -row 14 -column 0 -sticky w
        grid [label $wpgo.procolor_label -text "Color protein:"] \
          -row 14 -column 1 -sticky w
        grid [frame $wpgo.procolor_frame] \
          -row 14 -column 2 -sticky w
        tk_optionMenu $wpgo.procolor_frame.list ${ns}::proteincolor "Mobility"
        $wpgo.procolor_frame.list.menu delete 0
        $wpgo.procolor_frame.list.menu add radiobutton -label "Mobility" -variable ${ns}::proteincolor -command "${ns}::Protein_representation \$${ns}::molid"
        $wpgo.procolor_frame.list.menu add radiobutton -label "Bfactors" -variable ${ns}::proteincolor -command "${ns}::Protein_representation \$${ns}::molid"
        foreach acolor "Index Chain ResName ResType" {
          $wpgo.procolor_frame.list.menu add radiobutton -label $acolor \
              -variable ${ns}::proteincolor \
              -command "${ns}::Protein_representation \$${ns}::molid"
        }
        pack $wpgo.procolor_frame.list -side left -anchor w -fill x


        grid [button $wpgo.cutoffdistance_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Cutoff distance for placing bonds between nodes. Note that this distance only affects representation, not the normal modes contained in the ProDy GUI."}] \
          -row 16 -column 0 -sticky w
        grid [label $wpgo.cutoffdistance_label -text "Network cutoff distance:"] \
          -row 16 -column 1 -sticky w
        grid [frame $wpgo.cutoffdistance_frame] \
          -row 16 -column 2 -sticky w
        entry $wpgo.cutoffdistance_frame.entry -width 4 -textvariable ${ns}::cutoffdistance
        button $wpgo.cutoffdistance_frame.decr -text "-1.0" -padx 0 -pady 0 \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance - 1.0]]; ${ns}::Protein_representation \$${ns}::molid"
        button $wpgo.cutoffdistance_frame.incr -text "+1.0" -padx 0 -pady 0 \
          -command "set ${ns}::cutoffdistance \[::tcl::mathfunc::abs \[expr \$${ns}::cutoffdistance + 1.0]]; ${ns}::Protein_representation \$${ns}::molid"
        label $wpgo.cutoffdistance_frame.angstrom -text "A"
        pack $wpgo.cutoffdistance_frame.entry $wpgo.cutoffdistance_frame.decr $wpgo.cutoffdistance_frame.incr \
          $wpgo.cutoffdistance_frame.angstrom -side left -anchor w -fill x

        grid [button $wpgo.nodescale_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Scale the size of node spheres (or the vdW radius of alpha carbons)."}] \
          -row 17 -column 0 -sticky w
        grid [label $wpgo.nodescale_label -text "Scale node spheres:"] \
          -row 17 -column 1 -sticky w
        grid [frame $wpgo.nodescale_frame] \
          -row 17 -column 2 -sticky w
        entry $wpgo.nodescale_frame.entry -width 4 -textvariable ${ns}::spherescale
        button $wpgo.nodescale_frame.decr -text "-0.1" -padx 0 -pady 0 \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale - 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        button $wpgo.nodescale_frame.incr -text "+0.1" -padx 0 -pady 0 \
          -command "set ${ns}::spherescale \[::tcl::mathfunc::abs \[expr \$${ns}::spherescale + 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        pack $wpgo.nodescale_frame.entry $wpgo.nodescale_frame.decr $wpgo.nodescale_frame.incr \
          -side left -anchor w -fill x
          
        grid [button $wpgo.bondradius_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Radius of the artificial bonds representing springs for connected nodes in the network."}] \
          -row 18 -column 0 -sticky w
        grid [label $wpgo.bondradius_label -text "Network bond radius:"] \
          -row 18 -column 1 -sticky w
        grid [frame $wpgo.bondradius_frame] \
          -row 18 -column 2 -sticky w
        entry $wpgo.bondradius_frame.entry -width 4 -textvariable ${ns}::bondradius
        button $wpgo.bondradius_frame.decr -text "-0.1" -padx 0 -pady 0 \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius - 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        button $wpgo.bondradius_frame.incr -text "+0.1" -padx 0 -pady 0 \
          -command "set ${ns}::bondradius \[::tcl::mathfunc::abs \[expr \$${ns}::bondradius + 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        label $wpgo.bondradius_frame.angstrom -text "A"
        pack $wpgo.bondradius_frame.entry $wpgo.bondradius_frame.decr $wpgo.bondradius_frame.incr \
          $wpgo.bondradius_frame.angstrom -side left -anchor w -fill x


        grid [button $wpgo.tuberadius_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Radius of the tube representation used to depict the protein backbone structure."}] \
          -row 19 -column 0 -sticky w
        grid [label $wpgo.tuberadius_label -text "Backbone tube radius:"] \
          -row 19 -column 1 -sticky w
        grid [frame $wpgo.tuberadius_frame] \
          -row 19 -column 2 -sticky w
        entry $wpgo.tuberadius_frame.entry -width 4 -textvariable ${ns}::tuberadius
        button $wpgo.tuberadius_frame.decr -text "-0.1" -padx 0 -pady 0 \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius - 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        button $wpgo.tuberadius_frame.incr -text "+0.1" -padx 0 -pady 0 \
          -command "set ${ns}::tuberadius \[::tcl::mathfunc::abs \[expr \$${ns}::tuberadius + 0.1]]; ${ns}::Protein_representation \$${ns}::molid"
        label $wpgo.tuberadius_frame.angstrom -text "A"
        pack $wpgo.tuberadius_frame.entry $wpgo.tuberadius_frame.decr $wpgo.tuberadius_frame.incr \
          $wpgo.tuberadius_frame.angstrom -side left -anchor w -fill x


        grid [button $wpgo.material_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "The material used for drawing arrow and protein graphics."}] \
          -row 20 -column 0 -sticky w
        grid [label $wpgo.material_label -text "Graphics material:"] \
          -row 20 -column 1 -sticky w
        grid [frame $wpgo.material_frame] \
          -row 20 -column 2 -sticky w
        tk_optionMenu $wpgo.material_frame.list ${ns}::material "Opaque"
        $wpgo.material_frame.list.menu delete 0
        foreach mtrl "Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy" {
          $wpgo.material_frame.list.menu add radiobutton -label $mtrl \
              -variable ${ns}::material \
              -command "${ns}::Protein_representation \$${ns}::molid; ${ns}::Auto_update"    
        }
        pack $wpgo.material_frame.list -side left -anchor w -fill x

        grid [button $wpgo.resolution_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "The quality of arrow and protein graphics."}] \
          -row 21 -column 0 -sticky w
        grid [label $wpgo.resolution_label -text "Graphics resolution:"] \
          -row 21 -column 1 -sticky w
        grid [frame $wpgo.resolution_frame] \
          -row 21 -column 2 -sticky w
        tk_optionMenu $wpgo.resolution_frame.list ${ns}::resolution 6 
        $wpgo.resolution_frame.list.menu delete 0
        foreach resol "6 10 15 20 25 30 35 40 45 50" {
          $wpgo.resolution_frame.list.menu add radiobutton -label $resol \
              -variable ${ns}::resolution \
              -command "${ns}::Protein_representation \$${ns}::molid; ${ns}::Auto_update"  
        } 
        pack $wpgo.resolution_frame.list -side left -anchor w -fill x

        set wao [labelframe $w.animation_options -text "Animation Options" -bd 2]
        
        grid [button $wao.auto_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Update animation automatically when active mode selection changes."}] \
          -row 0 -column 0 -sticky w
        grid [label $wao.auto_label -text "Auto animate:"] \
          -row 0 -column 1 -sticky w
        grid [checkbutton $wao.auto_check -text "" \
            -variable ${ns}::autoanimate] \
          -row 0 -column 2 -sticky w
        
        #grid [button $wao.overanim_help -text "?" -padx 0 -pady 0 \
        #    -command {tk_messageBox -type ok -title "HELP" \
        #      -message "If checked, animation will be loaded to the same molecule in VMD."}] \
        #  -row 2 -column 0 -sticky w
        #grid [label $wao.overanim_label -text "Overwrite animation:"] \
        #  -row 2 -column 1 -sticky w
        #grid [checkbutton $wao.overanim_check -text "" -variable ${ns}::overanim] \
        #  -row 2 -column 2 -sticky w
        
        grid [button $wao.keepfiles_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, pdb file that contains coordinates used in the animation will be saved."}] \
          -row 3 -column 0 -sticky w
        grid [label $wao.keepfiles_label -text "Keep trajectory file:"] \
          -row 3 -column 1 -sticky w
        grid [checkbutton $wao.keepfiles_check -text "" -variable ${ns}::keepfiles] \
          -row 3 -column 2 -sticky w

        grid [button $wao.autoplay_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, animation will be played continuously."}] \
          -row 4 -column 0 -sticky w
        grid [label $wao.autoplay_label -text "Continuous autoplay:"] \
          -row 4 -column 1 -sticky w
        grid [checkbutton $wao.autoplay_check -text "" -variable ${ns}::autoplay] \
          -row 4 -column 2 -sticky w
          
        grid [button $wao.nframes_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "The resulting animation will have this many plus one frames."}] \
          -row 9 -column 0 -sticky w
        grid [label $wao.nframes_label -text "Number of frames:"] \
          -row 9 -column 1 -sticky w
        grid [frame $wao.nframes_frame] \
          -row 9 -column 2 -sticky w
        entry $wao.nframes_frame.entry -width 4 -textvariable ${ns}::nframes
        button $wao.nframes_frame.decr10 -text "-10" -padx 0 -pady 0 \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes - 10]"
        button $wao.nframes_frame.decr -text "-2" -padx 0 -pady 0 \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes - 2]"
        button $wao.nframes_frame.incr -text "+2" -padx 0 -pady 0 \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes + 2]"
        button $wao.nframes_frame.incr10 -text "+10" -padx 0 -pady 0 \
          -command "set ${ns}::nframes \[expr \$${ns}::nframes + 10]"
        pack $wao.nframes_frame.entry $wao.nframes_frame.decr10 \
          $wao.nframes_frame.decr $wao.nframes_frame.incr \
          $wao.nframes_frame.incr10 -side left -anchor w -fill x

        set wpo [labelframe $w.plotting_options -text "Plotting Options" -bd 2]

        grid [button $wpo.overplot_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "If checked, previously generated canvas will be used for plotting."}] \
          -row 0 -column 0 -sticky w
        grid [label $wpo.overplot_label -text "Overplot:"] \
          -row 0 -column 1 -sticky w
        grid [checkbutton $wpo.overplot_check -text "" \
            -variable ${ns}::overplot] \
          -row 0 -column 2 -sticky w

        grid [button $wpo.plotwidth_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Width of the plot in pixels."}] \
          -row 1 -column 0 -sticky w
        grid [label $wpo.plotwidth_label -text "Plot width:"] \
          -row 1 -column 1 -sticky w
        grid [entry $wpo.plotwidth_entry -width 4 -textvariable ${ns}::plotwidth] \
          -row 1 -column 2 -sticky w
        grid [label $wpo.spacing_label -text "  "] \
          -row 1 -column 3 -sticky w

        grid [button $wpo.plotheight_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Height of the plot in pixels."}] \
          -row 1 -column 4 -sticky w
        grid [label $wpo.plotheight_label -text "Plot height:"] \
          -row 1 -column 5 -sticky w
        grid [entry $wpo.plotheight_entry -width 4 -textvariable ${ns}::plotheight] \
          -row 1 -column 6 -sticky w

        grid [button $wpo.line_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Connect (or not) datapoint with lines."}] \
          -row 3 -column 0 -sticky w
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

        grid [button $wpo.linewidth_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Width of the lines connecting datapoints."}] \
          -row 3 -column 4 -sticky w
        grid [label $wpo.linewidth_label -text "Line width:"] \
          -row 3 -column 5 -sticky w
        grid [entry $wpo.linewidth_entry -width 4 -textvariable ${ns}::linewidth] \
          -row 3 -column 6 -sticky w
          
        grid [button $wpo.marker_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Draw markers at datapoints."}] \
          -row 5 -column 0 -sticky w
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
        
        
        grid [button $wpo.radius_help -text "?" -padx 0 -pady 0 \
            -command {tk_messageBox -type ok -title "HELP" \
              -message "Data point marker radius of circle and point, size of square."}] \
          -row 5 -column 4 -sticky w
        grid [label $wpo.radius_label -text "Marker size:"] \
          -row 5 -column 5 -sticky w
        grid [entry $wpo.radius_entry -width 4 -textvariable ${ns}::mradius] \
          -row 5 -column 6 -sticky w

        #grid [button $wpo.dash_help -text "?" -padx 0 -pady 0 \
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
          
        
        ${ns}::Load_coordinates
        ${ns}::Draw_arrows

        return $w
      }
      nmwizgui

    }

    lappend openfiles $filename
    #lappend namespaces $ns
    #lappend nmwizguis [string range $ns 2 end]
    
    set w .nmwizgui
    
    set wgf [labelframe $w.{[string range $ns 2 end]}frame -text "[subst $${ns}::title]" -bd 2]
    
    grid [button $wgf.name -width 8 -pady 2 -text "Show GUI" \
        -command "${ns}::nmwizgui" ] \
      -row 0 -column 0
    grid [button $wgf.website -width 8 -pady 2 -text "Remove" \
        -command "lset ::nmwiz::openfiles $guicount NONE; pack forget $wgf; ${ns}::Close_mols; namespace delete $ns; destroy .[string range $ns 2 end]"] \
      -row 0 -column 1

    pack $wgf -side top -fill x -expand 1
    
  }
}
# Porcupine
set porcupine_cone_radius 0.5

proc viewModePorcupine {mode direction molname} {
    global refxyz mode_hide_shorter mode_scale_by 
    global mode_materials mode_material mode_color mode_resolution
    global porcupine_cone_radius
    set molid [mol new]
    graphics $molid color $mode_color
    graphics $molid materials $mode_materials
    graphics $molid material $mode_material
    foreach xyz $refxyz v $mode {
        set v [vecscale [expr $direction * $mode_scale_by] $v]
        if {$mode_hide_shorter < [veclength $v]} {
            graphics $molid cone $xyz [vecadd $xyz $v] radius $porcupine_cone_radius resolution $mode_resolution
        }
    }
    mol rename $molid $molname
}

proc nmwiz_tk {} {
  ::nmwiz::init_gui
  return $::nmwiz::w
}

proc nmwiz_load {filename} {
  nmwiz_tk
  ::nmwiz::load_nmd $filename
} 

#nmwiz_tk
