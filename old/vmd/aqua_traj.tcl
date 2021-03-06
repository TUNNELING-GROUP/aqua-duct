
source trajectory_path.tcl

########################################################################
# settings of scope and objects

#set traj_object "water within 4 of (resid 88 90 136)"
#set traj_object "water within 7 of (resid 375)"
set traj_object "water within 6 of ((resid 88 90 136) and name CA)"
set traj_scope "water within 4 of protein"

proc get_object {} {
	return "water within 6 of ((resid 88 90 136) and name CA)"
}

proc get_scope {} {
	return "water within 4 of protein"
}

########################################################################
# settings of frames and frame step 
# as for now only frame_step 1 is working

# analyze every frame_step frame
proc frame_step {} {
    return 1
}

# analyze frames 0 to frames
proc frames_total {} {
    set frames [molinfo top get numframes]
    #set frames 200
    
    incr frames
    set frames [expr {[expr {$frames / [frame_step]}] * [frame_step]}]
    incr frames -1
    return $frames
}

########################################################################
# vmd related helpers of sel and chres

# generate chres out of sel
# chres is chain:resid sel is atomselect 
proc sel_to_chres {sel} {
    set chres_out ""
    foreach chres [lsort -uniq [$sel get {chain resid}]] {
        lappend chres_out [chres_colonize $chres]
    }
    return $chres_out
}

# add colon beween chain and resid
proc chres_colonize {chres} {
    set colon ""
    append colon [lindex $chres 0]
    append colon ":"
    append colon [lindex $chres 1]
    return $colon
}

# is any chres of sel1 present in sel2?
proc is_sel1_in_sel2 {sel1 sel2} {
    set chres2 [sel_to_chres $sel2]
    foreach chres1 [sel_to_chres $sel1] {
        if {[lsearch $chres2 $chres1] != -1} {
            return 1
        }
    }
    return 0
}

# is chres in sel
proc chres_in_sel {chres sel} {
    set sel_chres [sel_to_chres $sel]
    if {[lsearch $sel_chres $chres] != -1} {
            return 1
        }
    return 0
}

# is not chres in sel
proc chres_not_in_sel {chres sel} {
    return [expr [chres_in_sel $chres $sel]==0]
}

########################################################################
# find all trackable objects in all frames

proc find_objects_in_traj {} {
    # create selections
    set object [atomselect top [get_object]]
    # objects to trace
    set tracable ""
    # loop over frames, increment
    for {set f 0} {$f <= [frames_total]} {incr f} {
        # update selections to current frame
        $object frame $f
        $object update
        append tracable " " [sel_to_chres $object]
    }
    lsort -uniq $tracable
}

########################################################################
# plot helpers

proc plot_all_chres_out_path {traj_object traj_scope} {


    set chres_wody [find_objects_in_traj $traj_object]
    
    foreach woda $chres_wody {
        #puts $woda
        set path [out_path_of_chres $woda $traj_object $traj_scope]
        #puts $path
        plot_chres_path $woda $path green
    }

}

proc plot_chres_path {chres path {color blue} {linewidth 1} {update 0}} {
    
    set chain [lindex [split $chres ":"] 0]
    set resid [lindex [split $chres ":"] 1]
    set selection [atomselect top "chain $chain and resid $resid"]
    
    set gr_mol [$selection molindex]
    
    set coords {}
    foreach i [split $path] {
        #puts "frame $i"
        $selection frame $i
        if {$update} { $selection update }
        # compute the center of mass and save it on the list
        lappend coords [measure center $selection weight mass]
    }
    ##### now make the graphics and store the respective graphic ids in a list.
    set gr_list {}
    #puts $coords
    #puts [llength $coords]
    #puts $coords
    #puts [llength $coords]
	set coords [smooth_coords $coords]
    set coords [lassign $coords prev]

    # use the color scale?
    if {$color == "scale"} {
        set count 0
        incr num_frames
        foreach coord $coords {
            set color [expr [colorinfo num] + int([colorinfo max] * ($count + 0.0) / ($num_frames + 0.0))]
            graphics $gr_mol color $color
            lappend gr_list [graphics $gr_mol line $prev $coord width $linewidth]
            set prev $coord
            incr count
        }
    } else {
        # constant color
        graphics $gr_mol color $color
        foreach coord $coords {
            lappend gr_list [graphics $gr_mol line $prev $coord width $linewidth]
            set prev $coord
        }
    }
} 


proc plot_chres_paths_ico {chres path_in path_core path_out} {
    
    set linewidth 4
    set update 0
    
    set chain [lindex [split $chres ":"] 0]
    set resid [lindex [split $chres ":"] 1]
    set selection [atomselect top "chain $chain and resid $resid"]
    
    set gr_mol [$selection molindex]
    
    set path [concat $path_in $path_core $path_out]
    
    set coords {}
    foreach i [split $path] {
        #puts "frame $i"
        $selection frame $i
        if {$update} { $selection update }
        # compute the center of mass and save it on the list
        lappend coords [measure center $selection weight mass]
    }
    ##### now make the graphics and store the respective graphic ids in a list.
    set gr_list {}
    #puts $coords
    #puts [llength $coords]
    #puts $coords
    #puts [llength $coords]
	set coords [smooth_coords $coords]
    set coords [lassign $coords prev]

    # constant color
    for {set i 0} {$i < [expr {[llength $path] - 1}]} {incr i} {
        set coord [lindex $coords $i]
        #puts $coord
        set p [lindex $path $i]
        if {$p in $path_in} {set color red}
        if {$p in $path_core} {set color blue}
        if {$p in $path_out} {set color green}
        graphics $gr_mol color $color
        lappend gr_list [graphics $gr_mol line $prev $coord width $linewidth]
        set prev $coord
    }
} 

########################################################################
# smooth coords

proc smooth_coords {coords {window 10}} {
    
    set n [llength $coords]
	set half [expr {$window / 2}]

    if {$n > $window} {
        set output {}
        for {set i 0} {$i < $n} {incr i} {
			set phalf $half
			if {$i < $half} {
				set phalf [expr {$i + 1}]
			}
			if {[expr {$n - $i}] < $half} {
				#set phalf [expr {[expr {$n - $i}] + 0}]
				set phalf [expr {$n - $i}]
			}
			#puts "i: $i half: $half phalf: $phalf"
			set lo [expr {$i - $phalf}]
			set hi [expr {$i + $phalf}]
			#puts "lo: $lo hi: $hi"
			#puts [lrange $coords $lo $hi]
			lappend output [mean_coords [lrange $coords $lo $hi]]
        }
		return $output
    } else {
        return $coords
    }
}

proc mean_coords {coords} {
    set n [llength $coords]
    set x 0.
    set y 0.
    set z 0.
    foreach coord $coords {
        set x [expr {$x + [lindex $coord 0]}]
        set y [expr {$y + [lindex $coord 1]}]
        set z [expr {$z + [lindex $coord 2]}]
    }
    set x [expr {$x / $n}]
    set y [expr {$y / $n}]
    set z [expr {$z / $n}]
    return "$x $y $z"
}

########################################################################
# get coords of chres for given frame

proc get_chres_frame_coord {chres f} {
	
    set chain [lindex [split $chres ":"] 0]
    set resid [lindex [split $chres ":"] 1]
    set selection [atomselect top "chain $chain and resid $resid"]
    $selection frame $f
	return [measure center $selection weight mass]
}

########################################################################



########################################################################
# boolean like list manipulations

# return overlap of A and B
proc union {A B} {
    set output {}
    foreach ael $A {
        if {$ael in $B } {
            lappend output $ael
        }
    }
    return $output
}

# glue if A and B overlaps of touch each other
proc glue {A B} {
    if {[lindex $A end] >= [expr {[lindex $B 0] - 1}]} {
        return [lsort -uniq [concat $A $B]]
    }
}

# xor of A and B
proc xor {A B} {
    set AB [union $A $B]
    set output {}
    foreach el [glue $A $B] {
        if {$el ni $AB} {
            lappend output $el
        }
    }
    return $output
}

# assumes A precedes B, returns A part only
proc left {A B} {
    return [union $A [xor $A $B]]
}

# assumes A precedes B, returns B part only
proc right {A B} {
    return [union $B [xor $A $B]]
}

########################################################################

proc plot_paths_of_all {traj_object traj_scope} {

    foreach chres [find_objects_in_traj $traj_object] {
        plot_paths_of_chres $chres $traj_object $traj_scope
    }
}

proc old_plot_paths_of_chres {chres traj_object traj_scope} {
    
    set paths_in [in_path_of_chres $chres $traj_object $traj_scope]
    set paths_out [out_path_of_chres $chres $traj_object $traj_scope]
    puts $chres
    foreach path_in $paths_in {
        set path_out [lindex $paths_out 0]
        set path_glue [glue $path_in $path_out]
        if {[llength $path_glue] > 0} {
            # pop 0
            set paths_out [lassign $paths_out path_out]
            set path_core [union $path_in $path_out]
            set path_in [left $path_in $path_out]
            set path_out [right $path_core $path_out]
        } else {
            # in only, no out
            # to detect core one has to check if it is in traj_object...
            # skipp it as for now
            set path_core {}
            set path_out {}
        }
        puts "\tin:   $path_in"
        puts "\tcore: $path_core"
        puts "\tout:  $path_out"
        plot_chres_paths_ico $chres $path_in $path_core $path_out
    }
    set path_in {}
    set path_core {}
    # lefovers?
    foreach path_out $paths_out {
        puts "\tin:   $path_in"
        puts "\tcore: $path_core"
        puts "\tout:  $path_out"
        plot_chres_paths_ico $chres $path_in $path_core $path_out
    }        
}

proc trace_all_paths {} {
	set chress [find_objects_in_traj]
	#set chress "X:12647"
    # create selection
    set object [atomselect top [get_object]]
    set scope [atomselect top [get_scope]]
    # containers
	array set paths {}
	# init container
	foreach chres $chress {
		set paths($chres) {{} {} {}}
	}
    # loop over frames, search of frames types and coords
    puts "read frames types coords"
    for {set f 0} {$f <= [frames_total]} {incr f [frame_step]} {
		puts "frame $f of [frames_total]"
        # update selections to current frame
        $object frame $f
        $object update
        $scope frame $f
        $scope update
		# loop over chres
		foreach chres $chress {
			#puts "chres $chres"
			if {[chres_in_sel $chres $scope]} {
				lassign $paths($chres) frames types coords
				# save frame
				lappend frames $f
				# type scope or object
				if  {[chres_in_sel $chres $object]} {
					lappend types object
				} else {
					lappend types scope
				}
				# get and save coords
				lappend coords [get_chres_frame_coord $chres $f]
				# store it back in paths
				set paths($chres) "{$frames} {$types} {$coords}"
			}
		}
	}
	# now, for each chres find in and out paths
	puts "find out and in paths"
	foreach chres $chress {
		puts $chres
		# get path for chres
		lassign $paths($chres) frames types coords
		find_and_plot_ico_paths $frames $types $coords
	}
}


proc find_and_plot_ico_paths {frames types coords} {
    
	set paths_out [get_out_paths $frames $types]
	set paths_in [get_in_paths $frames $types]

    foreach path_in $paths_in {
        set path_out [lindex $paths_out 0]
        set path_glue [glue $path_in $path_out]
        if {[llength $path_glue] > 0} {
            # pop 0
            set paths_out [lassign $paths_out path_out]
            set path_core [union $path_in $path_out]
            set path_in [left $path_in $path_out]
            set path_out [right $path_core $path_out]
        } else {
            # in only, no out
            # to detect core one has to check if it is in traj_object...
            # skipp it as for now
            set path_core {}
            set path_out {}
        }
        plot_paths_ico $path_in $path_core $path_out $frames $coords
    }
    set path_in {}
    set path_core {}
    # lefovers?
    foreach path_out $paths_out {
        plot_paths_ico $path_in $path_core $path_out $frames $coords
    }        
}

proc plot_paths_ico {path_in path_core path_out frames all_coords} {
    
    set linewidth 2
    set update 0
    
	set all_coords [smooth_coords $all_coords]

    set gr_mol [molinfo top]
    
    set path [concat $path_in $path_core $path_out]
    
    set coords {}
    foreach i [split $path] {
        #puts "frame $i"
        lappend coords [lindex $all_coords [lsearch $frames $i]]
    }
    ##### now make the graphics and store the respective graphic ids in a list.
    set gr_list {}
    set coords [lassign $coords prev]

    # constant color
    for {set i 0} {$i < [expr {[llength $path] - 1}]} {incr i} {
        set coord [lindex $coords $i]
        #puts $coord
        set p [lindex $path $i]
        if {$p in $path_in} {set color red}
        if {$p in $path_core} {set color blue}
        if {$p in $path_out} {set color green}
        graphics $gr_mol color $color
        lappend gr_list [graphics $gr_mol line $prev $coord width $linewidth]
        set prev $coord
    }
} 


########################################################################
# out and in path detection
proc get_out_paths {frames types} {
    # containers
    set paths {}
    set current_path {}
    # loop over frames, increment
    for {set f 0} {$f <= [frames_total]} {incr f [frame_step]} {
		# get type
		if {$f in $frames} {
			set t [lindex $types [lsearch $frames $f]]
		} else {
			set t out
		}
		#puts "$f $t"
		# is out of scope?
		if {$t == "out"} {
            if {[llength $current_path] > 0} {
                lappend paths $current_path
                set current_path ""
            }
		} else {
            # chres in scope
            # is already tracked?
            if {[llength $current_path] > 0} {
                lappend current_path $f
            } else {
                # is in object zone?
                if {$t == "object"} {
                    lappend current_path $f
                }
            }
        }
    }
    # is anything to add?
    if {[llength $current_path] > 0} {
        lappend paths $current_path
    }
    return $paths
}

proc get_in_paths {frames types} {
    # containers
    set paths {}
    set current_path {}
    # loop over frames, increment
    for {set f [frames_total]} {$f >= 0} {incr f [expr {[frame_step] * -1}]} {
		# get type
		if {$f in $frames} {
			set t [lindex $types [lsearch $frames $f]]
		} else {
			set t out
		}
		#puts "$f $t"
		# is out of scope?
		if {$t == "out"} {
            if {[llength $current_path] > 0} {
                lappend paths [lreverse $current_path]
                set current_path ""
            }
		} else {
            # chres in scope
            # is already tracked?
            if {[llength $current_path] > 0} {
                lappend current_path $f
            } else {
                # is in object zone?
                if {$t == "object"} {
                    lappend current_path $f
                }
            }
        }
    }
    # is anything to add?
    if {[llength $current_path] > 0} {
        lappend paths [lreverse $current_path]
    }
    return [lsort $paths]
}



proc out_path_of_chres {chres traj_object traj_scope} {
    # create selection
    set object [atomselect top $traj_object]
    set scope [atomselect top $traj_scope]
    # containers
    set paths ""
    set current_path ""
    # loop over frames, increment
    for {set f 0} {$f <= [frames_total]} {incr f [frame_step]} {
        # update selections to current frame
        $object frame $f
        $object update
        $scope frame $f
        $scope update
        # is out of scope?
        if {[chres_not_in_sel $chres $scope]} {
            # chres out of scope
            if {[llength $current_path] > 0} {
                lappend paths $current_path
                set current_path ""
            }
        } else {
            # chres in scope
            # is already tracked?
            if {[llength $current_path] > 0} {
                lappend current_path $f
            } else {
                # is in object zone?
                if {[chres_in_sel $chres $object]} {
                    lappend current_path $f
                }
            }
        }
    }
    # is anything to add?
    if {[llength $current_path] > 0} {
        lappend paths $current_path
    }
    return $paths
} 

proc in_path_of_chres {chres traj_object traj_scope} {
    # create selection
    set object [atomselect top $traj_object]
    set scope [atomselect top $traj_scope]
    # containers
    set paths ""
    set current_path ""
    # loop over frames, increment
    for {set f [frames_total]} {$f >= 0} {incr f [expr {[frame_step] * -1}]} {
        # update selections to current frame
        $object frame $f
        $object update
        $scope frame $f
        $scope update
        # is out of scope?
        if {[chres_not_in_sel $chres $scope]} {
            # chres out of scope
            if {[llength $current_path] > 0} {
                lappend paths [lreverse $current_path]
                set current_path ""
            }
        } else {
            # chres in scope
            # is already tracked?
            if {[llength $current_path] > 0} {
                lappend current_path $f
            } else {
                # is in object zone?
                if {[chres_in_sel $chres $object]} {
                    lappend current_path $f
                }
            }
        }
    }
    # is anything to add?
    if {[llength $current_path] > 0} {
        lappend paths [lreverse $current_path]
    }
    return [lsort $paths]
} 


