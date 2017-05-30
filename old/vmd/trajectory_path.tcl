
#VMD  --- start of VMD description block
#Name:
# Trajectory path
#Synopsis:
# Draws the path of the center of mass of a selection through an animation
#Version:
# 1.2
#Uses VMD version:
# 1.8
#Ease of use:
# 2
#Procedures:
# <li>trajectory_path selection {color blue} {linewidth 1} {update 0}
# -- follows the center of mass of the given selection.  
# 'color' is a solid color, or "scale" for a color scale,
# 'update' toggles calling "$selection update", and 'linewidth' is 
# the width of the line drawn.
#Description:
# For each step in the animation, the center of mass of the selection is
# calculated.  A line connecting the successive center of mass coordinates
# is then added to the molecule id of the selection. 
# The color is a solid color (default is blue) or they are mapped to the 
# color scale from lowest (= first trajectory frame) to highest (= last 
# trajectory frame).
# The third argument decides (0 = no, 1 = yes) whether the selection
# is updated during the course of following the C.O.M. (default is no).
# The fourth argument allows to specify the width of the line (default 1).
# The procedure returns the graphics ids of the drawn objects.
#Example:
# <pre>
# set water [atomselect top "resid 5243"]
# trajectory_path $water scale
#
# # follow solvation shell around an atom.
# set solv  [atomselect top "water and exwithin 3.8 of index 199"]
# trajectory_path $solv yellow 3 1
#Files: 
# <a href="trajectory_path.vmd">trajectory_path</a>
#Author: 
# Andrew Dalke &lt;dalke@ks.uiuc.edu&gt;
# Axel Kohlmeyer &lt;axel.kohlmeyer@rub.de&gt; (linewidth/update)
#\VMD  --- end of block

proc trajectory_path {selection {color blue} {linewidth 1} {update 0}} {

 # save the current selection frame number and get the molecule id.
 set sel_frame [$selection frame]
 set gr_mol [$selection molindex]

 # make the list of coordinates
 set num_frames [molinfo $gr_mol get numframes]
 set coords {}
 for {set i 0} {$i < $num_frames} {incr i} {
   $selection frame $i
   if {$update} { $selection update }
   # compute the center of mass and save it on the list
   lappend coords [measure center $selection weight mass]
 }
 ##### now make the graphics and store the respective graphic ids in a list.
 set gr_list {}
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

 # return the selection to its original state
 $selection frame $sel_frame
 if {$update} {$selection update}
 return $gr_list
}

