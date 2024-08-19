set file_path [lindex $argv 0]
set output_path [lindex $argv 1]

mol new $file_path
mol modstyle 0 0 NewCartoon

set id_matrix {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}
molinfo top set center_matrix [list $id_matrix]
molinfo top set rotate_matrix [list $id_matrix]
molinfo top set scale_matrix [list $id_matrix]

rendering_method stl
render STL $output_path

exit