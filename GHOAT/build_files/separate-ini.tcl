mol load pdb hhhh-gggg.pdb
set all [atomselect top "(resid FIRST to LAST) or (resname MMM)"]
set a [atomselect top "(not resname MMM) and (resid FIRST to LAST)"]
set b [atomselect top "resname MMM"]
$all writepdb hhhh-gggg.pdb
$a writepdb hhhh.pdb
$b writepdb gggg.pdb
exit

