export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`

$FLUPRO/flutil/rfluka -M1 -e rootfluka NA60plus-JPSI.inp
