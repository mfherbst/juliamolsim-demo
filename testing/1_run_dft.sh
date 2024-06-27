#!/bin/bash
THISDIR=`dirname ${BASH_SOURCE[0]}`
$THISDIR/run_Al_bulk_1.jl
$THISDIR/run_Al_defect_1.jl
$THISDIR/run_Al_supercells_1.jl
$THISDIR/run_Al_supercells_2.jl
