1) recompile the code with
   -DAIR_PRESSURE
2) make initial
   use mk_initial.f90
3) use mk_depth.F to generate depth.txt

4) cases

  case_01, static case, Mglob=128, dx=5m,xc=300m, width=100m, amplitude=1m
  case_02, dynamic case,Mglob=512, dx=5m,xc=1000m, width=100m, amplitude=1m, u=11m

