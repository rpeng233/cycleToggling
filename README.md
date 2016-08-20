First run

$./compile

to compile the data generators and solvers, then run

$./runall graphSize rhsType hopCount

to run each solver on each graph. The type of right hand side used is either
"01" (one unit of flow from one endpoint to the other) or "Rand" (random).
The graphs, the corresponding matrix market files, and the right hand side
files are stored under ./graphdata, and the experiment results are stored under
./results_{graphSize}_{rhsType}.
