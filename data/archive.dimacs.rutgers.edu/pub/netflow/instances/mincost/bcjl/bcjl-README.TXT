The file bcjl-data.tar.Z contains the following data files:
capt-TightLoose-in
capt-TightLoose-out
capt-TightLoose-sum
capt-big-in
capt-big-out
capt-big-sum
capt-rnet_param_var
capt-sparse-in
capt-sparse-out
capt-sparse-sum
goto-big-in
goto-big-out
goto-big-sum
grid-big-in
grid-big-out
grid-big-sum
netg-TightLoose-in
netg-TightLoose-out
netg-TightLoose-sum
netg-big-in
netg-big-out
netg-big-sum
netg-rnet_param_var
netg-sparse-in
netg-sparse-out
netg-sparse-sum

After sending bcjl-data.tar.Z by ftp to a machine and directory where you
have write-permission, you can retrieve the files above; give the following
commands:

		uncompress bcjl-data.tar.Z
		tar -xf bcjl-data.tar

The naming of the data files conforms to the following convention. The name
consists of two or three parts separated by hyphens. The first part denotes
the generator used, the second part identifies the test, and the third part
specifies the content.

Generators
	capt  (Bland and Jensen)

	netg = netgen (Napier, Klingmann and Stutz)
               
	grid = gridgen (Lee and Orlin)

	goto = grid-on-torus (Goldberg)

		capt, gridgen, and grid-on-torus can be found in the dimacs
		directory
			 /pub/netflow/generators/network
		The version of netgen that we used can be found in the
		dimacs directory 
			/pub/netflow/submit/netgen-bcjl
		This version of netgen contains modifications to avoid
		overflows when n is large and to allow the highest-
		numbered sink arc to be incident to non-skeleton arcs.

Tests
	big:  capt-big-* and netg-big-* are observations from the main
	      experiment (see section 3 of the paper). The first 100
	      observations are from the base (estimation) phase, the second
	      100 observations from the validation phase.  grid-big-* and
	      goto-big-* are not part of the main experiment (see section
	      4).  grid-big-* follows the same form as above.  In the case
	      of goto-big-* there are only 100 observations; no validation
	      observations were taken due to very long run-times and poor
	      fit in the estimation phase.

	TightLoose:  "re-runs" of t<he 100 observations from the estimation 
                      phase of the big test, once with tight capacities and
		      once with loose capacities (see the last paragraph of
	              section 3.2).

	sparse:  20 large sparse observations (not part of the main
		 experiment -- see section 5.2.

	rnet_param_var:  repeated trials of the network simplex code on 2
			 instances using RNET 3.61 routines with different
			 choices of the frequency parameter, frq, and the
			 parameter p_0; (not part of the main experiment --
			 see section 5.4).  These files are in a form
			 different from the others -- the values of the
			 control variables and generator inputs can be read
			 by referencing other files.
			 

Content
	in:   file contains generator inputs.

	out:  file contains detailed output from the four solvers.

	sum:  file contains summary output.



	The output files report all run times in seconds. The abbreviations
	`rlx' and `pmin' are used for `RELAXT-III' and `MIN3' respectively. 
	For each network the output is given for each of the four
	codes always in the order: RNET 3.61, RELAXT-III, MIN3, SCALE.
        The last four entries in each line of the summary output files give 
	the four `ok-bits' in the order above.	The ok_bit for a code
	on a given network is set to zero if the cost of the solution
        computed by that code matches the cost of RNET 3.61's computed
        solution; otherwise it is set to one. Since the first ok_bit is for 
	RNET 3.61 it is zero by definition. The other values of zero reveal 
	that all four codes are terminating with the same reported optimal 
	value. The number of saturated arcs in the computed solution is 
	reported for each code as `nsat'.



Please direct any questions or requests to 

	bland@orie.cornell.edu
or
	ladanyi@cs.cornell.edu
