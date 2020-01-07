# To unbundle, sh this file
echo min_to_flo 1>&2
sed 's/.//' >min_to_flo <<'//GO.SYSIN DD min_to_flo'
-
-
-case $# in
-0|1)	echo "Usage: min_to_flo input.min NETFLO.INP"; exit 1
-esac
-
-
-awk -f trans.a <$1 > tmp.0000
-grep a tmp.0000 > tmp.a
-grep b tmp.0000 > tmp.b
-grep c tmp.0000 > tmp.c
-grep d tmp.0000 > tmp.d
-grep e tmp.0000 > tmp.e
-grep f tmp.0000 > tmp.f
-cat tmp.[abcdef] | awk -f format.a >$2
-\rm tmp.0000 tmp.[abcdef]
//GO.SYSIN DD min_to_flo
echo trans.a 1>&2
sed 's/.//' >trans.a <<'//GO.SYSIN DD trans.a'
-#This awk program is the first part of a system to translate
-#files in the DIMACS .min format to ones readable by Helgassen 
-#and Kennington's NETFLO program.  
-#It does no error checking. 
-
-$1 == "c" { #ignore comment lines 
-	}
-
-$1 == "p"  {prob = $2;  nodes = $3; arcs = $4;
-	    print "a", nodes  
-	}
-
-$1 == "n"  {name = $2; demand = $3;
-	    print "b", name, demand 
-	}
-
-$1 == "a"  {from = $2; to =$3; low=$4; cap=$5; cost=$6; 
-	    acount++;
-
-	    if (cap == -1) cap = 0;
-            else if (cap == 0) cap = -1;
-
-            print "e", acount, from, to, cost, cap, low;  
-	
-		incident[to]++;
-	}	
-
-END     { for (i = 1; i <= nodes; i++) {
-		if (1== (i % 8) ) printf("d ");    #start new line 
-		printf("%10d",incident[i]);
-		if (0== (i % 8) ) printf("\n");    #terminate line
-	   } 
-          #finish last line 
-          while (0 != (i % 8)) {printf("%10d", 0); i++;}
-	  printf("\n");
-	  
- 	  # add zero lines
-	  print "c\n"; 
-          print "f\n" ; 
-	}
-	 
-	   
//GO.SYSIN DD trans.a
echo format.a 1>&2
sed 's/.//' >format.a <<'//GO.SYSIN DD format.a'
-$1 == "a" {printf("%10d\n", $2);}
-
-$1 == "b" {printf("%10d%10d\n", $2, $3); }
-
-$1 == "c" {printf("%-80s\n","00000000000000000000");}
-
-$1 == "d" {printf("%10d%10d%10d%10d%10d%10d%10d%10d\n", $2,$3,$4,$5,$6,$7,$8,$9);
-}
-
-$1 == "e" {printf("%10d%10d%10d%10d%10d\n",$2,$3,$4,$5,$6)}
-
-$1 == "f" {for (i=1; i<=5;i++) printf("%10s","0000000000");
-           printf("%30s\n", " ");
-    	}
-
-          
-
//GO.SYSIN DD format.a

