function Duplicates    (T:
                        N := 81,
                        n := 20,
			gen := [10 : i in [1..n]])
/**********************************************************************
   Generates T Times the generation of n random elements in B_N, each of them 
   generated by gen[i], and returns the sum of duplicates over all runs
**********************************************************************/
   results := [];
   B_N := BraidGroup(N);
   a := [];
   duplicatesSum := 0;

   for run in [1..T] do
	   for i in [1..n] do
	      results[i] := i;
	      a[i] := &*[Random(B_N):j in [1..gen[i]]];
	      for k in [1..i-1] do
		  if AreIdentical(NormalForm(a[i]),NormalForm(a[k])) then
		      results[k] := i;
		  end if;
	      end for;
	   end for;

	   duplicates := #results - #SequenceToSet(results);
           duplicatesSum +:= duplicates;
	   print run,"th run with ",duplicates," duplicates";
           if duplicates gt 0 then
		print results;
           end if;
	   results := [];
           duplicates := 0;
   end for;

   return duplicatesSum;
end function;
