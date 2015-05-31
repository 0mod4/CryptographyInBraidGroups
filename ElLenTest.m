/******************************************************************
  returns length distribution of random braid elements concerning 
  different length functions
*************************************************************************/
function ElLenTest(N,T)
	SetSeed(Truncate(Realtime()));

	SupLens := AssociativeArray();
	GarsideLens := AssociativeArray();
	CanLens := AssociativeArray();
	InfLens := AssociativeArray();
	CheckLens := AssociativeArray();

	B_N := BraidGroup(N);

	//calculates the lengths of T random braid elements from B_N
	for i in [1..T] do
		print i;
		k := Random([1000..5000]);
		elt:=Random(B_N, k);
		SupLen := Supremum(elt);
		GarsideLen := #NormalForm(elt);
		CanLen := CanonicalLength(elt);
		InfLen := Infimum(elt);
		CheckLen := CanLen - (SupLen-InfLen);

		if IsDefined(SupLens, SupLen) then
			SupLens[SupLen] +:= 1;
		else
			SupLens[SupLen] := 1;
		end if;

		if IsDefined(GarsideLens, GarsideLen) then
			GarsideLens[GarsideLen] +:= 1;
		else
			GarsideLens[GarsideLen] := 1;
		end if;	

       		if IsDefined(CanLens, CanLen) then
			CanLens[CanLen] +:= 1;
		else
			CanLens[CanLen] := 1;
		end if;

		if IsDefined(InfLens, InfLen) then
			InfLens[InfLen] +:= 1;
		else
			InfLens[InfLen] := 1;
		end if;

		if IsDefined(CheckLens, CheckLen) then
			CheckLens[CheckLen] +:= 1;
		else
			CheckLens[CheckLen] := 1;
		end if;

	end for;

	//Print result
	print "Supremum";
	for k in Keys(SupLens) do
		print "Length ",k,": ",SupLens[k];
	end for;

	print "Infimum";
	for k in Keys(InfLens) do
		print "Length ",k,": ",InfLens[k];
	end for;

	print "Canonical";
	for k in Keys(CanLens) do
		print "Length ",k,": ",CanLens[k];
	end for;	

	print "Garside";
	for k in Keys(GarsideLens) do
		print "Length ",k,": ",GarsideLens[k];
	end for;

	print "Check: Canonical - (Sup-Inf)";
	for k in Keys(CheckLens) do
		print "Length ",k,": ",CheckLens[k];
	end for;

	return 1;
end function;
