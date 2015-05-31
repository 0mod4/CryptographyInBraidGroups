load "len.m";
load "mygt.m";

/*--------------------------------------------------------------------*/
function MatrixOrder(a, b, A, my_len)
/**********************************************************************
 * monomial ordering on length vectors. returns 1 if a>b, -1 if a<b and 0 else
 * A square matrix of size #a=#b, a and b sequences containing braid elements,
 * my_len is used to get length vector
 * * *********************************************************************/
   n := #a;
   m := #b;
   assert n eq m;
   assert not IsSingular(A);
   
   a := Vector(n, my_len(a));
   b := Vector(m, my_len(b));
   Aa := a*Transpose(A);
   Ab := b*Transpose(A);

   for i in [1..n] do
	if Aa[i] gt Ab[i] then
		return 1;
        end if;
	if Aa[i] lt Ab[i] then
		return -1;
        end if;
   end for;

   return 0;
end function;

function PeelOff(a, xbx, my_len, matrix, B_N)
/**********************************************************************
   Tries to peel off one generator from xbx by searching the ai in a 
   or an inverse such that axbxa has minimal length.
   returns  the group element ai for which it is minimal,
            the number i of this element (-i if inverse)
            the element after 'peeling off' axbxa
**********************************************************************/
   try
      //run through all ai and inverses and check which gives the shortest vector
      axbxa := [xbx[j]^a[1]: j in [1..#xbx]];
      min_a := a[1];
      min_a_num := 1;
      for i in [2..#a] do
         test := [xbx[j]^a[i]: j in [1..#xbx]];
	 if matrix eq 0 then
         	if (gtAV(axbxa, test, my_len) eq 1) then
            		axbxa := test;
            		min_a := a[i];
            		min_a_num := i;
         	end if;
	 else
		if (MatrixOrder(axbxa, test, matrix, my_len) eq 1) then
			axbxa := test;
			min_a := a[i];
			min_a_num := i;
		end if;
	 end if;
         
         test := [xbx[j]^Inverse(a[i]): j in [1..#xbx]];
         if matrix eq 0 then
		if (gtAV(axbxa, test, my_len) eq 1) then
            		axbxa := test;
           	 	min_a := Inverse(a[i]);
            		min_a_num := -i;
         	end if;
	 else
		if (MatrixOrder(axbxa, test, matrix, my_len) eq 1) then
			axbxa := test;
			min_a := Inverse(a[i]);
			min_a_num := -i;
		end if;
	 end if;
      end for;

      return min_a, min_a_num, axbxa;
   catch e
      print "Error in PeelOff: ",e`Object;
      return Identity(B_N), 1, [Identity(B_N)];
   end try;
end function;

function Attack(a, b, xbx, my_len, matrix, gen_x, B_N)
/**********************************************************************
   Length based attack on xbx
   returns true if successfull, else false
***********************************************************************/
   try
      generators := [];
      for i in [1..gen_x] do
         min_a, min_a_num, axbxa := PeelOff(a, xbx, my_len, matrix, B_N);
         if matrix eq 0 then
         	if (gtAV(axbxa, xbx, my_len) eq 1) then
            		return 0;
         	else
            		xbx := axbxa;
            		Append(~generators, min_a_num);
         	end if;
         else
		if (MatrixOrder(axbxa, xbx, matrix, my_len) eq 1) then
			return 0;
		else
			xbx := axbxa;
			Append(~generators, min_a_num);
		end if;
         end if;
      end for;
      //Test if all generators have been peeled off correctly
      if (b eq axbxa) then
         return 1;
      else
         return 0;
      end if;
   catch e
      print "Error in Attack: ",e`Object;
      return -1;
   end try;
end function;

/*--------------------------------------------------------------------*/
function Test(N, n, m, gen_a, gen_b, gen_x, Matrices)
/**********************************************************************
   Sets up an environment for length based attack and then attacks
   
   returns  1 if was successfull
            0 if was not successfull
            -1 if error
**********************************************************************/
   try
      //Build Braid Group
      B_N := BraidGroup(N);
      //Generate a_i and b_i
      a := [Random(B_N, gen_a[i]): i in [1..#gen_a]];
      b := [Random(B_N, gen_b[i]): i in [1..#gen_b]];
      //Generate x
      x := &*[Random(a):j in [1..gen_x]];
      
      xbx := [b[i]^x: i in [1..#b]];
      
      results := [];
      result_can := Attack(a, b, xbx, LenCan, 0, gen_x, B_N);
      result_gar := Attack(a, b, xbx, LenGarside, 0, gen_x, B_N);//uses gtAV if given Matrix is 0
      Append(~results, result_can);
      Append(~results, result_gar);
      print "gtAV: LenCan: ",result_can," LenGarside: ",result_gar;

      //Try with the given matrices
      for i in [1..#Matrices] do
        result_can := Attack(a, b, xbx, LenCan, Matrices[i], gen_x, B_N);
	Append(~results, result_can);
        result_gar := Attack(a, b, xbx, LenGarside, Matrices[i], gen_x, B_N);
        Append(~results, result_gar);
        print "Matrix Nr.",i,": LenCan: ",result_can," LenGarside: ",result_gar;
      end for;      

      return results;
   catch e
      print "error in Test";
      error e`Object;
      return -1;
   end try;
   
end function;


/*--------------------------------------------------------------------*/
function NumSuccess    (T:
                        Q := 10,
                        n := 20,
                        m := 20,
                        gen_x := 30,
                        my_gt := gtAV, 
                        gen_a := [10 : i in [1..n]], gen_b := [10 : i in [1..m]])
/**********************************************************************
   Coordinates T Runs of the attack, returns howmany have been successfull. Uses
   Monomial orderings based on random matrices, compares go average based ordering
   N:       Working in BraidGroup B_N (N-1 Artin generators) Set in code
   n:       number of a_i, Generators of the subgroup where conjugator x is from
   m:       number of b_i, Elements that are conjugated with x
   gen_a:   number of generators for the a_i
   gen_b:   number of generators for the b_i
   gen_x:   number of a_i that generate x
   len:     length function for the elements in B_N
   gt:      linear ordering on vectors containing elements from B_N using len
**********************************************************************/
   print "n=",n," m=",m," LenCan vs LenGarside, gtAV, Matrix gt, gen_a=gen_b=10";
   print "gen_x=9, N=9, T=",T," tries and Q=",Q," different Matrices";
   gen_x := 9;
   N := 9;

   SetSeed(Round(Realtime()));

   //create random matrices
   print "matrices contain integer values between -20 and 20";
   Matrices := [];
   j := 1;
   while j le Q do
        Append(~Matrices, RandomGLnZ(m, 20, 20));
        j := j+1;
   end while;

   print Matrices;

   SetSeed(42);

   i:=1;
   results := []; //enthält tupel mit Ergebnissen
   while i le T do
	result := Test(N, n, m, gen_a, gen_b, gen_x, Matrices);
   	Append(~results, result); 
	print "Test ",i," finished";
	i := i+1;
   end while;

   //Ergebnisse verarbeiten
   MatGarResults := [];
   MatCanResults := [];
   for i in [1..Q] do //one entry for each matrix
	Append(~MatGarResults, [0,0,0,0]); //0eqAV, 1eqAV, 0neAV, 1neAV
        Append(~MatCanResults, [0,0,0,0]); //0eqAV, 1eqAV, 0neAV, 1neAV
   end for;

   for i in [1..#results] do
   	cur_row := results[i];
	AVCan := cur_row[1];
        AVGar := cur_row[2];
        matnr := 1;
        j := 3;
	while j lt #cur_row do
		matCan := cur_row[j];
		matGar := cur_row[j+1];

		if matCan ne AVCan then
			if matCan eq 1 then
				MatCanResults[matnr][4]+:=1;
			else
				MatCanResults[matnr][3]+:=1;
			end if;
		else
			if matCan eq 1 then
				MatCanResults[matnr][2]+:=1;
			else
				MatCanResults[matnr][1]+:=1;
			end if;
		end if;

                if matGar ne AVGar then
                        if matGar eq 1 then
                                MatGarResults[matnr][4]+:=1;
                        else
                                MatGarResults[matnr][3]+:=1;
                        end if;
                else
                        if matGar eq 1 then
                                MatGarResults[matnr][2]+:=1;
                        else
                                MatGarResults[matnr][1]+:=1;
                        end if;
                end if;
		matnr := matnr +1;
		j := j+2;
	end while;
   end for;

   //Print results
   for i in [1..#MatCanResults] do
	print "Matrix Nr.",i," Canonical length: ",MatCanResults[i][1],"x 0eqAV, ",MatCanResults[i][2],"x 1eqAV, ",MatCanResults[i][3],"x 0neAV, ",MatCanResults[i][4],"x 1neAV";
        print "Matrix Nr.",i," Garside length: ",MatGarResults[i][1],"x 0eqAV, ",MatGarResults[i][2],"x 1eqAV, ",MatGarResults[i][3],"x 0neAV, ",MatGarResults[i][4],"x 1neAV";
   end for;

   return 1;
end function;
