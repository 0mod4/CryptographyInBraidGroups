load "len.m";
load "mygt.m";
load "my_min.m";

/*--------------------------------------------------------------------*/
function GenRandElts(gen_num, B_N) 
/**********************************************************************
   Generates Sequence of Random elements in Braidgroup Bn of length l
   where Entry i consists of gen_num[i] generators
**********************************************************************/
   a := [Random(B_N, gen_num[i]): i in [1..#gen_num]];
   return a;
end function;

/*--------------------------------------------------------------------*/

function MyMin(LenVec)
	n := #LenVec;
	if n eq 0 then 
		return 0;
	elif n eq 1 then
		return LenVec[1];
	end if;
	m := #LenVec[1];

	avlen := [&+LenVec[j]: j in [1..#LenVec]];

	min, mink := Minimum(avlen);
	return min, mink;
end function;


function PeelOff(xbx, a, mygt, mylen)
/**********************************************************************
   Peels off the generator ai from xbx
   returns  the place of the resulting vector in 
**********************************************************************/
   try
      //run through all ai and inverses
      axbxa := [[xbx[i]^a[j]: i in [1..#xbx]]: j in [1..#a]] cat [[xbx[i]^(Inverse(a[j])): i in [1..#xbx]]: j in [1..#a]];

      LenVector := [mylen(axbxa[j]): j in [1..#axbxa]];
      
      min, mink :=  MyMin(LenVector);
      
      return min, mink;
   catch e
      print "Error in PeelOff: ",e`Object;
      return 0, 0;
   end try;
end function;


/*--------------------------------------------------------------------*/
function Test(N, n, m, gen_a, gen_b, gen_x, mygt, mylen)
/**********************************************************************
   Sets up an environment for length based attack and then peels off one generator, 
   checks if the correct generator is at the t-th place
   
   returns  1 if yes
            0 if no
            -1 if error
**********************************************************************/
   try
      //Generate Braid Group
      B_N := BraidGroup(N);
      //Generate a_i and b_i
      a := GenRandElts(gen_a, B_N);
      b := GenRandElts(gen_b, B_N);
      //Generate x
      x_gen_nums := [Random([x: x in [-#a..#a]| x ne 0]):j in [1..gen_x]];
      x_gens := [];

      for i in x_gen_nums do
         if i lt 0 then
            Append(~x_gens, Inverse(a[-i]));
         else
            Append(~x_gens, a[i]);
         end if;
      end for;
      x := &*x_gens;
      xbx := [b[i]^x: i in [1..#b]];    

      min, mink := PeelOff(xbx, a, mygt, mylen);

      if x_gen_nums[#x_gen_nums] lt 0 then //take last of the x_gens since x^-1 b x 
         place := -x_gen_nums[#x_gen_nums];
      else
         place := x_gen_nums[#x_gen_nums]+n;
      end if;

      if place eq mink then
         return 1;
      else
         return 0;
      end if;
   catch e
      print "error in Test";
      error e`Object;
      return -1;
   end try;
   
end function;


/*--------------------------------------------------------------------*/
function NumSuccess    (T:
			t := 9,
                        N := 81,
                        n := 20,
                        m := 20,
                        gen_x := 30,
                        mylen := LenCan,
                        mygt := gtAV)
/**********************************************************************
runtime measuring of 4_2_avbased and its variations
**********************************************************************/
   print "Uses N=",N,", n=",n,", m=",m,", gen_x=",gen_x," len=LenCan, mygt=gtAV";
   k := 1;
   t0 := 0;

   print "First Test: Constant length";
   for length_a in [5,15,20] do
	gen_a := [length_a: i in [1..n]];
	for length_b in [5,15,20] do
		gen_b := [length_b: i in [1..m]];
		print "|ai| = ", gen_a, "|bi| = ", gen_b;
		i := 1;
		results := [];
		while i le T do
			t0 := Realtime();
			Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mygt, mylen));
			print "time = ",Realtime(t0);
			i +:= 1;
		end while;
	end for;
   end for;

   print "Second Test: Random length choose between 5 and 50";
   k := 1;
   while k lt t do	
	gen_a := [Random([5..50]) : i in [1..n]];
        gen_b := [Random([5..50]) : i in [1..m]];
	print "|ai|= ",gen_a," |bi|= ",gen_b;
 	i:=1;
	results := [];
 	while i le T do
		t0 := Realtime();
 		Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mygt, mylen));
		print "time = ", Realtime(t0);
 		i +:= 1;
   	end while;  
	k := k+1; 
   end while;   

   print "Third Test: high differences";
   k := 1;
   while k lt t do	
	gen_a := [Random([5,50,100]) : i in [1..n]];
        gen_b := [Random([5,50,100]) : i in [1..m]];
	print "|ai|= ",gen_a," |bi|= ",gen_b;
 	i:=1;
	results := [];
 	while i le T do
		t0 := Realtime();
 		Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mygt, mylen));
	 	print "time = ", Realtime(t0);	
		i +:= 1;
   	end while;
	k := k+1;   
   end while;   

   return 1;
end function;
