load "len.m";


function MyMin(LenVec)
    n:= #LenVec;
    if n eq 0 then
        return 0;
    elif n eq 1 then
        return LenVec[1];
    end if;
    m := #LenVec[1];

    mu := [Minimum([LenVec[i][j]: i in [1..#LenVec]]): j in [1..m]];

    majlen := [#[1: i in [1..#LenVec[j]]|LenVec[j][i] eq mu[i]]:j in [1..#LenVec]];

    min, mink := Minimum([m-majlen[i]: i in [1..#majlen]]);
    return m-min, mink;

end function;

/*--------------------------------------------------------------------*/

function PeelOff(xbx, a, mylen)
	/**********************************************************************
	  Peels off the generator ai from xbx
	  returns  the place of the resulting vector in 
	 **********************************************************************/
	try
                axbxa := [[xbx[i]^a[j]: i in [1..#xbx]]: j in [1..#a]] cat [[xbx[i]^(Inverse(a[j])): i in [1..#xbx]]: j in [1..#a]];

                LenVector := [mylen(axbxa[j]): j in [1..#axbxa]];

                min, mink :=  MyMin(LenVector);

                return min, mink;
	catch e
		print "Error in PeelOff: ",e`Object;
		return 0,0;
	end try;
end function;


/*--------------------------------------------------------------------*/
function Test(N, n, m, gen_a, gen_b, gen_x, mylen)
/**********************************************************************
   Sets up an environment for length based attack and then peels off one generator, 
   checks if the correct generator is at the 1st place
   
   returns  1 if yes
            0 if no
            -1 if error
**********************************************************************/
   try
      B_N := BraidGroup(N);

      a := [Random(B_N, gen_a[j]): j in [1..n]];
      b := [Random(B_N, gen_b[j]): j in [1..m]];

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

      min, mink := PeelOff(xbx, a, mylen);

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
			t := 10,
                        N := 81,
                        n := 20,
                        m := 20,
                        gen_x := 30)
/**********************************************************************
   variation of 4_2_majbased with gen_a, gen_b chosen with high difference
**********************************************************************/
   print "Uses N=",N,", n=",n,", m=",m," len=LenCan, majbased t=10";
   mylen := LenCan;
   k := 1;
   while k le t do
	gen_a := [Random([5,50,100]): i in [1..n]];
	gen_b := [Random([5,50,100]): i in [1..m]];
	print "|ai| = ",gen_a," |bi| = ",gen_b;
	i:=1;
	results := [];
	while i le T do
		Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
		i +:= 1;
	end while;
	
	k := k+1;
   
	testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
  	testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
	testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

	print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
   end while;

   print "Uses N=",N,", n=",n,", m=",m," len=LenGarside, majbased t=10";
   mylen := LenGarside;
   k := 1;
   while k le t do
        gen_a := [Random([5,50,100]): i in [1..n]];
        gen_b := [Random([5,50,100]): i in [1..m]];
        print "|ai| = ",gen_a," |bi| = ",gen_b;
        i:=1;
        results := [];
        while i le T do
                Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
                i +:= 1;
        end while;

        k := k+1;

        testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
        testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
        testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

        print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
   end while;


   return 1;
end function;
