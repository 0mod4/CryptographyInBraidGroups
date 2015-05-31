load "len.m";
load "mygt.m";


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

function PeelOff(a, xbx, mylen, B_N)
/**********************************************************************
   Tries to peel off one generator from xbx by searching the ai in a 
   or an inverse such that axbxa has minimal length.
   returns  the group element ai for which it is minimal,
            the number i of this element (-i if inverse)
            the element after 'peeling off' axbxa
**********************************************************************/
   try
	axbxa := [[xbx[i]^a[j]: i in [1..#xbx]]: j in [1..#a]] cat [[xbx[i]^(Inverse(a[j])): i in [1..#xbx]]: j in [1..#a]];
	LenVector := [mylen(axbxa[j]):j in [1..#axbxa]];
	min_a, min_a_num := MyMin(LenVector);
    
	return min_a, min_a_num, axbxa[min_a_num];
   catch e
	print "Error in PeelOff: ",e`Object;
	return Identity(B_N), 1, [Identity(B_N)];
   end try;
end function;

function Attack(a, b, xbx, mylen, gen_x, B_N)
/**********************************************************************
   Length based attack on xbx
   returns true if successfull, else false
***********************************************************************/
   try
      generators := [];
      for i in [1..gen_x] do
         min_a, min_a_num, axbxa := PeelOff(a, xbx, mylen, B_N);
	 //get majority based length of axbxa and xbx
	 len_axbxa := mylen(axbxa);
	 len_xbx := mylen(xbx);
	 ord_axbxa := #[1: i in [1..#axbxa]|len_axbxa[i] le len_xbx[i]];
	 ord_xbx := #[1: i in [1..#xbx]|len_xbx[i] le len_axbxa[i]];
         if (ord_axbxa lt ord_xbx) then
            return 0;
         else
            xbx := axbxa;
            Append(~generators, min_a_num);
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
function Test(N, n, m, gen_a, gen_b, gen_x, my_len)
/**********************************************************************
   Sets up an environment for length based attack and then attacks
   
   returns  1 if was successfull
            0 if was not successfull
            -1 if error
**********************************************************************/
   try
      //Generate Braid Group
      B_N := BraidGroup(N);
      a := [Random(B_N, gen_a[j]): j in [1..n]];
      b := [Random(B_N, gen_b[j]): j in [1..m]];
      //Generate x
      x := &*[Random(a)^Random([-1,1]): j in [1..gen_x]];
      xbx := [b[i]^x: i in [1..#b]];

      return Attack(a, b, xbx, my_len, gen_x, B_N);
   catch e
      print "error in Test";
      error e`Object;
      return -1;
   end try;
   
end function;


/*--------------------------------------------------------------------*/
function NumSuccess    (T:
            //            N := 4,
                        n := 20,
                        m := 20,
                        gen_x := 2,
                        gen_a := [10 : i in [1..n]], gen_b := [10 : i in [1..m]])
/**********************************************************************
   Variation of 4_4.m using majority based linear ordering
**********************************************************************/
   print "n=",n," m=",m," LenGarside, majbased, gen_a=",gen_a,"gen_b=",gen_b;
   mylen := LenGarside;
   for gen_x in [2..18] do
	for N in [4..20] do
   		i:=1;
	   	results := [];
	  	while i le T do
      			Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
      			i +:= 1;
   		end while;
   
  		print "N = ",N," gen_x=",gen_x;

   		testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
   		testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
   		testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

   		print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
  	 end for;
   end for;

   print "n=",n," m=",m," LenGarside, majbased, gen_a=",gen_a, "gen_b=",gen_b;
   mylen := LenGarside;
   for gen_x in [10..18] do
        for N in [4..20] do
                i:=1;
                results := [];
                while i le T do
                        Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
                        i +:= 1;
                end while;
   
                print "N = ",N," gen_x=",gen_x;

                testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
                testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
                testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

                print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
         end for;
   end for;


   return 1;
end function;
