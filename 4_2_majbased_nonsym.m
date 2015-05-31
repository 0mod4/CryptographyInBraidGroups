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
		//print LenVector;
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
   Code for random generation in 'the old' and 'the new way' is available
   
   returns  1 if yes
            0 if no
            -1 if error
**********************************************************************/
   try
      //Generate Braid Group
      B_N := BraidGroup(N);
      //Generate a_i and b_i
     // a := [Random(B_N, gen_a[j]): j in [1..n]];
     // b := [Random(B_N, gen_b[j]): j in [1..m]];
      //Generate the old way: 
      a := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_a[j]]]: j in [1..n]];
      b := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_b[j]]]: j in [1..m]];
      //Generate x
      x_gen_nums := [Random([x: x in [-#a..#a]| x ne 0]):j in [1..gen_x]];
      //print "x generators: ", x_gen_nums;
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
                        N := 81,
                        n := 20,
                        m := 20,
                        gen_x := 30)
/**********************************************************************
   Variation of 4_2_majbased with gen_a and gen_b chosen from given sets
   for Garside and Canonical length
**********************************************************************/
   print "Uses N=",N,", n=",n,", m=",m," gen_x=",gen_x," len=LenGarside, majbased";
   mylen := LenGarside;
   for length_ai in [5,10,15,20,25] do
	gen_a := [length_ai: i in [1..n]];
	for length_bi in [5,10,15,20,25] do
		gen_b := [length_bi: i in [1..m]];
		print "|ai| = ",length_ai," |bi| = ",length_bi;
   		i:=1;
   		results := [];
   		while i le T do
      			Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
      			i +:= 1;
   		end while;
   		
		testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
	   	testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
   		testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

   		print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
	end for;
   end for;

   print "Uses N=",N,", n=",n,", m=",m," len=LenCan, majbased";
   mylen := LenCan;
   for length_ai in [5,10,15,20,25] do
        gen_a := [length_ai: i in [1..n]];
        for length_bi in [5,10,15,20,25] do
                gen_b := [length_bi: i in [1..m]];
                print "|ai| = ",length_ai," |bi| = ",length_bi;
                i:=1;
                results := [];
                while i le T do
                        Append(~results, Test(N, n, m, gen_a, gen_b, gen_x, mylen));
                        i +:= 1;
                end while; 
                
                testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
                testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
                testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];
    
                print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered"; 
        end for;
   end for;

   return 1;
end function;
