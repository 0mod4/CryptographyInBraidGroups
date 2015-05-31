load "len.m";
load "mygt.m";

/*--------------------------------------------------------------------*/
function GenRandElts(gen_num, B_N) 
/**********************************************************************
   Generates Sequence of Random elements in Braidgroup Bn of length l
   where Entry i consists of gen_num[i] generators
**********************************************************************/
   a := [Random(B_N, gen_num[i]): i in [1..#gen_num]];
   return a;
end function;

function GenX(a, gen_num)
/**********************************************************************
   Generates x out of Random a_i in Braidgroup Bn
**********************************************************************/
   x := &*[Random(a):j in [1..gen_num]];
   return x;
end function;

/*--------------------------------------------------------------------*/

function PeelOff(a, xbx, my_len, my_gt, B_N)
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
         if (my_gt(axbxa, test, my_len) eq 1) then
            axbxa := test;
            min_a := a[i];
            min_a_num := i;
         end if;
         
         test := [xbx[j]^Inverse(a[i]): j in [1..#xbx]];
         if (my_gt(axbxa, test, my_len) eq 1) then
            axbxa := test;
            min_a := Inverse(a[i]);
            min_a_num := -i;
         end if;
      end for;
      return min_a, min_a_num, axbxa;
   catch e
      print "Error in PeelOff: ",e`Object;
      return Identity(B_N), 1, [Identity(B_N)];
   end try;
end function;

function Attack(a, b, xbx, my_len, my_gt, gen_x, B_N)
/**********************************************************************
   Length based attack on xbx
   returns true if successfull, else false
***********************************************************************/
   try
      generators := [];
      for i in [1..gen_x] do
         min_a, min_a_num, axbxa := PeelOff(a, xbx, my_len, my_gt, B_N);
         if (my_gt(axbxa, xbx, my_len) eq 1) then
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
function Test(N, n, m, gen_a, gen_b, gen_x, my_gt)
/**********************************************************************
   Sets up an environment for length based attack and then attacks
   
   returns  1 if was successfull
            0 if was not successfull
            -1 if error
**********************************************************************/
   try
      //Generate Braid Group
      B_N := BraidGroup(N);
      //Generate a_i and b_i
      a := GenRandElts(gen_a, B_N);
      b := GenRandElts(gen_b, B_N);
      //Generate x
      x := &*[Random(a):j in [1..gen_x]];
      xbx := [b[i]^x: i in [1..#b]];
      
      result_can := Attack(a, b, xbx, LenCan, my_gt, gen_x, B_N);
      result_gar := Attack(a, b, xbx, LenGarside, my_gt, gen_x, B_N);
      return [result_can, result_gar];
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
                        gen_x := 30,
                        my_gt := gtAV, 
                        gen_a := [10 : i in [1..n]], gen_b := [10 : i in [1..m]])
/**********************************************************************
   Coordinates T Runs of the attack, returns howmany have been successfull, compares different length
function and their success
   N:       Working in BraidGroup B_N (N-1 Artin generators) filled from List
   n:       number of a_i, Generators of the subgroup where conjugator x is from
   m:       number of b_i, Elements that are conjugated with x
   gen_a:   number of generators for the a_i
   gen_b:   number of generators for the b_i
   gen_x:   number of a_i that generate x
   len:     length function for the elements in B_N
   gt:      linear ordering on vectors containing elements from B_N using len
**********************************************************************/
   SetSeed(Truncate(Realtime()));

   print "n=",n," m=",m," LenCan vs LenGarside, gtAV, gen_a=gen_b=10";
   print "In Probability of success test for the whole attack largest difference between AVGar and AVCan: gen_x=8, N=11";
   gen_x := 8;
   N := 11;
   i:=1;
   results := [];
   while i le T do
	result := Test(N, n, m, gen_a, gen_b, gen_x, my_gt);
   	Append(~results, result); //resultCan, resultGar
	print "Test ",i," finished: ",result;
	i := i+1;
   end while;
   
   print "N = ",N," gen_x=",gen_x," gen_a=10, gen_b=10";

   //Get Results
   testsequalsuc	:= #[results[i]:i in [1..#results]|results[i][1] eq results[i][2] and results[i][2] eq 1];
   testsequalfail	:= #[results[i]:i in [1..#results]|results[i][1] eq results[i][2] and results[i][2] eq 0];
   testsOnlyGar 	:= #[results[i]:i in [1..#results]|results[i][1] ne results[i][2] and results[i][2] eq 1];
   testsOnlyCan 	:= #[results[i]:i in [1..#results]|results[i][1] ne results[i][2] and results[i][1] eq 1];
   
   print T," runs: ",testsequalsuc,"x both successfull, ",testsequalfail,"x both fail, ",testsOnlyGar,"x Gar successfull, Can not, ",testsOnlyCan," x Can successfull, Gar not";

   return 1;
end function;
