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

/*--------------------------------------------------------------------*/

function PeelOff(xbx, a, mygt, mylen)
/**********************************************************************
   Peels off the generator ai from xbx
   returns  the place of the resulting vector
**********************************************************************/
   try
      //run through all ai and inverses
      
      vector := [];
      
      for i in [1..#a] do //peel off a[i] means conjugating with a[i]^-1
         Append(~vector, [xbx[j]^Inverse(a[i]): j in [1..#xbx]]);
         Append(~vector, [xbx[j]^a[i]: j in [1..#xbx]]);
      end for;
      //print "vector erstellt";
      p := Id (Sym(#vector));
      MySort(~vector, mygt, mylen, ~p);

      return p;
   catch e
      print "Error in PeelOff: ",e`Object;
      return Id (Sym(#vector));
   end try;
end function;


/*--------------------------------------------------------------------*/
function Test(t, N, n, m, gen_a, gen_b, gen_x, mygt, mylen)
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

      perm := PeelOff(xbx, a, mygt, mylen);

      if x_gen_nums[#x_gen_nums] lt 0 then //take last of the x_gens since x^-1 b x 
         place := -x_gen_nums[#x_gen_nums]*2;
      else
         place := x_gen_nums[#x_gen_nums]*2 -1;
      end if;

      if place^(perm^-1) eq t then
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
                        t := 1,
                        N := 81,
                        n := 20,
                        m := 20,
                        gen_x := 30,
                        gen_a := [10 : i in [1..n]], 
                        gen_b := [10 : i in [1..m]],
                        mylen := LenGarside,
                        mygt := gtAV)
/**********************************************************************
   Coordinates T Runs of 'Peeling off' one generator, returns howmany times
   the correct generator gave a vector at t-th place (ordert according to length)
   N:       Working in BraidGroup B_N (N-1 Artin generators)
   n:       number of a_i, Generators of the subgroup where conjugator x is from
   m:       number of b_i, Elements that are conjugated with x
   gen_a:   number of generators for the a_i
   gen_b:   number of generators for the b_i
   gen_x:   number of a_i that generate x
   mylen:   length function for the elements in B_N
   mygt:    linear ordering on vectors containing elements from B_N using len, here average based is implemented
**********************************************************************/
   i:=1;
   results := [];
   while i le T do
      Append(~results, Test(t, N, n, m, gen_a, gen_b, gen_x, mygt, mylen));
      print "Test Nr. ",i," finished";
      i +:= 1;
   end while;
   
   //Get Results
   testsuccess  := #[results[i]:i in [1..#results]|results[i] eq 1];
   testfail     := #[results[i]:i in [1..#results]|results[i] eq 0];
   testerror    := #[results[i]:i in [1..#results]|results[i] eq -1];

   print T," runs: ",testsuccess,"x success, ",testfail,"x fail, ",testerror," Errors encountered";
   
   return testsuccess;
end function;
