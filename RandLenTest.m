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

function lenW (a, w, N)
   lenGarA := [];
 
   for j in [1..#a] do
   	nf := CanonicalFactorRepresentation(NormalForm(a[j]));
   	r := nf[2]; //Delta exponent
   	p := nf[3];//positive braids as permutations
   
   	lenGar := -r * Binomial(N,2);
   	for i in [1..#p] do
   		pi := Eltseq(p[i]);
        	lenP := #[<i,j> : j in [i..N], i in [1..N] | pi[i] gt pi[j]];
		if IsDefined(w, i) then
			lenP := lenP * w[i];
		end if;
        	lenGar := lenGar + lenP;
   	end for;
	lenGarA[j] := lenGar;
   end for;
   return lenGarA;
end function;

function gtAVw (a,b,w,N)
   if(&+lenW(a, w, N) gt &+lenW(b, w, N)) then
      return 1;
   elif (&+lenW(a, w, N) lt &+lenW(b, w, N)) then
      return -1;
   else
      return 0;
   end if;
end function;

function MySortW(Q, w, N, permut)
   n := #Q;
   if n eq 0 then
      permut := Id(Sym(1));
      return Q, permut;
   end if;

   inc := [1];
   h := 1;
   repeat
      h := 3*h + 1;
      Append(~inc, h);
   until h ge n;

   perm := [1 .. n];
   t := Max(#inc - 2, 1);
   for s := t to 1 by -1 do
      h := inc[s];
      for j := h + 1 to n do
         i := j - h;
         k := Q[j];
         kpos := perm[j];
         
         while gtAVw(k, Q[i], w, N) lt 0 do
         	Q[i+h] := Q[i];
         	perm[i+h] := perm[i];
         	i -:= h;
         	if i le 0 then
             		break;
         	end if;
         end while;
         
         Q[i+h] := k;
         perm[i+h] := kpos;
      end for;
   end for;
   permut := Sym(n) ! perm;

   return Q, permut;
end function;

/*--------------------------------------------------------------------*/

function PeelOff(xbx, a, w, N)
/**********************************************************************
   Peels off the generator ai from xbx
   returns  the place of the resulting vector in 
**********************************************************************/
   try
      //run through all ai and inverses
      
      vector := [];
      
      for i in [1..#a] do //peel off a[i] means conjugating with a[i]^-1
         Append(~vector, [xbx[j]^Inverse(a[i]): j in [1..#xbx]]);
         Append(~vector, [xbx[j]^a[i]: j in [1..#xbx]]);
      end for;
      //print "vector erstellt";
      p1 := Id (Sym(#vector));
      fake_w := AssociativeArray();
      print "Start Sort1 using Gar_V";
      vector, p1 := MySortW(vector, fake_w, N, p1);
      print "End Sort1";  

      p2 := Id (Sym(#vector));
      print "Start Sort2";
      vector, p2 := MySortW(vector, w, N, p2);
      print "End Sort2";
      //Sort(~vector, wrapGT, ~p);
      //print "sorted. Permutation: ";
      //print [1..#vector]^(p^-1);

      return p1, p2;
   catch e
      print "Error in PeelOff: ",e`Object;
      return Id (Sym(#vector));
   end try;
end function;


/*--------------------------------------------------------------------*/
function Test(t, N, n, m, gen_a, gen_b, gen_x, w)
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

      print "Test setup done";      

      perm1, perm2 := PeelOff(xbx, a, w, N);

      if x_gen_nums[#x_gen_nums] lt 0 then //take last of the x_gens since x^-1 b x 
         place := -x_gen_nums[#x_gen_nums]*2;
      else
         place := x_gen_nums[#x_gen_nums]*2 -1;
      end if;
      //print "got perm, check for ", place, " is at place: ";
      //print place^(perm^-1);
      suc1 := -1;
      suc2 := -1;

      if place^(perm1^-1) eq t then
         suc1 := 1;
      else
         suc1 := 0;
      end if;

      if place^(perm2^-1) eq t then
         suc2 := 1;
      else
         suc2 := 0;
      end if;
      return suc1, suc2;
   catch e
      print "error in Test";
      error e`Object;
      return -1, -1;
   end try;
   
end function;


/*--------------------------------------------------------------------*/
function NumSuccess    (T:
                        t := 1,
                        N := 9,
                        n := 20,
                        m := 20,
                        gen_x := 9,
                        gen_a := [10 : i in [1..n]], 
                        gen_b := [10 : i in [1..m]]
                        )
/**********************************************************************
   Coordinates T Runs of 'Peeling off' one generator, returns howmany times
   the correct generator gave a vector at t-th place (ordert according to length)
   N:       Working in BraidGroup B_N (N-1 Artin generators)
   n:       number of a_i, Generators of the subgroup where conjugator x is from
   m:       number of b_i, Elements that are conjugated with x
   gen_a:   number of generators for the a_i
   gen_b:   number of generators for the b_i
   gen_x:   number of a_i that generate x
   length function is chosen as weighted Garside length
   ordering function is set to average based linear ordering
**********************************************************************/

   //Choose weights for weighted Garside length
   max_i := 100;
   min_wi := 1;
   max_wi := 10;
   num_wi := 60;
   w := AssociativeArray();
   print "weights: ";
   for k in [1..num_wi] do
        i := Random([1..max_i]);
        wi := Random([wi: wi in [min_wi..max_wi]| wi ne 0 and wi ne 1]);
	w[i] := wi;
        print "w[",i,"] = ",wi;
   end for;   
   results := [];
   for i in [1..T] do
      suc1, suc2 := Test(t, N, n, m, gen_a, gen_b, gen_x, w); //suc1: old Gar, suc2: weighted Gar
      Append(~results, [suc1, suc2]);
      print "Test Nr. ",i," finished";
   end for;
   
   //Get Results
   x1eqGar  := #[results[i]:i in [1..#results]|results[i][1] eq 1 and results[i][2] eq 1];
   x0eqGar  := #[results[i]:i in [1..#results]|results[i][1] eq 0 and results[i][2] eq 0];
   x1neGar  := #[results[i]:i in [1..#results]|results[i][1] eq 1 and results[i][2] eq 0];
   x0neGar  := #[results[i]:i in [1..#results]|results[i][1] eq 0 and results[i][2] eq 1];
   testerror    := #[results[i]:i in [1..#results]|results[i][1] eq -1 or results[i][2] eq -1];

   print T," runs: \n",x1eqGar,"x x=1=Gar, ";
   print x0eqGar,"x x=0=Gar";
   print x1neGar,"x x=1!=Gar";
   print x0neGar,"x x=0!=Gar";   
   return 1;
end function;
