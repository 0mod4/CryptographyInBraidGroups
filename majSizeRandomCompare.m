load "len.m";

/*--------------------------------------------------------------------*/
function Test(t, N, n, m, gen_a, gen_b, gen_x, mylen)
/**********************************************************************
   Sets up an environment for length based attack and then peels off one generator, 
   checks if the correct generator is at the 1st place.
   Uses 'old' and 'new' Random in comparison
   
   returns  1 if yes
            0 if no
            -1 if error
**********************************************************************/
      //Generate Braid Group
      B_N := BraidGroup(N);
      SetElementPrintFormat(~B_N, "Word");
      //Generate a_i and b_i
      a := [Random(B_N, gen_a[j]): j in [1..n]];
      b := [Random(B_N, gen_b[j]): j in [1..m]];
      oa := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_a[j]]]: j in [1..n]];
      ob := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_b[j]]]: j in [1..m]];
      //Generate x
      x := &*[Random(a)^Random([-1,1]): j in [1..gen_x]];
      ox := &*[Random(a)^Random([-1,1]): j in [1..gen_x]];
 
      xbx := [b[i]^x: i in [1..#b]];
      oxbx := [ob[i]^ox: i in [1..#ob]];

      axbxa := [[xbx[i]^a[j]: i in [1..#xbx]]: j in [1..#a]] cat [[xbx[i]^Inverse(a[j]): i in [1..#xbx]]: j in [1..#a]];
      oaxbxa := [[oxbx[i]^oa[j]: i in [1..#oxbx]]: j in [1..#oa]] cat [[oxbx[i]^Inverse(oa[j]): i in [1..#oxbx]]: j in [1..#oa]];
      LenVec := [mylen(axbxa[j]): j in [1..#axbxa]];
      mu := [Minimum([LenVec[i][j]: i in [1..#LenVec]]): j in [1..m]];
      len := [#[1: i in [1..#LenVec[j]]|LenVec[j][i] eq mu[i]]:j in [1..#LenVec]];

      oLenVec := [mylen(oaxbxa[j]): j in [1..#oaxbxa]];
      omu := [Minimum([oLenVec[i][j]: i in [1..#oLenVec]]): j in [1..n]];
      olen := [#[1: i in [1..#oLenVec[j]]|oLenVec[j][i] eq omu[i]]:j in [1..#oLenVec]];

      return len, olen;

end function;


/*--------------------------------------------------------------------*/
function VecLenTest    (T:
                        t := 1,
                        N := 81,
                        n := 20,
                        m := 20,
                        gen_x := 30,
                        gen_a := [10 : i in [1..n]], 
                        gen_b := [10 : i in [1..m]])
/**********************************************************************
   Variation of 4_2.m using majority based linear ordering
   Coordinates T Runs of 'Peeling off' one generator comparing different random methods
   and uniform random method, prints the lengths of the vectors, 
   N:       Working in BraidGroup B_N (N-1 Artin generators)
   n:       number of a_i, Generators of the subgroup where conjugator x is from
   m:       number of b_i, Elements that are conjugated with x
   gen_a:   number of generators for the a_i
   gen_b:   number of generators for the b_i
   gen_x:   number of a_i that generate x
   mylen:   length function for the elements in B_N
**********************************************************************/
print "Uses N=",N," n=",n," m=",m," gen_x=",gen_x," gen_a=",gen_a," gen_b=",gen_b," len=LenGarside"; 
  mylen := LenGarside;
  i:=1;
   while i le T do
      SetSeed(Truncate(Realtime()));

      majority, omajority := Test(t, N, n, m, gen_a, gen_b, gen_x, mylen);
      
      print  "majority sizes old random: ";
      for l in [0..m] do
		print "mosize ",l,": ",#[v:v in [1..#omajority]|omajority[v] eq l];
      end for;
      
      print  "majority sizes new random: ";
      for l in [0..m] do
      		print "msize ",l,": ",#[v:v in [1..#majority]|majority[v] eq l];
      end for;
      
      i +:= 1;
   end while;

print "Uses N=",N," n=",n," m=",m," gen_x=",gen_x," gen_a=",gen_a," gen_b=",gen_b," len=LenCan";
  mylen := LenCan;
  i:=1;
   while i le T do
      SetSeed(Truncate(Realtime()));

      majority, omajority := Test(t, N, n, m, gen_a, gen_b, gen_x, mylen);
      
      print  "majority sizes old random: ";
      for l in [0..m] do
                print "mosize ",l,": ",#[v:v in [1..#omajority]|omajority[v] eq l];
      end for;

      print  "majority sizes new random: ";
      for l in [0..m] do
                print "msize ",l,": ",#[v:v in [1..#majority]|majority[v] eq l];
      end for;
      
      i +:= 1;
   end while;

   
   return 1;
end function;
