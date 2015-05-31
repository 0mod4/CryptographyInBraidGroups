load "len.m";

function PeelOff(xbx, a, mylen)
/**********************************************************************
   Peels off the generator ai from xbx
   returns  the place of the resulting vector in 
**********************************************************************/

      axbxa := [[xbx[i]^a[j]: i in [1..#xbx]]: j in [1..#a]] cat [[xbx[i]^Inverse(a[j]): i in [1..#xbx]]: j in [1..#a]];    

      lenVector := [mylen(axbxa[j]): j in [1..#axbxa]];      
      mu := [Minimum([lenVector[i][j]: i in [1..#lenVector]]): j in [1..#lenVector[1]]];
      averageBasedLengths := AssociativeArray();
      majorityBasedLengths := AssociativeArray();
      
      for i in [1..#lenVector] do
         averageBasedLength := &+lenVector[i];
         majorityBasedLength := #[1: j in [1..#lenVector[i]]|lenVector[i][j] eq mu[j]];
         if IsDefined(averageBasedLengths, averageBasedLength+1) then
	 	averageBasedLengths[averageBasedLength+1] +:= 1;
         else
   		averageBasedLengths[averageBasedLength+1] := 1;
	 end if;
         
         if IsDefined(majorityBasedLengths, majorityBasedLength+1) then
	 	majorityBasedLengths[majorityBasedLength+1] +:= 1;
         else
   		majorityBasedLengths[majorityBasedLength+1] := 1;
	 end if;
      end for;
      
      return averageBasedLengths, majorityBasedLengths;
end function;


/*--------------------------------------------------------------------*/
function Test(t, N, n, m, gen_a, gen_b, gen_x, mylen)
/**********************************************************************
   Sets up an environment for length based attack and then peels off one generator, 
   returns length
**********************************************************************/
      //Generate Braid Group
      B_N := BraidGroup(N);
      //Generate a_i and b_i
      a := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_a[j]]]: j in [1..n]];
      b := [&*[Random([B_N.k: k in [1..N-1]]): i in [1..gen_b[j]]]: j in [1..m]];
      //Generate x
      x := &*[Random(a)^Random([-1,1]): j in [1..gen_x]];
      xbx := [b[i]^x: i in [1..#b]];
      
      len1, len2 := PeelOff(xbx, a, mylen);
      return len1, len2;  
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
   Variation of VecLenTest.m using old Random
**********************************************************************/
print "Uses N=",N," n=",n," m=",m," gen_x=",gen_x," gen_a=",gen_a," gen_b=",gen_b," len=LenCan";
   mylen := LenCan;
   i:=1;
   while i le T do
      SetSeed(Truncate(Realtime()));

      average, majority := Test(t, N, n, m, gen_a, gen_b, gen_x, mylen);
      
      print  "average sizes: ";
      for k in Keys(average) do
         print "asize ",k-1,": ",average[k];
      end for;
      
      print  "majority sizes: ";
      for k in Keys(majority) do
         print "msize ",k-1,": ",majority[k];
      end for;
      
      i +:= 1;
   end while;

print "Uses N=",N," n=",n," m=",m," gen_x=",gen_x," gen_a=",gen_a," gen_b=",gen_b," len=LenGarside";
   mylen := LenGarside;
   i:=1;
   while i le T do
      SetSeed(Truncate(Realtime()));

      average, majority := Test(t, N, n, m, gen_a, gen_b, gen_x, mylen);

      print  "average sizes: ";
      for k in Keys(average) do
         print "asize ",k-1,": ",average[k];
      end for;

      print  "majority sizes: ";
      for k in Keys(majority) do
         print "msize ",k-1,": ",majority[k];
      end for;

      i +:= 1;
   end while;

   
   return 1;
end function;
