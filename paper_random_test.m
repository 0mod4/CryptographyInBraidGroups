function RandomTest(N, k, T) 
/**********************************************************************
   Generates T pseudorandom elements of B_N as product of k randomly 
   chosen artin generators and then counts how often each element 
   occurs
**********************************************************************/
   B_N := BraidGroup(N);
   SetSeed(Truncate(Realtime()));
   SetElementPrintFormat(~B_N, "Word");
   arr := AssociativeArray();
   for i in [1..T] do
      rand := NormalForm(&*[Random([B_N.j :j in [1..N-1]]): i in [1..k]]);
      try
         arr[WordToSequence(rand)] +:= 1;
      catch e
         arr[WordToSequence(rand)] := 1;
      end try;
   end for;
   
   //printing
   for i in Keys(arr) do
      print i,";",arr[i];
   end for;
   return 1;
end function;
