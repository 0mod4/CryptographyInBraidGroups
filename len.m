/*********************************************************************
Different length functions helper
**********************************************************************/

function LenCan (a)
/**********************************************************************
   Length Function for Elements in B_N: Canonical Length
**********************************************************************/
   l:=[CanonicalLength(a[i]):i in [1..#a]];
   return l;
end function;

function LenSup (a)
/**********************************************************************
Length Function for Elements in B_N: Supremum
**********************************************************************/
l:=[Supremum(a[i]):i in [1..#a]];
return l;
end function;

function LenInf (a)
/**********************************************************************
Length Function for Elements in B_N: Supremum
**********************************************************************/
l:=[Infimum(a[i]):i in [1..#a]];
return l;
end function;

function LenTest (a)
/**********************************************************************
   Length Function for Elements in B_N: Test: gives constant back
**********************************************************************/
   return 5;
end function;

function LenGarside(a)
/**********************************************************************
   Length Function for Elements in B_N: gives the Garside length (see paper)
**********************************************************************/
   l:=[#NormalForm(a[i]):i in [1..#a]];
   return l;
end function;
