/*--------------------------------------------------------------------*/
function gtAV (a,b,len)
/**********************************************************************
   linear ordering on tupels of conjugates: Average based from Script
**********************************************************************/
   if(&+len(a) gt &+len(b)) then
      return 1;
   elif (&+len(a) lt &+len(b)) then
      return -1;
   else 
      return 0;
   end if;
end function;
