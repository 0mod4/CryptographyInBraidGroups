function MyMin(a, compare, len)
//compare(a,b,len) returns 1 if a>b -1 if a<b else 0
	n:= #a;
	if n eq 0 then
		return 0;
	elif n eq 1 then
		return a[1];
	end if;
        
	min := a[1];
        mink := 1;
	for k in [2..#a] do
		if compare(min, a[k], len) eq 1 then
			min := a[k];
                        mink := k;
		end if;
	end for;

	return min, mink;
end function;
