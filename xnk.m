/****************************************************************/
/* Helper classes for randBraid.m for details see paper		*/
/* 'Generating random Braids'					*/
/****************************************************************/

function wo (a, b, B_n)
/****************************************************************/
/* returns a\b as defined in the paper where a and b are members*/
/* of B_n as is the result and a flag if an error occured	*/
/****************************************************************/
	a_seq := WordToSequence(a);
	b_seq := WordToSequence(b);
	//B_n := BraidGroup(n);
	if #a_seq gt 1 then
		print "Format of a not supported";
		return B_n.1, 1;
	end if;
	j := a_seq[1];
	if j gt n-1 then
		print "too large values in a";
		return B_n.1, 1;
	end if;

//sj\si
	if #b_seq eq 1 then
		i := b_seq[1];
		if i gt n-1 then
			print "too large values in b";
			return B_n.1, 1;
		end if;
		if Abs(i-j) gt 1 then
			return B_n.i, 0;
		elif Abs(i-j) eq 1 then
			return B_n.i*B_n.j, 0;
		else
			return B_n.1, 0;
		end if;		
	end if;

//sj\(si-1 si)
	if (#b_seq eq 2) and (b_seq[1]+1 eq b_seq[2]) then
		i := b_seq[2];
		if i gt n-1 then
			print "too large values in b";
			return B_n.1, 1;
		end if;
		if j eq i-2 then
			return (B_n.(i-1)*B_n.i)*(B_n.(i-2)*B_n.(i-1)), 0;
		elif j eq i-1 then
			return B_n.i, 0;
		elif j eq i+1 then
			return B_n.(i-1)*B_n.i*B_n.(i+1), 0;
		else
			return B_n.(i-1)*B_n.i, 0;
		end if;
	end if;

//sj\(si si-1 ... sm)
	m := b_seq[#b_seq];
	for  i in [1..#b_seq-1] do
		if b_seq[i] ne b_seq[i+1]+1 then
			print "Format of b not supported";
			return B_n.1, 1;
		end if;
	end for;
	i := b_seq[1];
	if (i gt n-1) or (m gt n-1) then
		print "too large values in b";
		return B_n.1, 1;
	end if;
	if j eq m-1 then
		retval := B_n.i;
		k := i-1;
		while k ge m-1 do
			retval := retval*B_n.k;
			k := k-1;
		end while;
		return retval, 0;
	elif j eq i then
		retval := B_n.(i-1);
		k := i-2;
		while k ge m do
			retval := retval*B_n.k;
			k := k-1;
		end while;
		return retval, 0; 
	elif j eq i+1 then
		retval := B_n.i*B_n.(i+1);
		k := i-1;
		print k,",",m;
		while k ge m do
			retval := retval*(B_n.k*B_n.(k+1));
			k := k-1;
		end while;
		return retval, 0;
	else
		return b;
	end if;

//else
	print "Format not supported";
	return B_n.1, 1;
end function;


function Min (S)
/****************************************************************/
/* returns the minimal element in S consisting of braids	*/
/****************************************************************/
	min := S;
	for s in S do
		Exclude(~min, s);
		notmin := false;
		for c in min do
			if c le s then
				notmin := true;
				break;
			end if;
		end for;
		if notmin eq false then
			Include(~min, s);
		end if;
	end for;
	return min;
end function;


function Sab_eq_Test(F, a, b, a_, b_, n)
/****************************************************************/
/* returns {} containing the elements in Sa_b_\Sab 		*/
/* that is either nothing or sigmab, sigma a-1 or		*/ 
/*					       sigmab*..*sigma a*/
/****************************************************************/
B_n := Universe(F);
//Test: sigma b=skd
if (b-1 lt a) then
	if (a-1 gt 0) and (B_n.(a-1) in F) then
		return {B_n.(a-1)};
	elif (a gt 0) and (B_n.a in F) then
		return {B_n.a};
	end if;
else
	if (a_ eq a-1) and (b_ eq b) then
		if (a-1 gt 0) and (B_n.(a-1) in F) then
			return {B_n.(a-1)};
		elif (b lt n) and (B_n.b in F) then
			return {B_n.b};
		elif ( &*[B_n.(b-k) : k in [0..(b-a)]] in F ) then
			return {&*[B_n.(b-k) : k in [0..(b-a)]]};
		end if;
	elif (a_ eq a) and (b_ eq b+1) then
		if (b lt n) and (B_n.b in F) then
			return {B_n.b};
		elif ( &*[B_n.(b-k) : k in [0..(b-a)]] in F ) then
			return {&*[B_n.(b-k) : k in [0..(b-a)]]};
		end if;
	end if;
end if;

return {};

end function;


function xnk(n,k,w,m,xnj)

/****************************************************************/
/* returns the number xnk(w,m) (Algorithm 2 in paper)		*/
/****************************************************************/

if (n lt 2) or (k lt 0) or (m gt n-1) then
	print "invalid input";
	return -1;
end if;

B_n := BraidGroup(n);

if IsIdentity(w) then
	if m lt 1 then
		print "invalid input for w id";
		return -1;
	end if;
end if;

w_seq := WordToSequence(w);
if #w_seq eq 0 then
	j := 0;
else
	j := w_seq[#w_seq];
end if;

if m lt j-1 then
	print "invalid input for this w";
	return -1;
end if;

if #xnj lt k+1 then
	print "invalid input: xnj too short";
	return -1;
end if;

t := #w_seq;

if t gt k then
	return 0;
end if;

F := {};

for i in [1..t] do
	if B_n.w_seq[i] in F then
		return 0;
	else 
		S := {B_n.ai:ai in [1..w_seq[i]]};
		for b in F do
			awob, err := wo(B_n.w_seq[i], b, B_n);
			if err eq 0 then
				Include(~S,awob);
			else
				print "Error in wo(",B_n.w_seq[i],", ",b,", ",n,")";
			end if; 
		end for;
		F := Min(S);
	end if;
end for;

F := Min({B_n.k: k in [1..m]} join F);
//T[l][r][s] init with 0
Ts := [0: c in [0..n-1]];
Tr := [Ts: c in [0..n-1]];
T  := [Tr: c in [0..k-t]];
T_ := [Tr: c in [0..k-t]];

if t eq 0 then
	a := 1;
	b := 1;
	alpha := 0;
	T[1][1][1] := 1;
else
	a := w_seq[#w_seq];
	b := a;
	alpha := 0;
	T[1][1][1] := 1;
end if;
mycount := 0;
while not ((a eq 1) and (b eq n)) and (mycount lt 1) do
	//mycount := mycount+1;
	//check if sb...ss in F for an s<a
	foundSequenceInF := false;
	s := Maximum(a-1,1);
	c := b-1;
	if not (b eq n) then
		testSequence := B_n.b;
		while c ge s do
			testSequence := testSequence * B_n.c;
			c := c-1;
		end while;
		while (s lt a) and (s gt 1) and (s lt n) do
			if testSequence in F then
				foundSequenceInF := true;
				break;
			end if;
			s := s-1;
			testSequence := testSequence*B_n.s;
		end while;
	else
		foundSequenceInF := false;
	end if;
	if ((b eq n) or (foundSequenceInF)) then
		a_ := a-1;
		b_ := b;
		alpha := alpha+1;
		
		//Test Sa'b' = Sab#
		Sabeq := Sab_eq_Test(F,a,b,a_,b_,n);
		if Sabeq eq {} then
			for l in [0..k-t] do
				for r in [0..n-1] do
					for s in [0..n-1] do
						if r eq 0 then
							T_[l+1][r+1][s+1] := &+[T[l+1][u+1][s+1]: u in [0..b-a]];
						else
							T_[l+1][r+t][s+t] := 0;
						end if;
					end for;
				end for;
			end for;
		elif Sabeq eq {B_n.(a-1)} then
			for l in [0..k-t] do
				for r in [0..n-1] do
					for s in [0..n-1] do
						if r eq 0 then
							T_[l+1][r+1][s+1] := &*[T[l+1][u+1][s+1]: u in [0..b-a]];
						elif (1 le r) and (r le l) and (r+s lt alpha) then 
							T_[l+1][r+1][s+1] := -T[l-r+1][r][s+1];
						elif (1 le r) and (r le l) and (r+s gt alpha) and (s gt 0) then
							T_[l+1][r+1][s+1] := -T[l-r+1][r][s];
						else
							T_[l+1][r+1][s+1] := 0;
						end if;
					end for;
				end for;
			end for;
		end if;
	else
		a_ := a;
		b_ := b+1;
		alpha := alpha + 1;

		//Test Sa'b' = Sab#
		Sabeq := Sab_eq_Test(F,a,b,a_,b_,n);

		if Sabeq eq {} then
			for l in [0..k-t] do
				for r in [0..n-1] do
					for s in [0..n-1] do
						if s eq 0 then 
							T_[l+1][r+1][s+1] := &+[T[l+1][r+1][u+1]: u in [0..b-a]];
						else
							T_[l+1][r+1][s+1] := 0;
						end if;
					end for;
				end for;
			end for;
		elif Sabeq eq {B_n.b} then
			for l in [0..k-t] do
				for r in [0..n-1] do
					for s in [0..n-1] do
						if s eq 0 then
							T_[l+1][r+1][s+1] := &+[T[l+1][r+1][u+1]: u in [0..b-a]];
						elif (1 le s) and (s le l) and (r+s lt alpha) then
							T_[l+1][r+1][s+1] := -T[l-s+1][r+1][s];
						elif (1 le s) and (s le l) and (r+s gt alpha) and (r gt 0) then
							T_[l+1][r+1][s+1] := -T[l-s+1][r][s];
						else
							T_[l+1][r+1][s+1] := 0;
						end if;
					end for;
				end for;
			end for;
		elif Sabeq ne {B_n.(a-1)} then
			for l in [0..k-t] do
				for r in [0..n-1] do
					for s in [0..n-1] do
						if s eq 0 then
							T_[l+1][r+1][s+1] := &+[T[l+1][r+1][u+1]: u in [0..b-a]];
						elif (s eq alpha) and (l ge alpha) and (r ge 1) then
							T_[l+1][r+1][s+1] := - &+[T[l-alpha+1][r][u+1]: u in [0..b-a]];
						else
							T_[l+1][r+1][s+1] := 0;
						end if;
					end for;
				end for;
			end for;
		end if;
	end if;
	a := a_;
	b := b_;
	T := T_;
end while;
retval := &+[ (&+[T[l+1][r+1][s+1]:r in [0..n-1], s in [0..n-1] ])*xnj[k-t-l+1] : l in [0..k-t]  ];
return retval;
end function;
