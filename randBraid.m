load "xnk.m";

function randBraid(n, k)
/*---------------------------------------------------------------*/
/*Generates random elements in B_n+ of length k according to 	 */
/*'Generating random Braids' Paper				 */
/*A better implementation can be found as Random(B_N, k) in Magma*/
/*---------------------------------------------------------------*/

//Compute xn,k
//Compute xmj, hmj for 0<=m<=n, 0<=j<=k
////hmj:
h := [];
initvec := [1];
for i in [1..k] do
	Append(~initvec, 0);
end for;
for i in [0..n] do
	Append(~h, initvec);
end for;
//Calculate from third line on, first two lines are correctly filled
for m in [2..n] do
	for j in [1..k] do
		hmj := 0;
		//print "h",m,",",j ;
		for i in [1..m] do
			//print "i = ",i;
			x := m-i;
			y := j - Binomial(i, 2);
			if y lt 0 then
			//	print j," - ",Binomial(i,2)," = ",y;
				break;
			end if;
			hmj := hmj + (-1)^(i+1)*h[x+1][y+1];
			//print "+(-1)^",i+1,"*h",x,",",y," (=",h[x+1][y+1],")";
		end for;
		h[m+1][j+1] := hmj;
	end for;
end for;

////xmj
x := [1]; //xn0 = 1
for m in [n..n] do
	for j in [1..k] do
		xmj := 0;
		//print "x",m,",",j;
		for i in [1..j] do 
			y := j-i;
			if i le Binomial(m, 2) then
				xmj := xmj + x[y+1]*h[m+1][i+1];
			//	print "+ x",y,"(=",x[y+1],")*h",m,",",i," (=",h[m+1][i+1],")";
			end if;
		end for;
		Append(~x, xmj * (-1));	
	end for;
end for;

//print h;
//print x;

SetSeed(Truncate(Realtime()));
r := Random(1, x[k+1]);

//print "chose r = ",r;

B_n := BraidGroup(n);
w := Id(B_n);
v := x[k+1] - r;
a := 0;

for l in [1..k] do
	a := Max(a-1, 1);
	b := n-1;
	mu := 0;
	while a lt b do
		m := Floor((a+b)/2);
		xnk_wm := xnk(n,k,w,m,x);
		if xnk_wm le v then
			b:=m;
			mu := xnk_wm;
		else
			a := m+1;
		end if;
	end while;
	if (a gt 0) and (a lt n) then
		w := w*B_n.a;
	else
		print "Unexpected error: Can not get ",a,"-th artin generator";
	end if;
	v := v-mu;
end for;

return w;

end function;
