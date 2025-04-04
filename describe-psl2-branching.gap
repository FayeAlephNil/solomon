construct_mats := function (p,n)
	local A,B,C,e;
	e := One(GF(p));
	A := e*[[1,-1],[0,1]];
	B := e*[[1,0],[1,1]];
	if IsEvenInt(n) then
		C := (A*B)^(n/2);
	else
		C := (A*B)^((n-1)/2) * A;
	fi;
	return [A,B,C];
end;

describe_branching_perm := function (p,n)
	local ABC, A,B,C, e,G, H, PG,
		hom1,
		Lst_PG, A_lst, B_lst,C_lst;
	e := One(GF(p));
	ABC := construct_mats(p,n);
	A := ABC[1];
	B := ABC[2];
	C := ABC[3];

	G := SL(2,p);
	H := Group(e*[[-1,0],[0,-1]]);
	hom1 := NaturalHomomorphismByNormalSubgroup(G,H);
	PG := Image(hom1);
	Lst_PG := List(PG);
	A_lst := List(Lst_PG,x -> Position(Lst_PG, hom1(A)*x));
	B_lst := List(Lst_PG,x -> Position(Lst_PG, hom1(B)*x));
	C_lst := List(Lst_PG,x -> Position(Lst_PG, hom1(C)*x));
	return [PermList(A_lst),PermList(B_lst),PermList(C_lst)];
end;

describe_branching_cycle := function (p,n)
	local perm_lst;
	perm_lst := describe_branching_perm(p,n);
	return List(perm_lst, x -> CycleStructurePerm(x));
end;

compute_cover_genus_bdry := function (p,n)
	local psl2fp_order,cycle_typ_tmp,cycle_typ,num_cyc_ABC,cone,euler_char,b,i,x,j;
	if p < 3 then
		psl2fp_order := 6;
	else
		psl2fp_order := p*(p^2-1)/2;
	fi;
	cycle_typ_tmp := describe_branching_cycle(p,n);
	cycle_typ := [];
	for j in [1..Length(cycle_typ_tmp)] do
		cycle_typ[j] := [];
		for i in [1..Length(cycle_typ_tmp[j])] do
			if (not IsBound(cycle_typ_tmp[j][i])) then
				cycle_typ[j][i] := 0;
			else
				cycle_typ[j][i] := cycle_typ_tmp[j][i];
			fi;
		od;
	od;
	num_cyc_ABC := List(cycle_typ,x -> Sum(x));
	b := num_cyc_ABC[3];
	cone := num_cyc_ABC[1]/psl2fp_order;
	euler_char := psl2fp_order*(1-n*(1-cone));
	return rec(genus := (euler_char + b - 2)/(-2), bound := b, fix := n*num_cyc_ABC[1]);
end;
