# Returns a generator for GL_{n-1}(Z[omega]/2)

R := GF(4);
omega := Z(4);

ReducedBurauGen := function(n, i)
    local M;

    if n < 2 then
        Error("n must be at least 2");
    fi;
	
    if i < 1 or i > n - 1 then
        Error("i must be in the range 1 to n - 1");
    fi;

    M := IdentityMat(n - 1, R);

    if i = 1 then
        # Top-left 2x2 block:
        # [ -t  1 ]
        # [  0  1 ]
        M[1][1] := -omega; # t = omega 
        M[1][2] := omega^0 * 1;
        M[2][2] := omega^0 * 1;

    elif i = n - 1 then
        # Bottom-right 2x2 block:
        # [ 1   0 ]
        # [ t  -t ]
        M[n - 2][n - 2] := omega^0 * 1;
        M[n - 1][n - 2] := omega; # t = omega;
        M[n - 1][n - 1] := -omega; # t = omega 

    else
        # Middle 3x3 block:
        # [ 1   0   0 ]
        # [t   -t   1 ]
        # [ 0   0   1 ]
        M[i - 1][i - 1] := omega^0;
        M[i][i - 1]     := omega;
        M[i][i]         := -omega;
        M[i][i + 1]     := omega^0;
        M[i + 1][i + 1] := omega^0;
    fi;

    return M;
end;

BurauFormedGroup := function(n)
	local G,U,gens,H;
	G := GL(n-1,R);
	gens := List([1..(n-1)], x -> ReducedBurauGen(n,x));
	H := Subgroup(G,gens);
	U := GU(n-1,2);
	return [G,U,H,gens];
end;
