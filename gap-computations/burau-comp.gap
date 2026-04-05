BurauGenerator := function(n, i)
    local M, B;
    M := IdentityMat(n);
    B := [ [2, -1], [1, 0] ];  # Burau block at t = -1

    M[i][i] := B[1][1];       # top-left of block
    M[i][i+1] := B[1][2];     # top-right
    M[i+1][i] := B[2][1];     # bottom-left
    M[i+1][i+1] := B[2][2];   # bottom-right

    return M;
end;

ReducedBurauGenerator := function(n, i)
    local M;

    if n < 2 then
        Error("n must be at least 2");
    fi;
	
    if i < 1 or i > n - 1 then
        Error("i must be in the range 1 to n - 1");
    fi;

    M := IdentityMat(n - 1, Integers);

    if i = 1 then
        # Top-left 2x2 block:
        # [ -t  1 ]
        # [  0  1 ]
        M[1][1] := 1; # t = -1 
        M[1][2] := 1;
        M[2][2] := 1;

    elif i = n - 1 then
        # Bottom-right 2x2 block:
        # [ 1   0 ]
        # [ t  -t ]
        M[n - 2][n - 2] := 1;
        M[n - 1][n - 2] := -1; # t = -1;
        M[n - 1][n - 1] := 1; # t = -1

    else
        # Middle 3x3 block:
        # [ 1   0   0 ]
        # [-1   1   1 ]
        # [ 0   0   1 ]
        M[i - 1][i - 1] := 1;
        M[i][i - 1]     := -1;
        M[i][i]         := 1;
        M[i][i + 1]     := 1;
        M[i + 1][i + 1] := 1;
    fi;

    return M;
end;

ReducedBurauOrbs := function (n)
	local G, m, F;
	m := n-1;
	F := GF(3);
	G := Group(List([1..(n-1)], x -> Z(3)^0 * ReducedBurauGenerator(n, x)));
	return Orbits(G, F^m, OnPoints);
end;

BurauOrbs := function (n)
	local G, m, F;
	F := GF(3);
	G := Group(List([1..(n-1)], x -> Z(3)^0 * BurauGenerator(n, x)));
	return Orbits(G, F^n, OnPoints);
end;
