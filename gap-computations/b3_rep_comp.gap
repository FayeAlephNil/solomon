LoadPackage("ACE");
TCENUM := ACETCENUM;;
LoadPackage("FR");

B3 := SurfaceBraidFpGroup(3,0,1);

anti := B3.1^3;
equi := B3.1^(-1) * B3.2 * B3.1;
conj := anti*equi;
anti_conj := conj*anti*conj^(-1);

S := anti*anti_conj * anti;
T := anti_conj^(-1);
G := Subgroup(B3,[anti,equi]);

M_anti := [ [1, 0], [1, 1] ];
M_equi := [ [-1, 1], [0, -1] ];

hom := GroupHomomorphismByImages(G, SL(2,Integers), [M_anti,M_equi]);
# H := NormalClosure(G,[S^2 * (S*T)^(-3), S^4, (S*T)^(6)]);
