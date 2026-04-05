LoadPackage("ACE");
TCENUM := ACETCENUM;;
LoadPackage("FR");

B4 := SurfaceBraidFpGroup(4,0,1);
x := B4.1 * B4.2 * B4.3;
y := B4.2 * B4.1;
Print(Index(B4,Subgroup(B4,[x,y])));

# F := FreeGroup(6);
#F := FreeGroup(4);
#homC2 := GroupHomomorphismByImages(F, List(GeneratorsOfGroup(F), x -> Z(3)));
#K := Kernel(homC2);
## nice_gens := [F.1^2, F.2^2, F.3^2, F.4^2, F.5^2, F.6^2, F.1*F.2, F.2*F.3, F.3*F.4,F.4*F.5, F.5*F.6];
## hom_gens := [F.1*F.2, F.2*F.3, F.3*F.4,F.4*F.5, F.5*F.6];
#nice_gens := [F.1^2, F.2^2, F.3^2, F.4^2, F.1*F.2, F.2*F.3, F.3*F.4];
#hom_gens := [F.1*F.2, F.2*F.3, F.3*F.4];
#
#KK := Subgroup(F,nice_gens);
## homS3 := GroupHomomorphismByImages(F, [(1,2), (1,3), (1,2), (1,3), (1,2), (1,3)]);
#homS3 := GroupHomomorphismByImages(F, [(1,2), (1,3), (1,2), (1,3)]);
#aut1 := GroupHomomorphismByImages(F, [F.2, F.2^(-1) * F.1 * F.2, F.3, F.4]);
#aut2 := GroupHomomorphismByImages(F, [F.1, F.3, F.3^(-1) * F.2 * F.3, F.4]);
#aut3 := GroupHomomorphismByImages(F, [F.1, F.2, F.4, F.4^(-1) * F.3 * F.4]);
## aut4 := GroupHomomorphismByImages(F, [F.1, F.2, F.3, F.5 ,F.5^(-1) * F.4 * F.5]);
## aut5 := GroupHomomorphismByImages(F, [F.1, F.2, F.3, F.4 ,F.6, F.6^(-1) * F.5 * F.6]);
#our_aut := Subgroup(AutomorphismGroup(F), [aut1,aut2,aut3]);
#
#OnLst := function (f,g)
#	local hom;
#	hom := GroupHomomorphismByImages(F, f);
#	return List([F.1, F.2, F.3, F.4], x -> hom(g(x)));
#end;
#
## OnLst([(1,2),(1,3),(1,2),(1,3),(1,2)], aut1);
#lst_orb := Orbit(our_aut, [(1,2),(1,3),(1,2),(1,3)], OnLst);
#hom_orb := Set(List(lst_orb, lst -> List(hom_gens, x -> (GroupHomomorphismByImages(F,lst))(x) )));
