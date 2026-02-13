load("lib.sage")

(Gr_3, gens3, reps3, timeout, comp_verts, _) = construct_orbit_hood(3,N=3000,M=3000,breadth=True, rho = rho_std(3))

(Gr_4, gens4, reps4, timeout, comp_verts, _) = construct_orbit_hood(4,N=3000,M=3000,breadth=True, rho = rho_std(4))


stab3_gens = build_stab_gens(Gr_3,gens3)
identity = gens3[0]*gens3[0]**(-1)
stab3_gens = [s for s in stab3_gens if s != identity]

stab4_gens = build_stab_gens(Gr_4,gens4)
identity = gens4[0]*gens4[0]**(-1)
stab4_gens = [s for s in stab4_gens if s != identity]

print("Generators for Br(q_3)")
print(stab3_gens,"\n")

print("Generators for Br(q_4)")
print(stab4_gens)
