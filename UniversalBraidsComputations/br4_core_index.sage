load("lib.sage")

(Gr_4, gens4, reps4, timeout, comp_verts, _) = construct_orbit_hood(4,N=3000,M=3000,breadth=True, rho = rho_std(4))

stab_gens = build_stab_gens(Gr_4,gens4)
identity = gens4[0]*gens4[0]**(-1)
stab_gens = [s for s in stab_gens if s != identity]
(huge_Gr, _, rreps, timeout, comp_verts2) = construct_huge_orbit(4,reps4,gens=stab_gens)

print("Computed index [Br(q_4) : Core(Br(q_4))] to be:", len(huge_Gr.vertices()))
