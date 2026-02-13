load("lib.sage")

(Gr_3, gens3, reps3, timeout, comp_verts, _) = construct_orbit_hood(3,N=3000,M=3000,breadth=True, rho = rho_std(3))

Gr_3.remove_loops()

orbit_3_img = display_orbit(Gr_3,gens3,edge_labels=True)
orbit_3_img.save('braid3_orbit.png')

(Gr_4, gens4, reps4, timeout, comp_verts, _) = construct_orbit_hood(4,N=3000,M=3000,breadth=True, rho = rho_std(4))

Gr_4.remove_loops()
orbit_4_img = display_orbit(Gr_4,gens4,edge_labels=True)
orbit_4_img.save('braid4_orbit.png',figsize=[2*6.4,2*4.8])

Gr_4_center = Gr_4.subgraph([0,1,2,5,6,7,8,11,12])
Gr_4_center.remove_loops()
orbit_4_center_img = display_orbit(Gr_4_center,gens4, edge_labels=True,layout='spring',iterations=10000,heights={
    0: [1,20,0],
    1: [7,8,6],
    2: [3,21,4]
})
orbit_4_center_img.save('braid4_orbit_center.png')


stab_gens = build_stab_gens(Gr_4,gens4)
identity = gens4[0]*gens4[0]**(-1)
stab_gens = [s for s in stab_gens if s != identity]
(huge_Gr, _, rreps, timeout, comp_verts2) = construct_huge_orbit(4,reps4,gens=stab_gens)

print("Computed index [Br(q_4) : Core(Br(q_4))] to be: ", len(huge_Gr.vertices()))
