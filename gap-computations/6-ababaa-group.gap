G_ababaa_gens := [
];

G_ababaa := Subgroup(B6, G6_ababaa_gens);
G_ababaa_core := Core(B6, G_ababaa);
Q_ababaa := B6/G_ababaa_core;
idx := Index(B6, G_ababaa);
core_idx := Index(B6, G_ababaa_core);
l := LowIndexSubgroups(Sp(6,3), Order(Sp(6,3))/core_idx);
