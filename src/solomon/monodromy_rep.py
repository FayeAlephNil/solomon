#from .surface_group import SurfaceGroup
from .surface_group import FreeGrpElement

from functools import reduce
from operator import mul
import math
import numpy as np
import flint

class MonodromyRep:
    domain = None
    underlying_lst = []
    contravariant = False
    homology_rep = False
    homology_mod_rep = -1
    integral = True

    @classmethod
    def action(cls, mclass, rep):
        assert isinstance(rep, MonodromyRep), "Action is Only for Monodromy Reps"
        return rep.twist(mclass)

    def __init__(self, source, map_gens, contra=False,hom_rep=False,hom_mod_rep=-1,integ=True):
        #assert isinstance(source, SurfaceGroup), "Source must be a surface group"
        self.domain = source
        self.underlying_lst = map_gens
        self.contravariant = contra
        self.homology_rep = hom_rep
        self.homology_mod_rep = hom_mod_rep
        self.integral = integ
        #if self.homology_mod_rep != -1:
        #    self.underlying_lst[:] = [mat % self.homology_mod_rep for mat in self.underlying_lst]

    def twist(self, mclass,debug=False):
        assert self.domain.g == 0, "Positive Genus Base not yet Implemented"
        assert isinstance(mclass,FreeGrpElement), "Currently for Free base Groups, braids must be passed in directly"

        # This workaround is because technically self.domain is a FP group, not a free group
        mclass_hom = mclass.artin_hom(self.domain)
        workaround_lst = [self.domain(mclass_hom(x).Tietze()) for x  in self.domain.gens]

        new_map = [self(a) for a in workaround_lst]
        if self.homology_rep and self.integral and not(self.homology_mod_rep):
            new_map[:] = [np.asmatrix(mat, dtype=int) for mat in new_map]
        return MonodromyRep(self.domain, new_map, self.contravariant,hom_rep=self.homology_rep,hom_mod_rep=self.homology_mod_rep,integ=self.integral)

    def homology_matrices(self):
        if self.homology_rep:
            return self.underlying_lst
        if self.integral:
            return [np.asmatrix(g.homology_matrix(),dtype=int) for g in self.underlying_lst]
        else:
            return [np.asmatrix(g.homology_matrix(),dtype=float) for g in self.underlying_lst]

    def homology_matrices_mod(self,N):
        if self.homology_mod_rep != -1:
            return self.underlying_list
        mats = self.homology_matrices()
        return [flint.nmod_mat(mat.shape[0], mat.shape[1], mat.flatten().tolist()[0],N) for mat in mats]

    def to_local_system(self):
        if self.homology_rep:
            return self
        return MonodromyRep(self.domain, self.homology_matrices(), self.contravariant,hom_rep=True)

    def to_local_system_mod(self,N):
        if self.homology_mod_rep != -1:
            return self
        return MonodromyRep(self.domain, self.homology_matrices_mod(N), self.contravariant,hom_rep=True,hom_mod_rep=N,integ=self.integral)

    def total_monodromy(self):
        return reduce(mul, self.underlying_lst)

    def __call__(self,other,power=1):
        assert isinstance(other, FreeGrpElement), "Call must be used on another group element"
        tietze_lst = other.Tietze()
        result = self.underlying_lst[0]**0

        for i in tietze_lst:
            to_multiply = self.underlying_lst[abs(i)-1]**int(math.copysign(1,i))
            if self.contravariant:
                result = to_multiply * result
            else:
                result = result * to_multiply
        return result

    def __str__(self):
        return str(self.underlying_lst)

    def __repr__(self):
        return repr(self.underlying_lst)

    def __hash__(self):
        if self.homology_mod_rep != -1:
            return hash(tuple([int(i) for mat in self.underlying_lst for i in mat.entries()]))
        elif self.homology_rep:
            return hash(tuple([tuple(mat.flatten().tolist()[0]) for mat in self.underlying_lst]))
        else:
            return hash(tuple(self.underlying_lst))

    def __eq__(self, other):
        doms = self.domain == other.domain
        contras = self.contravariant == other.contravariant
        if len(self.underlying_lst) != len(other.underlying_lst):
            return False
        elif self.homology_rep:
            lsts = True
            for i in range(len(self.underlying_lst)):
                this_one = self.underlying_lst[i]
                other_one = other.underlying_lst[i]
                if self.homology_mod_rep:
                    compare_them = this_one == other_one
                else:
                    compare_them = np.all(np.isclose(self.underlying_lst[i], other.underlying_lst[i],atol = 0.1))
                lsts = lsts and compare_them
        else:
            lsts = self.underlying_lst == other.underlying_lst
        return doms and contras and lsts
