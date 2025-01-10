#from .surface_group import SurfaceGroup
from .surface_group import FreeGrp, FreeGrpElement, FreeGrpHom

from functools import reduce
from operator import mul
import math

class MonodromyRep:
    domain = None
    underlying_lst = []
    contravariant = False

    @classmethod
    def action(cls, mclass, rep):
        assert isinstance(rep, MonodromyRep), "Action is Only for Monodromy Reps"
        return rep.twist(mclass)

    def __init__(self, source, map_gens, contra=False):
        #assert isinstance(source, SurfaceGroup), "Source must be a surface group"
        self.domain = source
        self.underlying_lst = map_gens
        self.contravariant = contra

    def twist(self, mclass):
        assert self.domain.g == 0, "Positive Genus Base not yet Implemented"
        assert isinstance(mclass,FreeGrpElement), "Currently for Free base Groups, braids must be passed in directly"

        # This workaround is because technically self.domain is a FP group, not a free group
        mclass_hom = mclass.artin_hom(self.domain)
        workaround_lst = [self.domain(mclass_hom(x).Tietze()) for x  in self.domain.gens]

        new_map = [self(a) for a in workaround_lst]
        return MonodromyRep(self.domain, new_map, self.contravariant)

    def homology_matrices(self):
        return [g.homology_matrix() for g in self.underlying_lst]

    def to_local_system(self):
        return MonodromyRep(self.domain, self.homology_matrices(), self.contravariant)

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
        return hash(tuple(self.underlying_lst))

    def __eq__(self, other):
        doms = self.domain == other.domain
        contras = self.contravariant == other.contravariant
        lsts = self.underlying_lst == other.underlying_lst
        return doms and contras and lsts
