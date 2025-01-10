import itertools
import math
import functools

class FreeGrp:
    gens_labels = []
    gens = []
    def __init__(self, inp,prefix="a"):
        if isinstance(inp, list):
            self.gens_labels = inp
            self.rank = len(self.gens_labels)
        elif isinstance(inp,int) and inp > 0:
            self.gens_labels = [prefix + str(i) for i in range(1,inp+1)]
            self.rank = inp
        else:
            assert False, "Input malformed"
        self.gens = [FreeGrpElement([i], self) for i in range(1,self.rank+1)]

    def single_artin(self, i):
        map_gens = self.gens.copy()
        if i > 0:
            map_gens[i-1] = self.gens[i]
            map_gens[i] = self.gens[i]**(-1) * self.gens[i-1] * self.gens[i]
        elif i < 0:
            j = abs(i)
            map_gens[j-1] = self.gens[j-1]*self.gens[j]* self.gens[j-1]**(-1)
            map_gens[j] = self.gens[j-1]
        return FreeGrpHom(map_gens)

    def identity_hom(self):
        map_gens = self.gens.copy()
        return FreeGrpHom(map_gens)

    def __call__(self,tietze_inp):
        return FreeGrpElement(tietze_inp, self)

class FreeGrpElement:
    def __init__(self,tietze_inp,par):
        this_run = tietze_inp
        any_removed = True
        while any_removed:
            any_removed = False
            for i in range(0,len(this_run)-1):
                if this_run[i] == -this_run[i+1]:
                    this_run.pop(i)
                    this_run.pop(i)
                    any_removed = True
                    break
        self.tietze = this_run
        self.parent = par

    def Tietze(self):
        return self.tietze

    def __mul__(self,other):
        assert isinstance(other,FreeGrpElement) and self.parent == other.parent, "Group Op Not Defined"
        return FreeGrpElement(self.Tietze() + other.Tietze(), self.parent)

    def inverse(self):
        return FreeGrpElement([-x for x in self.tietze[::-1]], self.parent)

    def __pow__(self,other):
        if other == 0:
            return FreeGrpElement([],self.parent)
        elif other > 0:
            return FreeGrpElement(self.tietze*other,self.parent)
        elif other < 0:
            return (self.inverse())**abs(other)

    def __str__(self):
        def str_single(i):
            res = str(self.parent.gens_labels[abs(i)-1])
            if i < 0:
                res += "^(-1)"
            return res

        return '*'.join([str_single(i) for i in self.tietze])

    def __repr__(self):
        return self.__str__()

    def artin_hom(self, grp):
        assert grp.rank == self.parent.rank + 1
        individual_actions = [grp.single_artin(i) for i in self.tietze]
        res = grp.identity_hom()
        for x in individual_actions:
            res = res.then(x)
        return res

    def artin_act(self,other):
        assert isinstance(other,FreeGrpElement) and other.parent.rank == self.parent.rank - 1
        other_hom = other.artin_hom(self.parent)
        return other_hom(self)

class FreeGrpHom:
    map_gens = []
    def __init__(self, send_gens,contra=False):
        for x in send_gens:
            assert hasattr(x,'__mul__')
        self.contravariant=contra
        self.map_gens = send_gens

    def __call__(self, other,power=1):
        assert isinstance(other, FreeGrpElement), "Call must be used on a group element"
        tietze_lst = other.Tietze()
        result = self.map_gens[0]**0

        for i in tietze_lst:
            to_multiply = self.map_gens[abs(i)-1]**int(math.copysign(1,i))
            if self.contravariant:
                result = to_multiply * result
            else:
                result = result * to_multiply
        return result

    def then(self,other):
        assert isinstance(other,FreeGrpHom)
        new_lst = [other(x) for x in self.map_gens]
        return FreeGrpHom(new_lst)

class SurfaceGroup(FreeGrp):
    genus_gens = []
    puncture_gens = []
    g = 0
    n = 0

    def __init__(self, genus, num_punctures):
        self.g = genus
        self.n = num_punctures

        genus_gen_names = [x for i in range(0,self.g) for x in ['a_' + str(i), 'b_' + str(i)]]
        if self.n > 0:
            puncture_gen_names = ['p_' + str(i) for i in range(1,self.n)]

            super().__init__(genus_gen_names + puncture_gen_names)
        else:
            free_grp = FreeGroup(genus_gen_names)
            commutators = [[2*i-1,2*i,1-2*i,-2*i] for i in range(1,self.g+1)]

            # Note: itertools.chain flattens the list
            word = free_grp(itertools.chain(*commutators))
            super().__init__(free_grp, tuple([word]))

        gens = super().gens
        self.genus_gens = gens[0:2*self.g-1]
        if self.n > 0:
            self.puncture_gens = gens[2*self.g:]

