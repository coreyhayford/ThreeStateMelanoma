from pysb import *  # import everything from pysb
from pysb.util import *  # for alias_model_components()

def declare_monomers():  # Monomer("name_of_monomer")
    Monomer("A")
    Monomer("B")
    Monomer("C")


def declare_parameters():  # Parameter("name_of_par", value)
    Parameter("A0")
    Parameter("B0")
    Parameter("C0")

    Parameter("k_div_A")#, .033)  # np.random.normal(.03, .005) )
    Parameter("k_div_B")#, .033)  # np.random.normal(.05,.005) )
    Parameter("k_div_C")#, .033)  # np.random.normal(.07,.005) )

    Parameter('k_death_A')#, .005)
    Parameter('k_death_B')#, .005)
    Parameter('k_death_C')#, .005)

    Parameter("k_AB",     0.025)
    Parameter("k_BA",     0.00004)
    Parameter("k_CB",     0.025)
    Parameter("k_BC",     0.00004)

    alias_model_components()  # "Make all model components visible as symbols in the caller's global namespace"

def declare_initial_conditions():  # Initial(call monomer's name, call parameter's name for IC value = number of initial cells)
    Initial(A, A0)
    Initial(B, B0)
    Initial(C, C0)


def declare_observables():  # Observables("name_of_obs", call monomer name you are observing)
    Observable("Obs_A", A)
    Observable("Obs_B", B)
    Observable("Obs_C", C)
    Observable("Obs_AB", A() + B())
    Observable("Obs_BC", B() + C())
    Observable("Obs_AC", A() + C())
    Observable("Obs_All", A() + B() + C())  ## Don't forget parentheses here


def declare_rules():  # Rule("rule_name", chemical reaction, call name of DIP rate/parameter)
    Rule("DIV_A", A() >> A() + A(), k_div_A)  ##Don't forget parentheses here
    Rule("DIV_B", B() >> B() + B(), k_div_B)
    Rule("DIV_C", C() >> C() + C(), k_div_C)
    Rule('Death_A', A() >> None, k_death_A)
    Rule('Death_B', B() >> None, k_death_B)
    Rule("Death_C", C() >> None, k_death_C)


def declare_transition_rules():
    Rule("A_to_B", A() >> B(), k_AB)
    Rule("B_to_A", B() >> A(), k_BA)
    Rule("B_to_C", B() >> C(), k_BC)
    Rule("C_to_B", C() >> B(), k_CB)