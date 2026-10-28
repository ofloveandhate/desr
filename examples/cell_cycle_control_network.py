"""
The "Cell cycle control network" example from the paper.

Adapted from: J. C. S. . J. J. Tyson, Mathematical modeling as a tool for investigating cell cycle control networks,
2007.


"""

from desr import ODESystem, ODETranslator


# the tex from our scaling symmetries paper, faithfully transcribed from: J. C. S. . J. J. Tyson, Mathematical modeling as a tool for investigating cell cycle control networks,
system_tex = \
      r'''\frac{d}{dt}\textrm{[Cyclin]} &= k_1 - k_2 \textrm{[Cyclin]} - k_3 \textrm{[Cyclin]} \textrm{[Cdk]} \\
\frac{d}{dt}\textrm{[MPF]} &= k_3 \textrm{[Cyclin]} \textrm{[Cdk]} - k_2 \textrm{[MPF]} - k_{\textrm{wee}} \textrm{[MPF]} + k_{25} \textrm{[preMPF]} \\
\frac{d}{dt}\textrm{[preMPF]} &= -k_2 \textrm{[preMPF]} + k_{\textrm{wee}} \textrm{[MPF]} - k_{25} \textrm{[preMPF]} \\
\frac{d}{dt}\textrm{[Cdc25P]} &= \frac{k_a \textrm{[MPF]}\left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right)}{K_a + \textrm{[total Cdc25]} - \textrm{[Cdc25P]}} - \frac{k_b \textrm{[PPase]} \textrm{[Cdc25P]}}{K_b + \textrm{[Cdc25P]}} \\
\frac{d}{dt}\textrm{[Wee1P]} &= \frac{k_e \textrm{[MPF]}\left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right)}{K_e + \textrm{[total Wee1]} - \textrm{[Wee1P]}} - \frac{k_f \textrm{[PPase]} \textrm{[Wee1P]}}{K_f + \textrm{[Wee1P]}} \\
\frac{d}{dt}\textrm{[IEP]} &= \frac{k_g \textrm{[MPF]}\left( \textrm{[total IE]} - \textrm{[IEP]} \right)}{K_g + \textrm{[total IE]} - \textrm{[IEP]}} - \frac{k_h \textrm{[PPase]} \textrm{[IEP]}}{K_h + \textrm{[IEP]}} \\
\frac{d}{dt}\textrm{[APC]} &= \frac{k_c \textrm{[MPF]}\left( \textrm{[total APC]} - \textrm{[APC]} \right)}{K_c + \textrm{[total APC]} - \textrm{[APC]}} - \frac{k_d \textrm{[PPase]} \textrm{[APC]}}{K_d + \textrm{[APC]}} \\
\textrm{[Cdk]} &= \textrm{[total Cdk]} - \textrm{[MPF]} - \textrm{[preMPF]} \\
k_{25} &= V_{25}' \left( \textrm{[total Cdc25]} - \textrm{[Cdc25P]} \right) + V_{25}'' \textrm{[Cdc25P]} \\
k_{\textrm{wee}} &= V_{\textrm{wee}}' \textrm{[Wee1P]} + V_{\textrm{wee}}'' \left( \textrm{[total Wee1]} - \textrm{[Wee1P]} \right) \\
k_2 &= V_2' \left( \textrm{[total APC]} - \textrm{[APC]} \right) + V_2'' \textrm{[APC]}'''


# The tex above is the paper's, unaltered.  `from_tex` sorts out which lines are equations
# of the system and which are definitions -- the last four give Cdk, k_25, k_wee and k_2 in
# terms of species that move -- and `sub_defns` asks for those to be substituted away.
# The primes in V_{25}' and V_{25}'' are parts of names here, not derivatives, which is
# `from_tex`'s default reading.
system = ODESystem.from_tex(system_tex, sub_defns=True)

print(system)
# silviana has faith in this code, up to this point.  she manually checked the ODEs in `system`, and they match what's in the 2007 paper.
# (that check was made on the system before the definitions were substituted in; to see
#  that form again, pass only the seven derivative lines to from_tex.)


# the commented-out ones have been substituted away
variable_order = \
(
't', # time
# the seven variables
'Cyclin',
'MPF',
'preMPF',
'Cdc25P',
'Wee1P',
'IEP',
'APC', 
# start the constants
'PPase',
'K_a',
'K_b',
'K_c',
'K_d',
'K_e',
'K_f',
'K_g',
'K_h',
'k_1',
# 'Cdk',
# 'k_2',
# 'k_25',
# 'k_wee',
'k_3',
'k_a',
'k_b',
'k_c',
'k_d',
'k_e',
'k_f',
'k_g',
'k_h',
'V_25prime', 
'V_25primeprime',
'V_weeprime', 
'V_weeprimeprime', 
# 'APCstar', # this was removed, because APC in the paper was recognized as a typo. APC* is supposed to be APC, as evidenced in the code in Table 1.
'V_2prime', 
'V_2primeprime',
'totalCdc25',
'totalIE',
'totalAPC',
'totalCdk', 
'totalWee1'
)

system.reorder_variables(variable_order)

print('\n\n---------\n\nsystem in the chosen variable order:')
print(system)

print(f'\n{system.variables=}')

print(f'\n{system.parameters=}')
print(f'{len(system.parameters)} parameters')

print(f'\n{system.dependent_variables=}')
print(f'{len(system.dependent_variables)} dependent variables')

print('\nVariable order: ', system.variables)

translator = ODETranslator.from_ode_system(system, new_indices_start_at=1)

print('\nScaling matrix',translator.scaling_matrix.__repr__())
print('Scaling matrix size',translator.scaling_matrix.shape)

print('\nHermite multiplier', translator.herm_mult.__repr__())

print('\nInvariants: ', translator.invariants())
print(f'{translator.invariants().shape} invariants')

print('\nsubstitutions:',translator.translate_parameter_substitutions(system))

system_translated = translator.translate(system)

# i really want to factor the system, cuz it's too long for a line in the paper.
system_translated.simplify_derivatives()

print('\nReduced system:')
print(system_translated.to_tex())#.to_tex()





print(f'{system_translated.variables=}')
print(f'{system_translated.parameters=}')
print(f'{len(system_translated.parameters)}  parameters')

print(f'{system_translated.dependent_variables=}')
print(f'{len(system_translated.dependent_variables)}  dependent variables')




# notes for work to do

# make when init cond is 0 and have divide by 0, raise and say to use a different variable order
# forbid to make substitutions with variables -- this will cause a bunch of unit tests to fail..
# assume all new symbols are constants