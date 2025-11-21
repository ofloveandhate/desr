"""
The "Cell cycle control network" example from the paper.

Adapted from: J. C. S. . J. J. Tyson, Mathematical modeling as a tool for investigating cell cycle control networks,
2007.


"""

import sympy
from desr.matrix_normal_forms import smf
from desr.ode_system import ODESystem
from desr.ode_system import maximal_scaling_matrix, rational_expr_to_power_matrix, hnf_col, hnf_row, normal_hnf_col
from desr.matrix_normal_forms import normal_hnf_row
from desr.ode_translation import ODETranslation, scale_action
from desr.tex_tools import expr_to_tex, matrix_to_tex
import desr.tex_tools as tex_tools



def pre_process_latex(system_tex):
	import re

	# get rid of textrm
	pattern = r"\\textrm\{([^}]*)\}"
	text = re.sub(pattern,r'\1',system_tex)


	# replaces spaces inside square braces.  
	# this is scary to me, because sometimes [] is used in place of parentheses in tex, but not in these eqns so it's ok
	text = re.sub(r"\[([^\]]+)\]", lambda m: "[" + m.group(1).replace(" ", "") + "]", text)


	# this codeblock removed because the APC* variable was recognized as a typo.
	# # replace '*' with 'star' only inside [...] 
	# text = re.sub(r"\[([^\]]*?)\]", lambda m: "[" + m.group(1).replace("*", "star") + "]", text)



	# replace [] with ()
	text = text.replace('[','(').replace(']',')')

	# strip left/right
	text = re.sub(r"\\left\s*\(", "(", text)
	text = re.sub(r"\\right\s*\)", ")", text)


	# combine variables with spaces between []
	pattern = re.compile(r"\(\s*([^()+=\-*/]+?)\s*\)")

	while True:
	    new_text = pattern.sub(r"\1", text)
	    if new_text == text:
	        break
	    text = new_text


	# turn d/dt x into dx/dt
	# pattern: \frac{d}{dt}Variable
	pattern = r"\\frac\{d\}\{dt\}([a-zA-Z0-9]+)"

	# replace with \frac{dVariable}{dt}
	text = re.sub(pattern, r"\\frac{d\1}{dt}", text)
	# print('after moving d/dt x --> dx/dt:\n',text,'\n\n')


	# repnace Variable2( with Variable2*(
	pattern = r"([a-zA-Z0-9]+)\("
	text = re.sub(pattern, r"\1*(", text)
	# print('after stars before ():\n',text,'\n')

	# Make multi-char subscipts just underscores.  Flatten all other {...} after underscore
	text = re.sub(r"_\{([^}]*)\}", r"_\1", text)

	# replace ' with prime
	text = re.sub(r"'", r"prime", text)

	# finally done doing regex stuff to preprocess so it makes the system from tex ok.
	return text



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


text = pre_process_latex(system_tex)


num_dynamic_eqns = 7
num_constants = 4

lines = text.split(r'\\')
assert(len(lines) == num_dynamic_eqns + num_constants)

system_tex_dynamic = r'\\'.join(lines[:num_dynamic_eqns])


system = ODESystem.from_tex(system_tex_dynamic)

print(system)
# silviana has faith in this code, up to this point.  she manually checked the ODEs in `system`, and they match what's in the 2007 paper.


# now we need to make the substitions.
# first we construct a dictionary of them
subme = {}
for line in lines[num_dynamic_eqns:]:
	lhs,rhs = line.split('&=')
	rhs = tex_tools._tex_to_sympy_one_line(rhs)
	print(f'{lhs} = {rhs}')
	subme[lhs] = rhs

# then we actually make the substitutions (using Sympy behind the scene)
system = system.diff_subs(subme)


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

print('\n\n---------\n\nsystem after substitutions:')
print(system)

print(f'\n{system.variables=}')

print(f'\n{system.constant_variables=}')
print(f'{len(system.constant_variables)} constant variables')

print(f'\n{system.non_constant_variables=}')
print(f'{len(system.non_constant_variables)} non-constant variables')

print('\nVariable order: ', system.variables)

translation = ODETranslation.from_ode_system(system, new_indices_start_at=1)

print('\nScaling matrix',translation.scaling_matrix.__repr__())
print('Scaling matrix size',translation.scaling_matrix.shape)

print('\nHermite multiplier', translation.herm_mult.__repr__())

print('\nInvariants: ', translation.invariants())
print(f'{translation.invariants().shape} invariants')

print('\nsubstitutions:',translation.translate_parameter_substitutions(system))

system_translated = translation.translate(system)

# i really want to factor the system, cuz it's too long for a line in the paper.
system_translated.simplify_derivatives()

print('\nReduced system:')
print(system_translated.to_tex())#.to_tex()





print(f'{system_translated.variables=}')
print(f'{system_translated.constant_variables=}')
print(f'{len(system_translated.constant_variables)}  constant variables')

print(f'{system_translated.non_constant_variables=}')
print(f'{len(system_translated.non_constant_variables)}  non-constant variables')




# notes for work to do

# make when init cond is 0 and have divide by 0, raise and say to use a different variable order
# forbid to make substitutions with variables -- this will cause a bunch of unit tests to fail..
# assume all new symbols are constants