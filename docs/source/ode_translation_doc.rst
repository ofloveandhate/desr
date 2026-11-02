ODETranslator and Reduction of ODESystems
==========================================

The :class:`~desr.ode_translation.ODETranslator` is the class used to represent the scaling symmetries of a system.

Creating ODE Translations
-------------------------

.. autoclass:: desr.ode_translation.ODETranslator
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.from_ode_system
    :noindex:

Reduction Methods
-----------------

.. automethod:: desr.ode_translation.ODETranslator.translate
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.translate_parameter
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.translate_dep_var
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.translate_general
    :noindex:

Reverse Translation
-------------------

Reverse translation is the process of taking solutions of the reduced system and recovering solutions of the original system.

.. automethod:: desr.ode_translation.ODETranslator.reverse_translate
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.reverse_translate_parameter
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.reverse_translate_dep_var
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.reverse_translate_general
    :noindex:


Extending from a set of Invariants
----------------------------------

.. automethod:: desr.ode_translation.ODETranslator.extend_from_invariants
    :noindex:


Useful Attributes
-----------------

.. automethod:: desr.ode_translation.ODETranslator.invariant_monomials
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.invariant_names
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.auxiliary_monomials
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.auxiliary_names
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.rewrite_rules
    :noindex:


Finding a column to operate on
------------------------------

.. automethod:: desr.ode_translation.ODETranslator.column_of
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.column_monomials
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.column_names
    :noindex:



Output Functions
----------------

.. automethod:: desr.ode_translation.ODETranslator.to_tex
    :noindex:

Advanced Methods
----------------

These methods will be familiar to those who use Lie groups to analyse more general symmetries of differential equations.
For more information, see :cite:`Fels1999` or :cite:`Hubert2007`.

.. automethod:: desr.ode_translation.ODETranslator.moving_frame
    :noindex:

.. automethod:: desr.ode_translation.ODETranslator.rational_section
    :noindex:
