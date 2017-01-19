
Original Sequence: {{dna_input}}
Peptide: {{ori_pep}}

Restriction map for original sequence


{{original_map}}

=========================




{% for item in bakpeps_out %}

Original sequence enzymes: {{oname}}

Proposed sequence enzymes: {{item.pames}}
Enzimes only in original sequence: {{item.o}}
Enzimes only in proposed sequence: {{item.p}}
=========================



{% endfor %}
