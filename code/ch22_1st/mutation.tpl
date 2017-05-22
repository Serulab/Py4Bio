Original Sequence: {{dna_input}}
Peptide: {{ori_pep}}

Restriction map for original sequence:

{{ori_map}}

=========================
{% for item in bakpeps_out %}

Original sequence enzymes: {{oname}}

{{item.graph}}

Proposed sequence enzymes: {{item.pames}}
Enzimes only in original sequence: {{item.o}}
Enzimes only in proposed sequence: {{item.p}}
=========================



{% endfor %}
