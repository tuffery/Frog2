(S'47fb6ed7b7cf741395b1915db1acc2bf'
p1
(ihappydoclib.parseinfo.moduleinfo
ModuleInfo
p2
(dp3
S'_namespaces'
p4
((dp5
(dp6
S'name_ring_bond'
p7
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p8
(dp9
g4
((dp10
(dp11
tp12
sS'_exception_info'
p13
(dp14
sS'_parameter_names'
p15
(S'bond'
p16
tp17
sS'_parameter_info'
p18
(dp19
g16
(NNNtp20
ssS'_filename'
p21
S'../python/frowns/Fingerprint/SmartsLinearPaths.py'
p22
sS'_docstring'
p23
S''
sS'_name'
p24
g7
sS'_parent'
p25
g2
sS'_comment_info'
p26
(dp27
sS'_configuration_values'
p28
(dp29
sS'_class_info'
p30
g10
sS'_function_info'
p31
g11
sS'_comments'
p32
S''
sbsS'generatePaths'
p33
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p34
(dp35
g4
((dp36
(dp37
tp38
sg13
(dp39
sg15
(S'molecule'
p40
S'maxdepth'
p41
S'name_atom'
p42
S'name_bond'
p43
S'name_ring_atom'
p44
S'name_ring_bond'
p45
S'make_rpath'
p46
tp47
sg18
(dp48
g45
(I1
S'name_ring_bond'
Ntp49
sg46
(I1
S'0'
Ntp50
sg40
(NNNtp51
sg42
(I1
S'name_atom'
Ntp52
sg41
(I1
S'5'
Ntp53
sg43
(I1
S'name_bond'
Ntp54
sg44
(I1
S'name_ring_atom'
Ntp55
ssg21
g22
sg23
S'(molecule, maxdepth, *name_atom, *name_bond) -> linear paths\n    Generate all linear paths through a molecule up to maxdepth\n    change name_atom and name_bond to name the atoms and bonds\n    in the molecule\n\n    name_atom and name_bond must return a stringable value'
p56
sg24
g33
sg25
g2
sg26
g27
sg28
(dp57
sg30
g36
sg31
g37
sg32
S''
sbsS'countRings'
p58
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p59
(dp60
g4
((dp61
(dp62
tp63
sg13
(dp64
sg15
(S'object'
p65
tp66
sg18
(dp67
g65
(NNNtp68
ssg21
g22
sg23
S'count the number of rings associated with an object,\n    return the number of rings and the number of aromatic rings'
p69
sg24
g58
sg25
g2
sg26
g27
sg28
(dp70
sg30
g61
sg31
g62
sg32
S''
sbsS'_bfswalk'
p71
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p72
(dp73
g4
((dp74
(dp75
tp76
sg13
(dp77
sg15
(S'atom'
p78
S'visitedAtoms'
p79
S'path'
p80
S'rpath'
p81
S'paths'
p82
S'depth'
p83
S'maxdepth'
p84
S'name_atom'
p85
S'name_bond'
p86
S'name_ring_atom'
p87
S'name_ring_bond'
p88
S'make_rpath'
p89
tp90
sg18
(dp91
g82
(NNNtp92
sg88
(NNNtp93
sg81
(NNNtp94
sg85
(NNNtp95
sg87
(NNNtp96
sg83
(NNNtp97
sg84
(NNNtp98
sg86
(NNNtp99
sg78
(NNNtp100
sg79
(NNNtp101
sg80
(NNNtp102
sg89
(I1
S'1'
Ntp103
ssg21
g22
sg23
S''
sg24
g71
sg25
g2
sg26
g27
sg28
(dp104
sg30
g74
sg31
g75
sg32
S''
sbsS'name_atom'
p105
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p106
(dp107
g4
((dp108
(dp109
tp110
sg13
(dp111
sg15
(S'atom'
p112
tp113
sg18
(dp114
g112
(NNNtp115
ssg21
g22
sg23
S''
sg24
g105
sg25
g2
sg26
g27
sg28
(dp116
sg30
g108
sg31
g109
sg32
S''
sbsS'name_bond'
p117
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p118
(dp119
g4
((dp120
(dp121
tp122
sg13
(dp123
sg15
(S'bond'
p124
tp125
sg18
(dp126
g124
(NNNtp127
ssg21
g22
sg23
S''
sg24
g117
sg25
g2
sg26
g27
sg28
(dp128
sg30
g120
sg31
g121
sg32
S''
sbsS'name_ring_atom'
p129
(ihappydoclib.parseinfo.functioninfo
FunctionInfo
p130
(dp131
g4
((dp132
(dp133
tp134
sg13
(dp135
sg15
(S'atom'
p136
tp137
sg18
(dp138
g136
(NNNtp139
ssg21
g22
sg23
S''
sg24
g129
sg25
g2
sg26
g27
sg28
(dp140
sg30
g132
sg31
g133
sg32
S''
sbstp141
sS'_import_info'
p142
(ihappydoclib.parseinfo.imports
ImportInfo
p143
(dp144
S'_named_imports'
p145
(dp146
sS'_straight_imports'
p147
(lp148
sbsg21
g22
sg23
S'"""Linear Paths\n\nGenerate linear paths for a molecule.\n\nFor example, generate all linear paths up to depth 5\npaths = generatePaths(molecule, maxdepth=5)\n\nThese paths can be used for a variety of cases, but we are using\nthem for the purposes of fingerprinting molecules.\nSee Fingerprint.py\n\nSmarts molecules are a little different than standard\nmolecule representations.  Sometimes, Smarts molecules\nhave complicated logical expressions that are used to\ngenerate the smarts atoms and bonds.\n\nThe following code is pretty much a hack.  It tries\nto determine whether smarts atoms and bonds rely\non any logical operators.\n\nIt an atom or bond does not have this property, then\na path is not made\n\n"""'
p149
sg24
S'SmartsLinearPaths'
p150
sg25
Nsg26
g27
sg28
(dp151
S'include_comments'
p152
I1
sS'cacheFilePrefix'
p153
S'.happydoc.'
p154
sS'useCache'
p155
I1
sS'docStringFormat'
p156
S'StructuredText'
p157
ssg30
g5
sg31
g6
sg32
S''
sbt.