#include "colors.inc"
 background {color White} 
 camera { location <0,0,-20> look_at <0,0,0>}
 light_source { <0,0,-20> color White}
// the following declares set global parameters 
#declare AtomRadius = 1.0 ; 
#declare BondRadius = 0.2; 
#declare CarbonColor = Black; 
#declare OxygenColor = Red; 
#declare NitroColor = Blue; 
#declare HydrogenColor = White; 
#declare SulfurColor = Yellow; 
#declare PhosphoColor = Green; 
#declare UnkColor = Pink; 
#declare CarbonRadius = 2.4; 
#declare OxygenRadius = 2.3; 
#declare NitroRadius = 2.35; 
#declare HydrogenRadius = 1.8; 
#declare SulfurRadius = 3; 
#declare PhosphoRadius = 3; 
#declare UnkRadius = 2.4; 
 union {
// 0 one.n1
 sphere { <2.344861,4.629998,-4.464782>, AtomRadius*NitroRadius texture { pigment {color NitroColor} finish { phong 1 } } }
 cylinder{ <2.344861,4.629998,-4.464782> , <2.334625,5.008974,-4.142460>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <2.344861,4.629998,-4.464782> , <2.336834,4.731190,-4.951095>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <2.344861,4.629998,-4.464782> , <2.338577,3.980053,-4.235807>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
// 1 one.h1a
 sphere { <2.324389,5.387950,-3.820138>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.324389,5.387950,-3.820138> , <2.334625,5.008974,-4.142460>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 2 one.h1b
 sphere { <2.328808,4.832383,-5.437408>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.328808,4.832383,-5.437408> , <2.336834,4.731190,-4.951095>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 3 one.c2
 sphere { <2.332294,3.330107,-4.006832>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.332294,3.330107,-4.006832> , <2.338577,3.980053,-4.235807>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.332294,3.330107,-4.006832> , <2.311799,2.791461,-4.455437>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.332294,3.330107,-4.006832> , <2.352220,3.193667,-3.278372>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 4 one.c3
 sphere { <2.291303,2.252816,-4.904042>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.291303,2.252816,-4.904042> , <2.311799,2.791461,-4.455437>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.291303,2.252816,-4.904042> , <2.284677,2.342708,-5.439901>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.291303,2.252816,-4.904042> , <2.272136,1.557014,-4.652644>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 5 one.h3a
 sphere { <2.278050,2.432600,-5.975762>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.278050,2.432600,-5.975762> , <2.284677,2.342708,-5.439901>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 6 one.c4
 sphere { <2.252969,0.861212,-4.401245>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.252969,0.861212,-4.401245> , <2.272136,1.557014,-4.652644>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.252969,0.861212,-4.401245> , <2.229519,0.443619,-4.749037>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.252969,0.861212,-4.401245> , <2.258801,0.736058,-3.714194>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 7 one.h4a
 sphere { <2.206070,0.026027,-5.096830>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.206070,0.026027,-5.096830> , <2.229519,0.443619,-4.749037>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 8 one.c5
 sphere { <2.264633,0.610905,-3.027143>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.264633,0.610905,-3.027143> , <2.258801,0.736058,-3.714194>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.264633,0.610905,-3.027143> , <2.314652,1.176209,-2.553520>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.264633,0.610905,-3.027143> , <2.198855,-0.229659,-2.721703>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 9 one.c6
 sphere { <2.364671,1.741513,-2.079896>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.364671,1.741513,-2.079896> , <2.314652,1.176209,-2.553520>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.364671,1.741513,-2.079896> , <2.407497,1.644502,-1.547584>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.364671,1.741513,-2.079896> , <2.368408,2.399370,-2.314904>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 10 one.h6a
 sphere { <2.450323,1.547491,-1.015273>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.450323,1.547491,-1.015273> , <2.407497,1.644502,-1.547584>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 11 one.c7
 sphere { <2.372144,3.057228,-2.549911>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <2.372144,3.057228,-2.549911> , <2.368408,2.399370,-2.314904>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.372144,3.057228,-2.549911> , <2.352220,3.193667,-3.278372>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <2.372144,3.057228,-2.549911> , <2.395429,3.468917,-2.194811>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 12 one.h7a
 sphere { <2.418713,3.880606,-1.839711>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.418713,3.880606,-1.839711> , <2.395429,3.468917,-2.194811>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 13 one.s8
 sphere { <2.133078,-1.070223,-2.416262>, AtomRadius*SulfurRadius texture { pigment {color SulfurColor} finish { phong 1 } } }
 cylinder{ <2.133078,-1.070223,-2.416262> , <2.198855,-0.229659,-2.721703>, BondRadius texture { pigment { color SulfurColor } finish { phong 1 } } }
 cylinder{ <2.133078,-1.070223,-2.416262> , <2.174139,-1.581211,-3.027258>, BondRadius texture { pigment { color SulfurColor } finish { phong 1 } } }
 cylinder{ <2.133078,-1.070223,-2.416262> , <2.748620,-1.228513,-1.937687>, BondRadius texture { pigment { color SulfurColor } finish { phong 1 } } }
 cylinder{ <2.133078,-1.070223,-2.416262> , <1.377066,-1.197138,-1.986336>, BondRadius texture { pigment { color SulfurColor } finish { phong 1 } } }
// 14 one.o9
 sphere { <2.215201,-2.092200,-3.638253>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <2.215201,-2.092200,-3.638253> , <2.174139,-1.581211,-3.027258>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 15 one.o10
 sphere { <3.364162,-1.386804,-1.459113>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <3.364162,-1.386804,-1.459113> , <2.748620,-1.228513,-1.937687>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 16 one.n11
 sphere { <0.621054,-1.324053,-1.556411>, AtomRadius*NitroRadius texture { pigment {color NitroColor} finish { phong 1 } } }
 cylinder{ <0.621054,-1.324053,-1.556411> , <1.377066,-1.197138,-1.986336>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <0.621054,-1.324053,-1.556411> , <0.250563,-1.938017,-1.744483>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <0.621054,-1.324053,-1.556411> , <0.604704,-1.117448,-0.840236>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
// 17 one.c12
 sphere { <-0.119929,-2.551980,-1.932556>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.119929,-2.551980,-1.932556> , <0.250563,-1.938017,-1.744483>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.119929,-2.551980,-1.932556> , <0.246185,-2.964736,-2.017710>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.119929,-2.551980,-1.932556> , <-0.465522,-2.709621,-1.524624>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.119929,-2.551980,-1.932556> , <-0.530034,-2.409934,-2.566806>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 18 one.h12a
 sphere { <0.612299,-3.377491,-2.102864>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <0.612299,-3.377491,-2.102864> , <0.246185,-2.964736,-2.017710>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 19 one.h12b
 sphere { <-0.811114,-2.867261,-1.116693>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.811114,-2.867261,-1.116693> , <-0.465522,-2.709621,-1.524624>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 20 one.c13
 sphere { <-0.940138,-2.267888,-3.201057>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.940138,-2.267888,-3.201057> , <-0.530034,-2.409934,-2.566806>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.940138,-2.267888,-3.201057> , <-0.578712,-2.125443,-3.603395>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.940138,-2.267888,-3.201057> , <-1.316011,-2.896895,-3.426568>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.940138,-2.267888,-3.201057> , <-1.416457,-1.673591,-3.093296>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 21 one.h13a
 sphere { <-0.217285,-1.982998,-4.005733>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.217285,-1.982998,-4.005733> , <-0.578712,-2.125443,-3.603395>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 22 one.c14
 sphere { <-1.691885,-3.525902,-3.652079>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.691885,-3.525902,-3.652079> , <-1.316011,-2.896895,-3.426568>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.691885,-3.525902,-3.652079> , <-2.051224,-3.693917,-3.260720>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.691885,-3.525902,-3.652079> , <-1.970099,-3.428048,-4.124322>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.691885,-3.525902,-3.652079> , <-1.328120,-3.937377,-3.748951>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 23 one.h14a
 sphere { <-2.410563,-3.861933,-2.869362>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.410563,-3.861933,-2.869362> , <-2.051224,-3.693917,-3.260720>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 24 one.h14b
 sphere { <-2.248314,-3.330193,-4.596564>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.248314,-3.330193,-4.596564> , <-1.970099,-3.428048,-4.124322>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 25 one.h14c
 sphere { <-0.964355,-4.348852,-3.845822>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.964355,-4.348852,-3.845822> , <-1.328120,-3.937377,-3.748951>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 26 one.c15
 sphere { <-1.892776,-1.079294,-2.985535>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.892776,-1.079294,-2.985535> , <-1.416457,-1.673591,-3.093296>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.892776,-1.079294,-2.985535> , <-1.592834,-0.620020,-2.882691>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.892776,-1.079294,-2.985535> , <-2.226645,-1.168607,-2.548693>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.892776,-1.079294,-2.985535> , <-2.200622,-0.983308,-3.438774>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 27 one.h15a
 sphere { <-1.292892,-0.160745,-2.779848>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-1.292892,-0.160745,-2.779848> , <-1.592834,-0.620020,-2.882691>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 28 one.h15b
 sphere { <-2.560514,-1.257920,-2.111851>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.560514,-1.257920,-2.111851> , <-2.226645,-1.168607,-2.548693>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 29 one.h15c
 sphere { <-2.508468,-0.887322,-3.892014>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.508468,-0.887322,-3.892014> , <-2.200622,-0.983308,-3.438774>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 30 one.c16
 sphere { <0.588354,-0.910843,-0.124062>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <0.588354,-0.910843,-0.124062> , <0.604704,-1.117448,-0.840236>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.588354,-0.910843,-0.124062> , <0.852542,-0.427625,-0.078578>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.588354,-0.910843,-0.124062> , <0.058406,-0.830028,0.006148>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.588354,-0.910843,-0.124062> , <0.890450,-1.376909,0.403310>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 31 one.h16a
 sphere { <1.116731,0.055592,-0.033093>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <1.116731,0.055592,-0.033093> , <0.852542,-0.427625,-0.078578>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 32 one.h16b
 sphere { <-0.471541,-0.749214,0.136359>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.471541,-0.749214,0.136359> , <0.058406,-0.830028,0.006148>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 33 one.c17
 sphere { <1.192545,-1.842974,0.930682>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <1.192545,-1.842974,0.930682> , <0.890450,-1.376909,0.403310>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.192545,-1.842974,0.930682> , <1.699347,-1.973303,0.736251>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.192545,-1.842974,0.930682> , <0.824870,-2.464361,0.961079>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.192545,-1.842974,0.930682> , <1.285985,-1.509910,1.637724>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 34 one.h17a
 sphere { <2.206150,-2.103632,0.541821>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.206150,-2.103632,0.541821> , <1.699347,-1.973303,0.736251>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 35 one.o18
 sphere { <0.457195,-3.085747,0.991476>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <0.457195,-3.085747,0.991476> , <0.824870,-2.464361,0.961079>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
 cylinder{ <0.457195,-3.085747,0.991476> , <0.628886,-3.332416,0.569816>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 36 one.h18a
 sphere { <0.800577,-3.579084,0.148157>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <0.800577,-3.579084,0.148157> , <0.628886,-3.332416,0.569816>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 37 one.c19
 sphere { <1.379426,-1.176845,2.344766>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <1.379426,-1.176845,2.344766> , <1.285985,-1.509910,1.637724>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.379426,-1.176845,2.344766> , <1.925954,-1.053143,2.328051>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.379426,-1.176845,2.344766> , <1.076446,-0.540885,2.490308>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.379426,-1.176845,2.344766> , <1.344392,-1.663585,2.946055>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 38 one.h19a
 sphere { <2.472482,-0.929440,2.311335>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <2.472482,-0.929440,2.311335> , <1.925954,-1.053143,2.328051>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 39 one.n20
 sphere { <0.773465,0.095076,2.635850>, AtomRadius*NitroRadius texture { pigment {color NitroColor} finish { phong 1 } } }
 cylinder{ <0.773465,0.095076,2.635850> , <1.076446,-0.540885,2.490308>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <0.773465,0.095076,2.635850> , <1.058985,0.432495,2.860378>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
 cylinder{ <0.773465,0.095076,2.635850> , <0.121814,0.251247,2.652008>, BondRadius texture { pigment { color NitroColor } finish { phong 1 } } }
// 40 one.h20a
 sphere { <1.344504,0.769914,3.084906>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <1.344504,0.769914,3.084906> , <1.058985,0.432495,2.860378>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 41 one.c21
 sphere { <-0.529838,0.407418,2.668167>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.529838,0.407418,2.668167> , <0.121814,0.251247,2.652008>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.529838,0.407418,2.668167> , <-0.978369,0.006837,2.454640>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.529838,0.407418,2.668167> , <-0.732770,0.999488,3.000515>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 42 one.o22
 sphere { <-1.426901,-0.393743,2.241113>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <-1.426901,-0.393743,2.241113> , <-0.978369,0.006837,2.454640>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 43 one.o23
 sphere { <-0.935702,1.591559,3.332864>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <-0.935702,1.591559,3.332864> , <-0.732770,0.999488,3.000515>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
 cylinder{ <-0.935702,1.591559,3.332864> , <-0.921332,2.220509,2.960934>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 44 one.c24
 sphere { <-0.906963,2.849459,2.589003>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.906963,2.849459,2.589003> , <-0.921332,2.220509,2.960934>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.906963,2.849459,2.589003> , <-1.117138,3.228426,2.941233>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.906963,2.849459,2.589003> , <-0.229639,3.125816,2.383569>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.906963,2.849459,2.589003> , <-1.318489,2.914559,1.950745>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 45 one.h24a
 sphere { <-1.327312,3.607393,3.293462>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-1.327312,3.607393,3.293462> , <-1.117138,3.228426,2.941233>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 46 one.c25
 sphere { <0.447685,3.402172,2.178135>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <0.447685,3.402172,2.178135> , <-0.229639,3.125816,2.383569>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.447685,3.402172,2.178135> , <0.768462,3.008577,1.953310>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.447685,3.402172,2.178135> , <0.712173,3.618633,2.617331>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <0.447685,3.402172,2.178135> , <0.322916,3.912909,1.677412>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 47 one.h25a
 sphere { <1.089239,2.614983,1.728485>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <1.089239,2.614983,1.728485> , <0.768462,3.008577,1.953310>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 48 one.h25b
 sphere { <0.976661,3.835094,3.056527>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <0.976661,3.835094,3.056527> , <0.712173,3.618633,2.617331>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 49 one.o26
 sphere { <0.198146,4.423647,1.176690>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <0.198146,4.423647,1.176690> , <0.322916,3.912909,1.677412>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
 cylinder{ <0.198146,4.423647,1.176690> , <-0.486601,4.344824,0.951333>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 50 one.c27
 sphere { <-1.171349,4.266001,0.725977>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.171349,4.266001,0.725977> , <-0.486601,4.344824,0.951333>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.171349,4.266001,0.725977> , <-1.470767,4.707217,0.891261>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.171349,4.266001,0.725977> , <-1.179174,4.224661,0.004707>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.171349,4.266001,0.725977> , <-1.450682,3.622830,1.019231>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 51 one.h27a
 sphere { <-1.770185,5.148434,1.056544>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-1.770185,5.148434,1.056544> , <-1.470767,4.707217,0.891261>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 52 one.o28
 sphere { <-1.187000,4.183321,-0.716563>, AtomRadius*OxygenRadius texture { pigment {color OxygenColor} finish { phong 1 } } }
 cylinder{ <-1.187000,4.183321,-0.716563> , <-1.179174,4.224661,0.004707>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
 cylinder{ <-1.187000,4.183321,-0.716563> , <-1.404369,3.517773,-0.888241>, BondRadius texture { pigment { color OxygenColor } finish { phong 1 } } }
// 53 one.c29
 sphere { <-1.621737,2.852224,-1.059919>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.621737,2.852224,-1.059919> , <-1.404369,3.517773,-0.888241>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.621737,2.852224,-1.059919> , <-2.152905,2.870776,-1.226560>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.621737,2.852224,-1.059919> , <-1.306945,2.654393,-1.475416>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.621737,2.852224,-1.059919> , <-1.548267,2.426618,-0.438686>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 54 one.h29a
 sphere { <-2.684072,2.889327,-1.393200>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.684072,2.889327,-1.393200> , <-2.152905,2.870776,-1.226560>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 55 one.h29b
 sphere { <-0.992152,2.456560,-1.890913>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.992152,2.456560,-1.890913> , <-1.306945,2.654393,-1.475416>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 56 one.c30
 sphere { <-1.474796,2.001012,0.182547>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.474796,2.001012,0.182547> , <-1.548267,2.426618,-0.438686>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.474796,2.001012,0.182547> , <-0.949762,1.826438,0.203391>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.474796,2.001012,0.182547> , <-1.820415,1.565285,0.197973>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.474796,2.001012,0.182547> , <-1.602406,2.490335,0.747517>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 57 one.h30a
 sphere { <-0.424727,1.651865,0.224236>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.424727,1.651865,0.224236> , <-0.949762,1.826438,0.203391>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 58 one.h30b
 sphere { <-2.166035,1.129558,0.213399>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.166035,1.129558,0.213399> , <-1.820415,1.565285,0.197973>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 59 one.c31
 sphere { <-1.730015,2.979659,1.312486>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-1.730015,2.979659,1.312486> , <-1.318489,2.914559,1.950745>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.730015,2.979659,1.312486> , <-1.450682,3.622830,1.019231>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.730015,2.979659,1.312486> , <-1.602406,2.490335,0.747517>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-1.730015,2.979659,1.312486> , <-2.274561,3.020916,1.420450>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 60 one.h31a
 sphere { <-2.819108,3.062173,1.528413>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.819108,3.062173,1.528413> , <-2.274561,3.020916,1.420450>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 61 one.c32
 sphere { <1.309359,-2.150325,3.547345>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <1.309359,-2.150325,3.547345> , <1.344392,-1.663585,2.946055>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.309359,-2.150325,3.547345> , <1.590470,-2.611348,3.410603>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.309359,-2.150325,3.547345> , <1.581031,-1.921264,3.975817>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <1.309359,-2.150325,3.547345> , <0.615966,-2.336296,3.785311>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 62 one.h32a
 sphere { <1.871581,-3.072372,3.273861>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <1.871581,-3.072372,3.273861> , <1.590470,-2.611348,3.410603>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 63 one.h32b
 sphere { <1.852703,-1.692204,4.404288>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <1.852703,-1.692204,4.404288> , <1.581031,-1.921264,3.975817>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 64 one.c33
 sphere { <-0.673140,-3.736975,3.647542>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.673140,-3.736975,3.647542> , <-0.375283,-3.129622,3.835409>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.673140,-3.736975,3.647542> , <-0.407925,-4.083767,3.322839>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.673140,-3.736975,3.647542> , <-1.357721,-3.905675,3.871102>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 65 one.h33a
 sphere { <-0.142710,-4.430558,2.998134>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.142710,-4.430558,2.998134> , <-0.407925,-4.083767,3.322839>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 66 one.c34
 sphere { <-2.042302,-4.074375,4.094662>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-2.042302,-4.074375,4.094662> , <-1.357721,-3.905675,3.871102>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.042302,-4.074375,4.094662> , <-2.271523,-4.540972,3.933864>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.042302,-4.074375,4.094662> , <-2.396128,-3.631912,4.503674>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 67 one.h34a
 sphere { <-2.500744,-5.007570,3.773067>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.500744,-5.007570,3.773067> , <-2.271523,-4.540972,3.933864>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 68 one.c35
 sphere { <-2.749953,-3.189449,4.912687>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-2.749953,-3.189449,4.912687> , <-2.396128,-3.631912,4.503674>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.749953,-3.189449,4.912687> , <-3.256869,-3.308480,5.071164>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.749953,-3.189449,4.912687> , <-2.436052,-2.555931,5.127369>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 69 one.h35a
 sphere { <-3.763785,-3.427510,5.229642>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-3.763785,-3.427510,5.229642> , <-3.256869,-3.308480,5.071164>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 70 one.c36
 sphere { <-2.122151,-1.922412,5.342050>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-2.122151,-1.922412,5.342050> , <-2.436052,-2.555931,5.127369>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.122151,-1.922412,5.342050> , <-2.398582,-1.572205,5.653739>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-2.122151,-1.922412,5.342050> , <-1.473518,-1.766124,5.134718>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 71 one.h36a
 sphere { <-2.675014,-1.221998,5.965427>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-2.675014,-1.221998,5.965427> , <-2.398582,-1.572205,5.653739>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 72 one.c37
 sphere { <-0.824885,-1.609834,4.927385>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.824885,-1.609834,4.927385> , <-1.473518,-1.766124,5.134718>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.824885,-1.609834,4.927385> , <-0.603433,-1.135337,5.076685>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.824885,-1.609834,4.927385> , <-0.451157,-2.066051,4.475330>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
// 73 one.h37a
 sphere { <-0.381981,-0.660840,5.225984>, AtomRadius*HydrogenRadius texture { pigment {color HydrogenColor} finish { phong 1 } } }
 cylinder{ <-0.381981,-0.660840,5.225984> , <-0.603433,-1.135337,5.076685>, BondRadius texture { pigment { color HydrogenColor } finish { phong 1 } } }
// 74 one.c38
 sphere { <-0.077427,-2.522268,4.023276>, AtomRadius*CarbonRadius texture { pigment {color CarbonColor} finish { phong 1 } } }
 cylinder{ <-0.077427,-2.522268,4.023276> , <-0.375283,-3.129622,3.835409>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.077427,-2.522268,4.023276> , <-0.451157,-2.066051,4.475330>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
 cylinder{ <-0.077427,-2.522268,4.023276> , <0.615966,-2.336296,3.785311>, BondRadius texture { pigment { color CarbonColor } finish { phong 1 } } }
  matrix < -0.276498, -0.308172, 0.910263, -0.735246, -0.542101, -0.406866, 0.618839, -0.781765, -0.076693, 0,0,0 >
 rotate <180,0,0>
}
