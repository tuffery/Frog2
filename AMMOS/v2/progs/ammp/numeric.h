/* numeric.h
*
*  definitions of constants as variables
* this can speed up code.
*/
/* just the extern definition for routines */
#ifdef NUMERIC_DEFINED

extern float one, two,three,four,five,six,seven,eight,nine,twelve;
extern float third,eightth,fourteenth;

#else
/* do the explicit definitions */
#ifdef DEFINE_NUMERIC
#define NUMERIC_DEFINED
float one = 1.;
float two = 2.;
float three = 3.;
float four = 4.;
float five = 5.;
float six = 6.;
float seven = 7.;
float eight = 8.;
float nine = 9.;
float twelve = 12.;
float half = 0.5 ;
float third = 1./3;
float eightth = .125;
float fourteenth = 1./14;


#else
/*  do the explicit variable substitution */
#define one 1.
#define two 2.
#define three 3.
#define four 4.
#define five 5.
#define six 6.
#define seven 7.
#define eight 8.
#define nine 9.
#define twelve 12. 
#define half  0.5
#define third 1./3
#define eightth 1./8
#define fourteenth 1./14
#endif
#endif

#define ROOTHALF 0.7071067811865475244008444
#define ROOT2    1.4142135623730950488016887
#define PI 3.1415926535897932384626433
#define TWOPI 6.2831853071795864769252866
#define ROOTPI 1.7724538509055160272981674
/* BOHR is the BOHR radius */
#define BOHR 0.5291771
#define INVBOHR 1.8897265
/*
#define INVBOHR 1.8897265*1.8897265
*/
/* HARTREE is 332.17752/(2*BOHR) */
/*#define HARTREE 627.725 */
#define HARTREE 313.8625



