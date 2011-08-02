"""Turn incomplete SMILES into something more likely to be valid.

Daylight does not accept SMILES which are not well formed.  Suppose
you are typing in a new SMILES string, like "CC(=O)[0-]".  Part way
through you have terms like "CC(=" or "CC(=O)[O".  Daylight doesn't
like these, so you can't, for example, use the depiction code to see
what the partial string looks like.

This module tries to fix up the string to something Daylight does
accept.  In the above case, the fixed forms of the partial terms are
"CC" and "CC(=O)[O]".

Incomplete ring closures are turned into "*" atoms with atomic number
equal to the ring closure number.  For example, "C5CC" becomes
"C[5*]CC".


Bugs:

Daylight doesn't accept some aromatic notations like "ccc" but will
accept others like "cccc".  If that's a problem, as a workaround you
might try passing the string to the toolkit and if it fails convert
all the aromatic terms to uppercase and try again.

"""

import Smiles, Handler


def cleanup_closure(s):
    if s[:1] == "%":
        return int(s[1:])
    return int(s)

class SilentSaveTokens(Handler.SilentErrors, Handler.SaveTokens):
    pass

def cleanup(s):
    save_h = SilentSaveTokens()
    Smiles.tokenize(s, save_h)
    save = []
    for name, pos, text in save_h:
        save.append( (name, text) )

    closures = {}

    in_quote = 0
    i = 0
    for name, text in save:
        if name == "open_bracket":
            in_quote = 1
            has_element = 0
            pos = i
        elif name == "close_bracket":
            in_quote = 0
        elif name == "element" and in_quote:
            has_element = 1
        i = i + 1
    if in_quote:
        if has_element:
            save.append( ("close_bracket", "]") )
        else:
            del save[pos:]
    
    i = 0
    for name, text in save:
        if name == "closure":
            val = cleanup_closure(text)
            if closures.has_key(val):
                del closures[val]
            else:
                closures[val] = i
        i = i + 1

    for val, i in closures.items():
        save[i] = ("fake", "[%d*]" % val)

    # trim trailing '(', '.' and bond symbols
    while save and (save[-1][0] == "close_branch" or
                    save[-1][0] == "dot" or
                    save[-1][0] == "bond"):
        save.pop()


    paren_count = 0
    for name, text in save:
        if name == "open_branch":
            paren_count = paren_count + 1
        elif name == "close_branch":
            paren_count = paren_count - 1
        assert paren_count >= 0, paren_count
    save.extend( [("close_branch", ")")] * paren_count )

    smi = ""
    for name, text in save:
        smi = smi + text
    return smi

def test_show(s):
    for i in range(1, len(s)+1):
        print repr(s[:i]), repr(cleanup(s[:i]))

def test():
    test_show("[Xe].[6Li]")
    test_show("c1cccc1Oc%34ccccc%34")
    test_show("C(CC)[NH4+]")
    test_show("C(CC(=O)[O-2]")

if __name__ == "__main__":
    test()

