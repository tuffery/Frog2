"""A set of objects that can create logical expressions via
combining logical operators and Smarts Primitives.
"""
class NotMatch:
    def __init__(self, child):
        self.child = child
    def __eq__(self, obj):
        return not self.child == obj
    def __str__(self):
        return "not (%s)" % (self.child,)

class AndMatch:
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def __eq__(self, obj):
        return self.left == obj and self.right == obj
    def __str__(self):
        return "AND(%s, %s)" % (self.left, self.right)
    
class OrMatch:
    def __init__(self, left, right):
        self.left = left
        self.right = right
    def __eq__(self, obj):
        return (self.left == obj) or (self.right == obj)
    def __str__(self):
        return "OR(%s, %s)" % (self.left, self.right)

bool_unary_not = 76
bool_strong_and = 77
bool_or = 78
bool_weak_and = 79
binary_operators = [bool_strong_and, bool_or, bool_weak_and]
boolean_operators = binary_operators + [bool_unary_not]
text_to_bool = {
    "&": bool_strong_and,
    ",": bool_or,
    ";": bool_weak_and,
    "!": bool_unary_not,
    }

class ExpressionList:
    def __init__(self):
        self.matchers = []
        
    def __nonzero__(self):
        return len(self.matchers) != 0
    
    def add_matcher(self, obj):
        # see if the matcher is not in the boolean list
        if self.matchers:
            target = self.matchers[-1]
            for _matcher in boolean_operators:
                if target is _matcher:
                    break
            else:
                self.matchers.append(bool_strong_and)
        #if self.matchers and self.matchers[-1] not in boolean_operators:
        #    self.matchers.append(bool_strong_and)
        self.matchers.append(obj)
        
    def add_operator(self, op):
        assert op in binary_operators or op == bool_unary_not
        if __debug__:
            if self.matchers:
                if op in binary_operators:
                    s = self.matchers[-1]
                    for _matcher in binary_operators:
                        assert s is not _matcher
            else:
                for _matcher in binary_operators:
                    assert op is not _matcher

        self.matchers.append(op)

    def make_matcher(self):
        matchers = self.matchers[:]
        i = 0
        while i < len(matchers):
            if matchers[i] is bool_unary_not:
                matchers[i:i+2] = [NotMatch(matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] is bool_strong_and:
                matchers[i-1:i+2] = [AndMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] is bool_or:
                matchers[i-1:i+2] = [OrMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        i = 1
        while i < len(matchers):
            if matchers[i] is bool_weak_and:
                matchers[i-1:i+2] = [AndMatch(matchers[i-1], matchers[i+1])]
            else:
                i = i + 1
        assert len(matchers) == 1, matchers
        return matchers[0]

class AtomExpression(ExpressionList):
    pass

class BondExpression(ExpressionList):
    pass
        

