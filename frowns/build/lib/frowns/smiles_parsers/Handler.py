import UserList, sys

class TokenHandler:
    def begin(self):
        """Called before any other callbacks"""
        pass
    def add_token(self, name, pos, text):
        """Called for each token; state name, character position, token text"""
        pass
    def error(self, errmsg, pos, text):
        """Called when an parse error occurs; error position, unparsed text"""
        raise "%s at position %d: %s" % (errmsg, pos, repr(text))
    def end(self):
        pass

class WriteHandler(TokenHandler):
    def __init__(self, outfile = sys.stdout):
        self.outfile = outfile
    def add_token(self, name, pos, text):
        self.outfile.write(" %d --> %s %s\n" % (pos, name, text))

class SaveTokens(TokenHandler, UserList.UserList):
    def begin(self):
        self.data[:] = []
    def add_token(self, name, pos, text):
        self.data.append( (name, pos, text) )

# A mixin class
class SilentErrors:
    def error(self, errmsg, pos, text):
        pass
