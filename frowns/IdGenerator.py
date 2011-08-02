"""IdGenerator
A simple class to provide a sequence of integers that can be used
for identification purposes

generator = IdGenerator()
generator() -> returns 0
generator() -> returns 1
...

generator = IdGenerator(1000)
generator() -> returns 1000
generator() -> returns 1001
"""

class IdGenerator:
    def __init__(self, start=1):
        self.start = start-1

    def __call__(self):
        self.start += 1
        return self.start

defaultGenerator = IdGenerator()

def test():
    generator = IdGenerator()
    assert generator() == 1
    assert generator() == 2

    generator = IdGenerator(start=1000)
    assert generator() == 1000
    assert generator() == 1001

if __name__ == "__main__":
    test()


