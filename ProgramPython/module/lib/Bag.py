#filename = Bag.py

class Bag:
    def __init__(self):
        self.data = []
    def add(self, x):
        self.data.append(x)
    def addtwice(self, x):
        self.add(x)
        self.add(x)
    def show(self):
        print self.data
        print "Show Bag"

