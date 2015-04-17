from lib.Bag import Bag
from lib.Employee import Employee

newbag = Bag()
newbag.add('nike')
newbag.addtwice('adidas')
newbag.show()

emp1 = Employee("John", 2000)
emp2 = Employee("Jack", 3000)

emp1.displayEmployee()
emp1.displayCount()
emp2.displayEmployee()
emp2.displayCount()

