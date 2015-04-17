#Filename = dictionary.py
#methods that implement some hash function in perl

'''
Hash in perl can use like $hash{key1}{key2}, but in python such function is not functional.
This method dict_add can create dictionary like perl hash.
breakpoint = {}
breakpoint = dict_add(breakpoint, rils, chro, newbp)
print breakpoint[rils][chro] will give newbp
'''
def dict_add( dict0, key1, key2, value ):
   if dict0.has_key(key1):
       temp  = dict0.copy()
       dict1 = temp[key1]
       dict1.update({key2:value})
       temp[key1] = dict1
   else:
       temp  = dict0.copy()
       dict1 = {key2:value}
       temp.update({key1:{key2:value}})

   return temp


